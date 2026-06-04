#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <string>
#include <algorithm>

struct BodyState {
    double x, y, vx, vy;
};

struct SystemState {
    BodyState b1, b2, b3;

    SystemState operator+(const SystemState& other) const {
        return {
            {b1.x + other.b1.x, b1.y + other.b1.y, b1.vx + other.b1.vx, b1.vy + other.b1.vy},
            {b2.x + other.b2.x, b2.y + other.b2.y, b2.vx + other.b2.vx, b2.vy + other.b2.vy},
            {b3.x + other.b3.x, b3.y + other.b3.y, b3.vx + other.b3.vx, b3.vy + other.b3.vy}
        };
    }
    SystemState operator-(const SystemState& other) const {
        return {
            {b1.x - other.b1.x, b1.y - other.b1.y, b1.vx - other.b1.vx, b1.vy - other.b1.vy},
            {b2.x - other.b2.x, b2.y - other.b2.y, b2.vx - other.b2.vx, b2.vy - other.b2.vy},
            {b3.x - other.b3.x, b3.y - other.b3.y, b3.vx - other.b3.vx, b3.vy - other.b3.vy}
        };
    }
    SystemState operator*(double scalar) const {
        return {
            {b1.x * scalar, b1.y * scalar, b1.vx * scalar, b1.vy * scalar},
            {b2.x * scalar, b2.y * scalar, b2.vx * scalar, b2.vy * scalar},
            {b3.x * scalar, b3.y * scalar, b3.vx * scalar, b3.vy * scalar}
        };
    }
};

void add_gravity(const BodyState& A, const BodyState& B, double& ax, double& ay) {
    double dx = B.x - A.x;
    double dy = B.y - A.y;
    double r2 = dx * dx + dy * dy;
    double r3 = r2 * std::sqrt(r2);
    if (r3 < 1e-6) r3 = 1e-6;
    ax += dx / r3;
    ay += dy / r3;
}

SystemState derivatives(double t, const SystemState& s) {
    SystemState d;
    
    d.b1.x = s.b1.vx; d.b1.y = s.b1.vy;
    d.b2.x = s.b2.vx; d.b2.y = s.b2.vy;
    d.b3.x = s.b3.vx; d.b3.y = s.b3.vy;

    d.b1.vx = 0; d.b1.vy = 0;
    add_gravity(s.b1, s.b2, d.b1.vx, d.b1.vy);
    add_gravity(s.b1, s.b3, d.b1.vx, d.b1.vy);

    d.b2.vx = 0; d.b2.vy = 0;
    add_gravity(s.b2, s.b1, d.b2.vx, d.b2.vy);
    add_gravity(s.b2, s.b3, d.b2.vx, d.b2.vy);

    d.b3.vx = 0; d.b3.vy = 0;
    add_gravity(s.b3, s.b1, d.b3.vx, d.b3.vy);
    add_gravity(s.b3, s.b2, d.b3.vx, d.b3.vy);

    return d;
}

// 1. Метод RK4 (Классический Рунге-Кутта 4-го порядка).
SystemState stepRK4(double t, const SystemState& s, double h) {
    SystemState k1 = derivatives(t, s);
    SystemState k2 = derivatives(t + h/2.0, s + k1 * (h/2.0));
    SystemState k3 = derivatives(t + h/2.0, s + k2 * (h/2.0));
    SystemState k4 = derivatives(t + h, s + k3 * h);
    return s + (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (h / 6.0);
}

// 2. Метод RK3/8 (Правило 3/8 Рунге-Кутта 4-го порядка)
SystemState stepRK38(double t, const SystemState& s, double h) {
    SystemState k1 = derivatives(t, s);
    SystemState k2 = derivatives(t + h/3.0, s + k1 * (h/3.0));
    SystemState k3 = derivatives(t + 2.0*h/3.0, s + k1 * (-h/3.0) + k2 * h);
    SystemState k4 = derivatives(t + h, s + k1 * h - k2 * h + k3 * h);
    return s + (k1 + k2 * 3.0 + k3 * 3.0 + k4) * (h / 8.0);
}

// 3. Метод Рунге-Кутта 3-го порядка (Метод Кутты)
SystemState stepRK3(double t, const SystemState& s, double h) {
    SystemState k1 = derivatives(t, s);
    SystemState k2 = derivatives(t + h/2.0, s + k1 * (h/2.0));
    SystemState k3 = derivatives(t + h, s + k1 * (-h) + k2 * 2.0 * h);
    return s + (k1 + k2 * 4.0 + k3) * (h / 6.0);
}

// 4. DOPRI5 (Дорман-Принс 5(4)) с адаптивным шагом
SystemState stepDOPRI5(double t, const SystemState& s, double& h, double tolerance) {
    const double c2 = 1.0/5.0, c3 = 3.0/10.0, c4 = 4.0/5.0, c5 = 8.0/9.0, c6 = 1.0;
    while (true) {
        SystemState k1 = derivatives(t, s);
        SystemState k2 = derivatives(t + c2*h, s + k1 * (h * (1.0/5.0)));
        SystemState k3 = derivatives(t + c3*h, s + k1 * (h * (3.0/40.0)) + k2 * (h * (9.0/40.0)));
        SystemState k4 = derivatives(t + c4*h, s + k1 * (h * (44.0/45.0)) - k2 * (h * (56.0/15.0)) + k3 * (h * (32.0/9.0)));
        SystemState k5 = derivatives(t + c5*h, s + k1 * (h * (19372.0/6561.0)) - k2 * (h * (25360.0/2187.0)) + k3 * (h * (64448.0/6561.0)) - k4 * (h * (212.0/729.0)));
        SystemState k6 = derivatives(t + c6*h, s + k1 * (h * (9017.0/3168.0)) - k2 * (h * (355.0/33.0)) + k3 * (h * (46732.0/5247.0)) + k4 * (h * (49.0/176.0)) - k5 * (h * (5103.0/18656.0)));
        
        SystemState s_next = s + (k1 * (35.0/384.0) + k3 * (500.0/1113.0) + k4 * (125.0/192.0) - k5 * (2187.0/6784.0) + k6 * (11.0/84.0)) * h;
        SystemState k7 = derivatives(t + h, s_next);
        
        SystemState error = (k1 * (71.0/57600.0) - k3 * (71.0/16695.0) + k4 * (71.0/1920.0) - k5 * (17253.0/339200.0) + k6 * (22.0/525.0) - k7 * (1.0/40.0)) * h;
        
        double err_norm = std::sqrt(
            error.b1.x*error.b1.x + error.b1.y*error.b1.y + error.b2.x*error.b2.x + error.b2.y*error.b2.y + error.b3.x*error.b3.x + error.b3.y*error.b3.y
        );
        
        if (err_norm <= tolerance) {
            double scale = 0.9 * std::pow(tolerance / (err_norm + 1e-16), 0.2);
            h *= std::max(0.1, std::min(5.0, scale));
            return s_next;
        }
        h *= std::max(0.1, std::min(0.5, 0.9 * std::pow(tolerance / err_norm, 0.25)));
    }
}

int main() {
    std::cout << std::scientific << std::setprecision(17);

    double t_end = 6.3259;
    double h_fixed = 0.005;
    double ref_tolerance = 1e-13;

    std::ifstream config("config.txt");
    if (config.is_open()) {
        config >> t_end >> h_fixed >> ref_tolerance;
        config.close();
    }

    double x3 = 0.97000436, y3 = -0.24308753;
    double vx3 = -0.46620531, vy3 = -0.43236573;

    SystemState init_state = {
        { -x3, -y3, vx3, vy3 },
        { 0.0, 0.0, -2.0 * vx3, -2.0 * vy3 },
        { x3, y3, vx3, vy3 }
    };

    struct ReferencePoint { double t; SystemState state; };
    std::vector<ReferencePoint> ref_path;
    
    double t = 0.0;
    SystemState s_ref = init_state;
    double h_adaptive = h_fixed;
    
    ref_path.push_back({t, s_ref});
    while (t < t_end) {
        if (t + h_adaptive > t_end) h_adaptive = t_end - t;
        s_ref = stepDOPRI5(t, s_ref, h_adaptive, ref_tolerance);
        t += h_adaptive;
        ref_path.push_back({t, s_ref});
    }

    std::cout << "t,b1_x,b1_y,b2_x,b2_y,b3_x,b3_y\n";
    for (const auto& pt : ref_path) {
        std::cout << pt.t << ","
                  << pt.state.b1.x << "," << pt.state.b1.y << ","
                  << pt.state.b2.x << "," << pt.state.b2.y << ","
                  << pt.state.b3.x << "," << pt.state.b3.y << "\n";
    }
    
    return 0;
}
