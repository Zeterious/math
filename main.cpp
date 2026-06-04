#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <string>
#include <algorithm>

struct State {
    double x, y, vx, vy;

    State operator+(const State& other) const {
        return { x + other.x, y + other.y, vx + other.vx, vy + other.vy };
    }
    State operator-(const State& other) const {
        return { x - other.x, y - other.y, vx - other.vx, vy - other.vy };
    }
    State operator*(double scalar) const {
        return { x * scalar, y * scalar, vx * scalar, vy * scalar };
    }
};


State derivatives(double t, const State& s) {
    double r3 = std::pow(s.x * s.x + s.y * s.y, 1.5);
    if (r3 < 1e-6) r3 = 1e-6;

    State d;
    d.x = s.vx;
    d.y = s.vy;
    d.vx = -s.x / r3;
    d.vy = -s.y / r3;
    return d;
}

// 1. Метод RK4 (Классический Рунге-Кутта 4-го порядка)
State stepRK4(double t, const State& s, double h) {
    State k1 = derivatives(t, s);
    State k2 = derivatives(t + h / 2.0, s + k1 * (h / 2.0));
    State k3 = derivatives(t + h / 2.0, s + k2 * (h / 2.0));
    State k4 = derivatives(t + h, s + k3 * h);
    return s + (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (h / 6.0);
}

// 2. Метод RK3/8 (Правило 3/8 Рунге-Кутта 4-го порядка)
State stepRK38(double t, const State& s, double h) {
    State k1 = derivatives(t, s);
    State k2 = derivatives(t + h / 3.0, s + k1 * (h / 3.0));
    State k3 = derivatives(t + 2.0 * h / 3.0, s + k1 * (-h / 3.0) + k2 * h);
    State k4 = derivatives(t + h, s + k1 * h - k2 * h + k3 * h);
    return s + (k1 + k2 * 3.0 + k3 * 3.0 + k4) * (h / 8.0);
}

// 3. Метод Рунге-Кутта 3-го порядка (Метод Кутты)
State stepRK3(double t, const State& s, double h) {
    State k1 = derivatives(t, s);
    State k2 = derivatives(t + h / 2.0, s + k1 * (h / 2.0));
    State k3 = derivatives(t + h, s + k1 * (-h) + k2 * 2.0 * h);
    return s + (k1 + k2 * 4.0 + k3) * (h / 6.0);
}

// 4. DOPRI5 (Дорман-Принс 5(4)) с адаптивным шагом
State stepDOPRI5(double t, const State& s, double& h, double tolerance) {
    const double c2 = 1.0 / 5.0, c3 = 3.0 / 10.0, c4 = 4.0 / 5.0, c5 = 8.0 / 9.0, c6 = 1.0;

    while (true) {
        State k1 = derivatives(t, s);
        State k2 = derivatives(t + c2 * h, s + k1 * (h * (1.0 / 5.0)));
        State k3 = derivatives(t + c3 * h, s + k1 * (h * (3.0 / 40.0)) + k2 * (h * (9.0 / 40.0)));
        State k4 = derivatives(t + c4 * h, s + k1 * (h * (44.0 / 45.0)) - k2 * (h * (56.0 / 15.0)) + k3 * (h * (32.0 / 9.0)));
        State k5 = derivatives(t + c5 * h, s + k1 * (h * (19372.0 / 6561.0)) - k2 * (h * (25360.0 / 2187.0)) + k3 * (h * (64448.0 / 6561.0)) - k4 * (h * (212.0 / 729.0)));
        State k6 = derivatives(t + c6 * h, s + k1 * (h * (9017.0 / 3168.0)) - k2 * (h * (355.0 / 33.0)) + k3 * (h * (46732.0 / 5247.0)) + k4 * (h * (49.0 / 176.0)) - k5 * (h * (5103.0 / 18656.0)));
        
        State s_next = s + (k1 * (35.0 / 384.0) + k3 * (500.0 / 1113.0) + k4 * (125.0 / 192.0) - k5 * (2187.0 / 6784.0) + k6 * (11.0 / 84.0)) * h;
        
        State k7 = derivatives(t + h, s_next);
        
        State error = (k1 * (71.0 / 57600.0) - k3 * (71.0 / 16695.0) + k4 * (71.0 / 1920.0) - k5 * (17253.0 / 339200.0) + k6 * (22.0 / 525.0) - k7 * (1.0 / 40.0)) * h;
        
        double err_norm = std::sqrt(error.x * error.x + error.y * error.y + error.vx * error.vx + error.vy * error.vy);
        if (err_norm <= tolerance) {
                    
                    double scale = 0.9 * std::pow(tolerance / (err_norm + 1e-16), 0.2);
                    scale = std::max(0.1, std::min(5.0, scale));
                    h *= scale;
                    return s_next;
                }
                else {
                    
                    double scale = 0.9 * std::pow(tolerance / err_norm, 0.25);
                    scale = std::max(0.1, std::min(0.5, scale));
                    h *= scale;
                }
            }
        }

        
        double compute_abs_error(const State& current, const State& reference) {
            double dx = current.x - reference.x;
            double dy = current.y - reference.y;
            return std::sqrt(dx * dx + dy * dy);
        }

int main() {
    
    std::cout << std::scientific << std::setprecision(17);
    
    
    double t_end = 5.0;
    double h_fixed = 0.01;
    double ref_tolerance = 1e-12;
    
    std::ifstream config("config.txt");
    if (config.is_open()) {
        config >> t_end >> h_fixed >> ref_tolerance;
        config.close();
    }
    else {
        std::cout << "# Предупреждение: Файл config.txt не найден. Исполняются параметры по умолчанию.\n";
    }
    
    
    State init_state = { 1.0, 0.0, 0.0, 0.5 };
    
    
    struct ReferencePoint {
        double t;
        State state;
    };
    std::vector<ReferencePoint> ref_path;
    
    double t = 0.0;
    State s_ref = init_state;
    double h_adaptive = h_fixed;
    
    ref_path.push_back({ t, s_ref });
    while (t < t_end) {
        if (t + h_adaptive > t_end) h_adaptive = t_end - t;
        s_ref = stepDOPRI5(t, s_ref, h_adaptive, ref_tolerance);
        t += h_adaptive;
        ref_path.push_back({ t, s_ref });
    }
    
    std::cout << "# Сгенерировано эталонных точек: " << ref_path.size() << "\n";
    std::cout << "------------------------------------------------------------------------------------------------------------------------\n";
    std::cout << "Time \t\t Method \t Absolute Error \t Relative Error\n";
    std::cout << "------------------------------------------------------------------------------------------------------------------------\n";
    
   
    State s_rk4 = init_state;
    State s_rk38 = init_state;
    State s_rk3 = init_state;
    State s_dopri5_test = init_state;
    
    double h_test = h_fixed;
    
    for (size_t i = 0; i < ref_path.size() - 1; ++i) {
        double current_t = ref_path[i].t;
        double dt = ref_path[i + 1].t - current_t;
        State ref_next = ref_path[i + 1].state;
        double ref_norm = std::sqrt(ref_next.x * ref_next.x + ref_next.y * ref_next.y);
        
        s_rk4 = stepRK4(current_t, s_rk4, dt);
        s_rk38 = stepRK38(current_t, s_rk38, dt);
        s_rk3 = stepRK3(current_t, s_rk3, dt);
        s_dopri5_test = stepDOPRI5(current_t, s_dopri5_test, h_test, 1e-5);
        
        
        if (i % 10 == 0 || i == ref_path.size() - 2) {
            double abs_rk4 = compute_abs_error(s_rk4, ref_next);
            double abs_rk38 = compute_abs_error(s_rk38, ref_next);
            double abs_rk3 = compute_abs_error(s_rk3, ref_next);
            
            std::cout << current_t + dt << "\t RK4 \t " << abs_rk4 << "\t " << (abs_rk4 / ref_norm) << "\n";
        }
    }
}
