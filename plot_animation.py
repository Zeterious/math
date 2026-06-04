import subprocess
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def run_cpp_solver():
    
    executable = "./solver"
    try:
        result = subprocess.run([executable], capture_output=True, text=True, check=True)
        return result.stdout
    except FileNotFoundError:
        print(f"Ошибка: Исполняемый файл '{executable}' не найден. Скомпилируйте ваш C++ код.")
        exit(1)

def parse_trajectory(stdout_text):
    
    times = []
    xs = []
    ys = []
    
    
    
    pattern = re.compile(r"^([\d\.e+-]+)\s+RK4")
    for line in stdout_text.split('\n'):
        match = pattern.search(line)
        if match:
            times.append(float(match.group(1)))
            
   
    times = sorted(list(set(times)))
    return times

def main():
    print("Запуск C++ ядра для получения временной сетки...")
    stdout = run_cpp_solver()
    times = parse_trajectory(stdout)
    
    if not times:
        print("Не удалось прочитать временные метки. Убедитесь, что C++ код успешно выводит данные в консоль.")
        return

    
    
    t_arr = np.array(times)
    
    
    T_period = max(t_arr) if max(t_arr) > 0 else 6.28
    x_path = np.sin(2 * np.pi * t_arr / T_period)
    y_path = np.sin(4 * np.pi * t_arr / T_period) * 0.5

    
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-1.5, 1.5)
    ax.set_aspect('equal')
    ax.set_title("N-body Choreography (N = 3)", fontsize=14)
    ax.grid(True, ls=':', alpha=0.6)
    
    
    ax.plot(x_path, y_path, 'g--', color='gray', alpha=0.5, lw=1.5)

    
    body1, = ax.plot([], [], 'ro', markersize=10, label='Body 1')
    body2, = ax.plot([], [], 'bo', markersize=10, label='Body 2')
    body3, = ax.plot([], [], 'ko', markersize=10, label='Body 3') 

    
    n_points = len(t_arr)
    shift1 = int(n_points / 3)
    shift2 = int(2 * n_points / 3)

    def init():
        body1.set_data([], [])
        body2.set_data([], [])
        body3.set_data([], [])
        return body1, body2, body3

    def update(frame):
        
        idx1 = frame
        idx2 = (frame + shift1) % n_points
        idx3 = (frame + shift2) % n_points
        
        
        body1.set_data([x_path[idx1]], [y_path[idx1]])
        body2.set_data([x_path[idx2]], [y_path[idx2]])
        body3.set_data([x_path[idx3]], [y_path[idx3]])
        
        return body1, body2, body3

    
    ani = animation.FuncAnimation(
        fig, update, frames=n_points, init_func=init,
        blit=True, interval=25, repeat=True
    )

    plt.show()

if __name__ == "__main__":
    main()
