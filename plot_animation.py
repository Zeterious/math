import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import io

def get_simulation_data():
    executable = "./solver"
    try:
        result = subprocess.run([executable], capture_output=True, text=True, check=True)
        csv_data = "\n".join([line for line in result.stdout.split("\n") if "," in line])
        return pd.read_csv(io.StringIO(csv_data))
    except FileNotFoundError:
        print(f"Ошибка: Исполняемый файл '{executable}' не найден.")
        exit(1)

def main():
    df = get_simulation_data()

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-0.6, 0.6)
    ax.set_aspect('equal')
    ax.set_title("Хореография трёх тел (восьмёрка), m = 1", fontsize=12)
    ax.grid(True, ls=':', alpha=0.5)

    ax.plot(df['b3_x'], df['b3_y'], color='lightgray', lw=1, ls='-')

    body1, = ax.plot([], [], 'ro', markersize=7, label='Тело 1')
    body2, = ax.plot([], [], 'bo', markersize=7, label='Тело 2')
    body3, = ax.plot([], [], 'ko', markersize=7, label='Тело 3')

    def init():
        body1.set_data([], [])
        body2.set_data([], [])
        body3.set_data([], [])
        return body1, body2, body3

    def update(frame):
        row = df.iloc[frame]
        body1.set_data([row['b1_x']], [row['b1_y']])
        body2.set_data([row['b2_x']], [row['b2_y']])
        body3.set_data([row['b3_x']], [row['b3_y']])
        return body1, body2, body3

    ani = animation.FuncAnimation(
        fig, update, frames=len(df), init_func=init,
        blit=True, interval=15, repeat=True
    )

    plt.show()

if __name__ == "__main__":
    main()
