import numpy as np
import matplotlib.pyplot as plt

def explicit_method_heat_equation(L, T, alpha, x_points, t_points, u0, u_left, u_right):
    dx = L / (x_points - 1)
    dt = T / (t_points - 1)
    r = alpha * dt / dx**2

    if r > 0.5:
        raise ValueError(f"Stability condition not met: r = {r} > 0.5")

    x = np.linspace(0, L, x_points)
    t = np.linspace(0, T, t_points)
    u = np.zeros((t_points, x_points))

    u[0, :] = u0(x)
    u[:, 0] = u_left
    u[:, -1] = u_right

    for n in range(0, t_points - 1):
        for i in range(1, x_points - 1):
            u[n + 1, i] = u[n, i] + r * (u[n, i - 1] - 2 * u[n, i] + u[n, i + 1])

    return x, t, u

L, T, alpha = 1.0, 0.5, 0.01
x_points, t_points = 50, 500
u0 = lambda x: np.sin(np.pi * x)
u_left, u_right = 0, 0

x, t, u = explicit_method_heat_equation(L, T, alpha, x_points, t_points, u0, u_left, u_right)

plt.imshow(u, extent=[0, L, 0, T], origin='lower', aspect='auto', cmap='hot')
plt.colorbar(label="Temperature")
plt.title("Explicit Method Solution")
plt.xlabel("Position (x)")
plt.ylabel("Time (t)")
plt.show()
