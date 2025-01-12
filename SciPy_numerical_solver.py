from scipy.integrate import solve_ivp
import numpy as np

L = 1.0
alpha = 0.01
x_points = 100
T = 1.0
dx = L / (x_points - 1)
dt = 0.01

x = np.linspace(0, L, x_points)
u0 = np.sin(np.pi * x)  # Initial condition

def heat_eq(t, u):
    dudt = np.zeros_like(u)
    dudt[1:-1] = alpha * (u[:-2] - 2 * u[1:-1] + u[2:]) / dx**2
    return dudt

sol = solve_ivp(heat_eq, [0, T], u0, method='RK45', t_eval=np.linspace(0, T, int(T / dt)))
