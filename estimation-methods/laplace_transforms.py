import numpy as np
import matplotlib.pyplot as plt
from sympy import symbols, laplace_transform, inverse_laplace_transform, sin, pi, exp, Function

def solve_heat_equation_laplace(f, L=1.0, alpha=0.01, x_points=100, t_points=100, T=1.0):
    """
    Solve the 1D heat equation using Laplace transform and separation of variables.

    Parameters:
        f (function): Initial condition function f(x) defined on [0, L].
        L (float): Length of the rod.
        alpha (float): Thermal diffusivity.
        x_points (int): Number of spatial points.
        t_points (int): Number of time points.
        T (float): Total time.
    """
    x, t, s = symbols('x t s')
    U = Function('U')(x, s)

    u_init = f(x)

    # Laplace-transformed equation: sU(x, s) - f(x) = alpha * d^2U(x, s)/dx^2
    U_general = u_init / (s + (pi ** 2 * alpha))  # Derived analytical solution for U(x, s)

    # Perform the inverse Laplace transform to get u(x, t)
    u_solution = inverse_laplace_transform(U_general, s, t)

    x_vals = np.linspace(0, L, x_points)
    t_vals = np.linspace(0, T, t_points)

    u_numeric = np.zeros((t_points, x_points))
    for i, t_val in enumerate(t_vals):
        for j, x_val in enumerate(x_vals):
            u_numeric[i, j] = u_solution.subs({x: x_val, t: t_val}).evalf()

    plt.figure(figsize=(8, 6))
    for i in range(0, t_points, max(1, t_points // 5)):
        plt.plot(x_vals, u_numeric[i, :], label=f"t = {t_vals[i]:.2f}")

    plt.title("Heat Equation Solution via Laplace Transform")
    plt.xlabel("x")
    plt.ylabel("Temperature u(x, t)")
    plt.legend()
    plt.grid(True)
    plt.show()

# Example 1: Initial condition f(x) = sin(pi * x)
solve_heat_equation_laplace(lambda x: sin(pi * x), L=1.0, alpha=0.01, x_points=100, t_points=100, T=1.0)

# Example 2: Initial condition f(x) = x * (1 - x)
solve_heat_equation_laplace(lambda x: x * (1 - x), L=1.0, alpha=0.01, x_points=100, t_points=100, T=1.0)
