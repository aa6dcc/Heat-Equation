import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

def heat_equation_solution(f, L=1.0, alpha=0.01, N=50, x_points=100, t_points=100, T=1.0) -> None:
    """
    Solve the 1D heat equation using Fourier series with a user-defined initial condition.

    Parameters:
        f (function): Initial condition function f(x) defined on [0, L].
        L (float): Length of the rod.
        alpha (float): Thermal diffusivity.
        N (int): Number of terms in the Fourier series.
        x_points (int): Number of spatial points.
        t_points (int): Number of time points.
        T (float): Total time.
    """
    x = np.linspace(0, L, x_points)
    t = np.linspace(0, T, t_points)

    # Compute Fourier coefficients b_n
    def compute_bn(n):
        # Integral of f(x) * sin(n*pi*x/L) over [0, L]
        integrand = lambda x: f(x) * np.sin(n * np.pi * x / L)
        bn, error_estimate = quad(integrand, 0, L)  # Integrate using scipy's quad
        return (2 / L) * bn

    b = [compute_bn(n) for n in range(1, N + 1)]

    # Compute the solution u(x, t)
    def u_xt(x, t):
        solution = np.zeros_like(x)
        for n in range(1, N + 1):
            term = b[n - 1] * np.sin(n * np.pi * x / L) * np.exp(-alpha * (n * np.pi / L)**2 * t)
            solution += term
        return solution

    plt.figure(figsize=(8, 6))
    for time in np.linspace(0, T, 5):  # Plot for 5 time steps
        u = u_xt(x, time)
        plt.plot(x, u, label=f"t = {time:.2f}")

    plt.title("Heat Equation Solution via Fourier Series")
    plt.xlabel("x")
    plt.ylabel("Temperature u(x, t)")
    plt.legend()
    plt.grid(True)
    plt.show()


# Examples with different initial conditions

# Example 1: Initial condition f(x) = sin(pi * x)
heat_equation_solution(lambda x: np.sin(np.pi * x), L=1.0, alpha=0.01, N=50, x_points=100, t_points=100, T=1.0)

# Example 2: Initial condition f(x) = x * (1 - x)
heat_equation_solution(lambda x: x * (1 - x), L=1.0, alpha=0.01, N=50, x_points=100, t_points=100, T=1.0)

# Example 3: Initial condition f(x) = 1
heat_equation_solution(lambda x: 1, L=1.0, alpha=0.01, N=50, x_points=100, t_points=100, T=1.0)
