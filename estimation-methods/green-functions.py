import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

def heat_equation_green(f, L=1.0, alpha=0.01, x_points=100, t_points=100, T=1.0)->None:
    """
    Solve the 1D heat equation using the Green's function method with a user-defined initial condition.

    Parameters:
        f (function): Initial condition function f(x) defined on [0, L].
        L (float): Length of the rod.
        alpha (float): Thermal diffusivity.
        x_points (int): Number of spatial points.
        t_points (int): Number of time points.
        T (float): Total time.
    """
    x = np.linspace(0, L, x_points) 
    t = np.linspace(0, T, t_points) 

    # Green's function for the heat equation in 1D
    def G(x, xi, t):
        if t <= 0:
            return 0
        return np.sqrt(1 / (4 * np.pi * alpha * t)) * np.exp(-((x - xi)**2) / (4 * alpha * t))

    # Compute the solution u(x, t) using Green's function
    def u_xt(x, t):
        solution = np.zeros_like(x)
        for i, xi in enumerate(x):  # Evaluate integral for each point in x
            integrand = lambda x_prime: f(x_prime) * G(xi, x_prime, t)
            integral, _ = quad(integrand, 0, L)  # Integrate over the domain [0, L]
            solution[i] = integral
        return solution

    plt.figure(figsize=(8, 6))
    for time in np.linspace(0, T, 5):  
        u = u_xt(x, time)
        plt.plot(x, u, label=f"t = {time:.2f}")

    plt.title("Heat Equation Solution via Green's Function")
    plt.xlabel("x")
    plt.ylabel("Temperature u(x, t)")
    plt.legend()
    plt.grid(True)
    plt.show()


# Examples with different initial conditions

# Example 1: Initial condition f(x) = sin(pi * x)
heat_equation_green(lambda x: np.sin(np.pi * x), L=1.0, alpha=0.01, x_points=100, t_points=100, T=1.0)

# Example 2: Initial condition f(x) = x * (1 - x)
heat_equation_green(lambda x: x * (1 - x), L=1.0, alpha=0.01, x_points=100, t_points=100, T=1.0)

# Example 3: Initial condition f(x) = 1
heat_equation_green(lambda x: 1, L=1.0, alpha=0.01, x_points=100, t_points=100, T=1.0)
