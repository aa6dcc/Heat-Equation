import numpy as np
import matplotlib.pyplot as plt

def solve_heat_equation(L=1.0, T=1.0, alpha=0.01, nx=50, nt=100, u0=None, boundary_conditions=(0, 0)):
    """
    Solve the 1D heat equation numerically using the finite difference method.
    
    Parameters:
        L (float): Length of the rod.
        T (float): Total simulation time.
        alpha (float): Thermal diffusivity.
        nx (int): Number of spatial points.
        nt (int): Number of time steps.
        u0 (function): Initial condition function f(x). If None, defaults to u(x, 0) = sin(pi * x).
        boundary_conditions (tuple): Dirichlet boundary conditions (u(0, t), u(L, t)).
    
    Returns:
        x (ndarray): Spatial grid points.
        t (ndarray): Time grid points.
        u (ndarray): Solution array of shape (nt, nx).
    """
    dx = L / (nx - 1)
    dt = T / (nt - 1)
    r = alpha * dt / dx**2

    if r > 0.5:
        raise ValueError("The solution is unstable. Choose smaller dt or larger dx.")
    
    x = np.linspace(0, L, nx)
    t = np.linspace(0, T, nt)
    u = np.zeros((nt, nx))

    # Initial condition
    if u0 is None:
        u[0, :] = np.sin(np.pi * x)
    else:
        u[0, :] = u0(x)
    
    # Boundary conditions
    u[:, 0] = boundary_conditions[0]  # u(0, t)
    u[:, -1] = boundary_conditions[1]  # u(L, t)
    
    for n in range(0, nt - 1):
        for i in range(1, nx - 1):
            u[n + 1, i] = u[n, i] + r * (u[n, i - 1] - 2 * u[n, i] + u[n, i + 1])
    
    return x, t, u

def verify_maximum_principle(u, x, t, u0=None, boundary_conditions=(0, 0)):
    """
    Verify the Maximum Principle for the heat equation solution.
    
    Parameters:
        u (ndarray): Solution array of shape (nt, nx).
        x (ndarray): Spatial grid points.
        t (ndarray): Time grid points.
        u0 (function): Initial condition function f(x).
        boundary_conditions (tuple): Boundary conditions (u(0, t), u(L, t)).
    """
    max_initial = np.max(u[0, :])  # Maximum at t = 0
    max_boundary = max(boundary_conditions)  # Maximum on spatial boundaries
    
    max_interior = np.max(u)  # Maximum in the entire domain
    max_allowed = max(max_initial, max_boundary)
    
    print(f"Maximum in domain: {max_interior}")
    print(f"Maximum on boundaries and initial condition: {max_allowed}")
    
    if max_interior > max_allowed:
        print("Maximum Principle violated!")
    else:
        print("Maximum Principle satisfied.")

# Example: Solve and verify Maximum Principle
L = 1.0
T = 1.0
alpha = 0.01
nx = 50
nt = 100

# Solve heat equation with default initial condition and zero boundary conditions
x, t, u = solve_heat_equation(L=L, T=T, alpha=alpha, nx=nx, nt=nt, boundary_conditions=(0, 0))

verify_maximum_principle(u, x, t)

# Plot the solution with Maximum Principle bounds
plt.figure(figsize=(8, 6))
for n in range(0, nt, max(1, nt // 5)):
    plt.plot(x, u[n, :], label=f"t = {t[n]:.2f}")

max_initial = np.max(u[0, :])
max_boundary = max(0, 0) 
max_allowed = max(max_initial, max_boundary)

plt.axhline(max_initial, color='red', linestyle='--', label=f"Max Initial = {max_initial:.2f}")
plt.axhline(max_boundary, color='green', linestyle='--', label=f"Max Boundary = {max_boundary:.2f}")

plt.title("Heat Equation Solution with Maximum Principle Bounds")
plt.xlabel("x")
plt.ylabel("u(x, t)")
plt.legend()
plt.grid(True)
plt.show()
