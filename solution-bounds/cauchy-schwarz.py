import numpy as np
import matplotlib.pyplot as plt

def solve_heat_equation_with_cauchy_schwarz(L=1.0, alpha=0.01, x_points=100, t_points=100, T=1.0, 
                                            u0=lambda x: np.sin(np.pi * x), u_left=0, u_right=0):
    """
    Solve the 1D heat equation and compute bounds using the Cauchy-Schwarz inequality.

    Parameters:
        L (float): Length of the rod.
        alpha (float): Thermal diffusivity.
        x_points (int): Number of spatial points.
        t_points (int): Number of time points.
        T (float): Total time.
        u0 (function): Initial condition function f(x) defined on [0, L].
        u_left (float): Boundary condition at x=0.
        u_right (float): Boundary condition at x=L.

    Returns:
        x (np.ndarray): Spatial points.
        t (np.ndarray): Time points.
        u (np.ndarray): Solution array u(x, t).
        energy (np.ndarray): Energy ||u(x, t)||^2 over time.
        bounds (np.ndarray): Array of bounds derived using Cauchy-Schwarz.
    """
    # Discretize space and time
    x = np.linspace(0, L, x_points)
    dx = x[1] - x[0]
    dt = T / (t_points - 1)
    t = np.linspace(0, T, t_points)

    # Stability condition
    r = alpha * dt / dx**2
    if r > 0.5:
        raise ValueError(f"Stability condition not met: r = {r} > 0.5")

    # Initialize the solution matrix
    u = np.zeros((t_points, x_points))

    # Set the initial condition
    u[0, :] = u0(x)

    # Set boundary conditions
    u[:, 0] = u_left
    u[:, -1] = u_right

    # Time-stepping to solve the heat equation
    for n in range(0, t_points - 1):
        for i in range(1, x_points - 1):
            u[n + 1, i] = u[n, i] + r * (u[n, i - 1] - 2 * u[n, i] + u[n, i + 1])

    # Compute energy and bounds using Cauchy-Schwarz
    energy = np.array([np.sum(u[n, :]**2) * dx for n in range(t_points)])  # Energy ||u(x, t)||^2
    bounds = np.sqrt(energy)  # Cauchy-Schwarz: ||u v|| ≤ ||u|| ||v||

    return x, t, u, energy, bounds


# Define parameters
L = 1.0
alpha = 0.01
x_points = 50
t_points = 10000
T = 100.0

# Initial condition: sin(pi * x)
initial_condition = lambda x: np.sin(np.pi * x)
u_left_boundary = 0  # Boundary condition at x=0
u_right_boundary = 0  # Boundary condition at x=L

# Solve the heat equation
x, t, u, energy, bounds = solve_heat_equation_with_cauchy_schwarz(L, alpha, x_points, t_points, T,
                                                                  u0=initial_condition,
                                                                  u_left=u_left_boundary,
                                                                  u_right=u_right_boundary)

# Plot the superposed graph
plt.figure(figsize=(10, 6))
plt.plot(t, energy, label="Energy ||u(x,t)||^2", color='blue')
plt.plot(t, bounds, label="Cauchy-Schwarz Bound", linestyle='--', color='orange')
plt.title("Energy and Cauchy-Schwarz Bounds Over Time")
plt.xlabel("Time (t)")
plt.ylabel("Magnitude")
plt.legend()
plt.grid(True)
plt.show()
