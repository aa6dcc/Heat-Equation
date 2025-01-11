import numpy as np
import matplotlib.pyplot as plt

def solve_heat_equation_energy_bounds(L, alpha, x_points, t_points, T, u0, u_left, u_right):
    """
    Solve the 1D heat equation and calculate energy-based bounds using the energy method.

    Parameters:
        L (float): Length of the domain.
        alpha (float): Thermal diffusivity.
        x_points (int): Number of spatial points.
        t_points (int): Number of time points.
        T (float): Total time duration.
        u0 (callable): Initial condition function, u(x, 0).
        u_left (float): Boundary condition at x=0.
        u_right (float): Boundary condition at x=L.

    Returns:
        x, t, u, energy: Spatial points, time points, solution matrix, and energy at each time step.
    """
    dx = L / (x_points - 1)
    dt = T / (t_points - 1)
    r = alpha * dt / dx**2

    if r > 0.5:
        raise ValueError(f"Stability condition not met: r = {r} > 0.5")

    # Discretize space and time
    x = np.linspace(0, L, x_points)
    t = np.linspace(0, T, t_points)
    u = np.zeros((t_points, x_points))

    # Initial and boundary conditions
    u[0, :] = u0(x)
    u[:, 0] = u_left
    u[:, -1] = u_right

    # Time stepping (FTCS scheme)
    for n in range(0, t_points - 1):
        for i in range(1, x_points - 1):
            u[n + 1, i] = u[n, i] + r * (u[n, i - 1] - 2 * u[n, i] + u[n, i + 1])

    # Calculate energy
    energy = np.array([np.trapz(u[n, :]**2, x) for n in range(t_points)])

    return x, t, u, energy

L = 1.0
alpha = 0.01
x_points = 100
t_points = 50000  
T = 100.0         

# Define initial and boundary conditions
initial_condition = lambda x: np.sin(np.pi * x)  # Initial condition: sin(pi * x)
u_left_boundary = 0.0  # Boundary condition at x=0
u_right_boundary = 0.0  # Boundary condition at x=L

# Solve the heat equation
x, t, u, energy = solve_heat_equation_energy_bounds(L, alpha, x_points, t_points, T,
                                                    u0=initial_condition,
                                                    u_left=u_left_boundary,
                                                    u_right=u_right_boundary)

fig, ax = plt.subplots(figsize=(10, 6))

ax.plot(t, energy, label="Energy ||u(x,t)||^2", color='blue')
ax.set_title("Energy Decay and Bounds (Extended Time Frame)")
ax.set_xlabel("Time (t)")
ax.set_ylabel("Energy")

ax.fill_between(t, 0, energy, color='blue', alpha=0.2, label="Energy Bounds")
ax.set_yscale('log')  # Log scale for exponential decay visualization
ax.set_ylim(1e-8, max(energy) * 1.1)  
ax.set_xlim(0, T)  
ax.legend()
ax.grid(True)

plt.tight_layout()
plt.show()
