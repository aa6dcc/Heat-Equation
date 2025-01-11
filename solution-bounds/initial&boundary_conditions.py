import numpy as np
import matplotlib.pyplot as plt

def solve_heat_equation_with_bounds(L=1.0, alpha=0.01, x_points=100, t_points=100, T=1.0, 
                                    u0=lambda x: np.sin(np.pi * x), u_left=0, u_right=0):
    """
    Solve the 1D heat equation and compute bounds based on initial and boundary conditions.
    Adjusts time step to meet stability condition if needed.
    """
    # Discretize space
    x = np.linspace(0, L, x_points)
    dx = x[1] - x[0]

    # Adjust time points to ensure stability
    r_target = 0.5
    dt = r_target * dx**2 / alpha
    t_points = int(T / dt) + 1
    t = np.linspace(0, T, t_points)

    # Compute the actual r value
    r = alpha * dt / dx**2

    # Initialize the solution matrix
    u = np.zeros((t_points, x_points))

    u[0, :] = u0(x)

    # Set boundary conditions
    u[:, 0] = u_left
    u[:, -1] = u_right

    for n in range(0, t_points - 1):
        for i in range(1, x_points - 1):
            u[n + 1, i] = u[n, i] + r * (u[n, i - 1] - 2 * u[n, i] + u[n, i + 1])

    max_initial = np.max(u[0, :])
    max_boundary = max(u_left, u_right)
    upper_bound = max(max_initial, max_boundary)
    lower_bound = min(np.min(u[0, :]), min(u_left, u_right))

    return x, t, u, (lower_bound, upper_bound)

L = 1.0
alpha = 0.01
x_points = 100
t_points = 100
T = 1.0

# Initial condition: sin(pi * x)
initial_condition = lambda x: np.sin(np.pi * x)
u_left_boundary = 0  # Boundary condition at x=0
u_right_boundary = 0  # Boundary condition at x=L

x, t, u, bounds = solve_heat_equation_with_bounds(L, alpha, x_points, t_points, T,
                                                  u0=initial_condition,
                                                  u_left=u_left_boundary,
                                                  u_right=u_right_boundary)

plt.figure(figsize=(8, 6))
for n in range(0, len(t), max(1, len(t) // 5)):
    plt.plot(x, u[n, :], label=f"t = {t[n]:.2f}")

plt.axhline(bounds[0], color='red', linestyle='--', label=f"Lower Bound = {bounds[0]:.2f}")
plt.axhline(bounds[1], color='green', linestyle='--', label=f"Upper Bound = {bounds[1]:.2f}")

plt.title("Heat Equation Solution with Bounding Using Initial and Boundary Conditions")
plt.xlabel("x")
plt.ylabel("Temperature u(x, t)")
plt.legend()
plt.grid(True)
plt.show()
