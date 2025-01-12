import numpy as np
import matplotlib.pyplot as plt

def heat_equation_runge_kutta(L, T, alpha, nx, nt, u0, u_left, u_right):
    """
    Solve the 1D heat equation using Runge-Kutta (RK2) time integration.

    Parameters:
        L (float): Length of the rod.
        T (float): Total time.
        alpha (float): Thermal diffusivity.
        nx (int): Number of spatial grid points.
        nt (int): Number of time steps.
        u0 (callable): Initial condition function u(x, 0).
        u_left (float): Boundary condition at x=0.
        u_right (float): Boundary condition at x=L.

    Returns:
        x (np.ndarray): Spatial grid points.
        t (np.ndarray): Time points.
        u (np.ndarray): Solution array u(x, t).
    """
    # Discretize space and time
    dx = L / (nx - 1)
    dt = T / (nt - 1)
    x = np.linspace(0, L, nx)
    t = np.linspace(0, T, nt)
    
    # Stability condition
    r = alpha * dt / dx**2
    if r > 0.5:
        raise ValueError(f"Stability condition not met: r = {r} > 0.5")

    # Initialize the solution matrix
    u = np.zeros((nt, nx))
    u[0, :] = u0(x)  
    u[:, 0] = u_left  # Boundary condition at x=0
    u[:, -1] = u_right  # Boundary condition at x=L

    # Helper function for the spatial derivative
    def laplacian(u):
        dudx2 = np.zeros_like(u)
        dudx2[1:-1] = (u[:-2] - 2 * u[1:-1] + u[2:]) / dx**2
        return alpha * dudx2

    # Runge-Kutta time stepping
    for n in range(nt - 1):
        k1 = dt * laplacian(u[n, :])
        k2 = dt * laplacian(u[n, :] + 0.5 * k1)
        u[n + 1, :] = u[n, :] + k2

    return x, t, u

L = 1.0  # Length of the rod
T = 1.0  # Total time
alpha = 0.01  # Thermal diffusivity
nx = 50  # Number of spatial points
nt = 200  # Number of time points

initial_condition = lambda x: np.sin(np.pi * x)  # Initial condition: sin(pi * x)
u_left_boundary = 0.0  # Boundary condition at x=0
u_right_boundary = 0.0  # Boundary condition at x=L

x, t, u = heat_equation_runge_kutta(L, T, alpha, nx, nt,
                                    u0=initial_condition,
                                    u_left=u_left_boundary,
                                    u_right=u_right_boundary)

X, T = np.meshgrid(x, t)

fig = plt.figure(figsize=(10, 6))
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, T, u, cmap="viridis")
ax.set_title("Heat Equation Solution Using Runge-Kutta (RK2)")
ax.set_xlabel("Rod Position (x)")
ax.set_ylabel("Time (t)")
ax.set_zlabel("Temperature (u)")
plt.show()
