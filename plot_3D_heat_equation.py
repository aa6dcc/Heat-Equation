import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Define the solution to the heat equation
def heat_solution_1(x, t, alpha):
    """
    Heat equation solution: u(x, t) = exp(-alpha * pi^2 * t) * sin(pi * x).

    Parameters:
        x (numpy.ndarray): Spatial points.
        t (numpy.ndarray): Time points.
        alpha (float): Thermal diffusivity.

    Returns:
        numpy.ndarray: Solution u(x, t).
    """
    return np.exp(-alpha * np.pi**2 * t) * np.sin(np.pi * x)

L = 1.0       # Length of the rod
alpha = 0.01  # Thermal diffusivity
x_points = 100
t_points = 100
T = 1.0       # Total time

# Discretize space and time
x = np.linspace(0, L, x_points)
t = np.linspace(0, T, t_points)
X, T_grid = np.meshgrid(x, t)  # Create grid for 3D plotting

U = heat_solution_1(X, T_grid, alpha)

fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')

surf = ax.plot_surface(X, T_grid, U, cmap='viridis', edgecolor='none')
ax.set_title("3D Heat Equation Solution", fontsize=16)
ax.set_xlabel("x (Position)", fontsize=12)
ax.set_ylabel("t (Time)", fontsize=12)
ax.set_zlabel("u(x, t) (Temperature)", fontsize=12)

fig.colorbar(surf, ax=ax, shrink=0.5, aspect=10) #color bar

plt.show()

# Define the heat equation solution
def heat_solution_2(x, t, alpha, L):
    return np.exp(-alpha * (np.pi**2) * t) * np.sin(np.pi * x / L)

L = 1.0          # Length of the rod
alpha = 1.0      # Diffusion coefficient
x_points = 100   # Number of spatial points
t_points = 100   # Number of time points
T = 0.5          # Maximum time

# Discretize space and time
x = np.linspace(0, L, x_points)
t = np.linspace(0, T, t_points)
X, T_grid = np.meshgrid(x, t)

U = heat_solution_2(X, T_grid, alpha, L)

fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d', facecolor='lightgrey')

# Surface plot with color gradient
surf = ax.plot_surface(T_grid, X, U, cmap='jet', edgecolor='k', linewidth=0.5)

ax.set_title("PDE Solution with a=L=1", fontsize=16)
ax.set_xlabel("Time", fontsize=12)
ax.set_ylabel("Rod", fontsize=12)
ax.set_zlabel("Heat", fontsize=12)

ax.set_xticks([0, 0.1, 0.2, 0.3, 0.4, 0.5])
ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
ax.set_zticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])

fig.colorbar(surf, shrink=0.5, aspect=10)

plt.show()
