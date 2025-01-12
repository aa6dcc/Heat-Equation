import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go

def monte_carlo_heat_equation(L, T, x_points, t_points, n_particles, n_steps):
    """
    Solve the heat equation using Monte Carlo simulations.

    Parameters:
        L (float): Length of the domain.
        T (float): Total simulation time.
        x_points (int): Number of spatial points.
        t_points (int): Number of time points.
        n_particles (int): Number of particles for the Monte Carlo simulation.
        n_steps (int): Number of steps per particle.
        
    Returns:
        x (np.ndarray): Spatial points.
        t (np.ndarray): Time points.
        u (np.ndarray): Solution array u(x, t).
    """
    # Discretize space and time
    x = np.linspace(0, L, x_points)
    t = np.linspace(0, T, t_points)
    dx = x[1] - x[0]
    dt = T / (t_points - 1)
    alpha = dx**2 / (2 * dt)  # Effective diffusivity

    # Initialize solution
    u = np.zeros((t_points, x_points))

    u[0, :] = np.exp(-100 * (x - L / 2)**2)  # Gaussian peak at the center

    # Monte Carlo simulation for each time point
    for n in range(1, t_points):
        for i in range(x_points):
            sum_temp = 0
            for _ in range(n_particles):
                pos = x[i]  # Start particle at position x[i]
                for _ in range(n_steps):
                    pos += np.random.choice([-dx, dx])  # Random walk step
                    if pos <= 0 or pos >= L:  # Reflective boundary conditions
                        break
                sum_temp += np.exp(-100 * (pos - L / 2)**2)  # Contribution
            u[n, i] = sum_temp / n_particles

    return x, t, u


L = 1.0          # Length of the domain
T = 0.1          # Total simulation time
x_points = 100   # Number of spatial points
t_points = 50    # Number of time points
n_particles = 50  # Number of particles for Monte Carlo
n_steps = 10    # Number of steps per particle

x, t, u = monte_carlo_heat_equation(L, T, x_points, t_points, n_particles, n_steps)

plt.figure(figsize=(8, 6))
plt.imshow(u.T, extent=[0, T, 0, L], origin='lower', aspect='auto', cmap='hot')
plt.colorbar(label='Temperature')
plt.title('Heat Equation Solution via Monte Carlo Simulation')
plt.xlabel('Time (t)')
plt.ylabel('Position (x)')
plt.show()

plt.figure(figsize=(8, 6))
time_index = int(t_points / 2)
plt.plot(x, u[time_index, :], label=f'Temperature at t={t[time_index]:.2f}')
plt.title('Temperature Distribution')
plt.xlabel('Position (x)')
plt.ylabel('Temperature')
plt.legend()
plt.grid(True)
plt.show()

def monte_carlo_heat_eq_plotly(num_particles=1000, num_steps=100, domain_length=1.0, dt=0.01, dx=0.01, alpha=0.01):
    """
    Monte Carlo simulation for the heat equation with Plotly visualization.
    
    Parameters:
        num_particles (int): Number of particles in the simulation.
        num_steps (int): Number of time steps.
        domain_length (float): Length of the spatial domain.
        dt (float): Time step size.
        dx (float): Spatial step size.
        alpha (float): Thermal diffusivity.
    """
    # Discretize the spatial domain
    x_points = int(domain_length / dx) + 1
    x = np.linspace(0, domain_length, x_points)
    
    # Particle positions, initialize uniformly in the domain
    particle_positions = np.random.uniform(0, domain_length, size=num_particles)
    
    # Precompute probabilities for random walk
    jump_prob = alpha * dt / dx**2  # Ensure this satisfies stability conditions
    if jump_prob > 0.5:
        raise ValueError("Jump probability too high! Reduce dt or increase dx for stability.")
    
    # Array to store particle counts
    particle_density = np.zeros((num_steps, x_points))
    
    # Monte Carlo simulation: track particle movements
    for step in range(num_steps):
        # Update particle positions using random walk
        random_moves = np.random.choice([-1, 0, 1], size=num_particles, p=[jump_prob / 2, 1 - jump_prob, jump_prob / 2])
        particle_positions += random_moves * dx
        
        # Reflective boundary conditions
        particle_positions = np.clip(particle_positions, 0, domain_length)
        
        # Count particles in spatial bins
        particle_counts, _ = np.histogram(particle_positions, bins=x_points, range=(0, domain_length))
        particle_density[step, :] = particle_counts / num_particles
    
    # Create an animated Plotly heatmap
    frames = []
    for step in range(num_steps):
        frames.append(go.Frame(
            data=[go.Heatmap(z=[particle_density[step]], x=x, y=[""], colorscale="Viridis", showscale=True)],
            name=f"t={step*dt:.2f}"
        ))
    
    fig = go.Figure(
        data=[go.Heatmap(z=[particle_density[0]], x=x, y=[""], colorscale="Viridis", showscale=True)],
        layout=go.Layout(
            title="Monte Carlo Simulation of Heat Equation",
            xaxis=dict(title="Spatial Domain (x)"),
            yaxis=dict(title="Temperature Density"),
            updatemenus=[
                dict(
                    type="buttons",
                    showactive=False,
                    buttons=[
                        dict(label="Play", method="animate", args=[None, dict(frame=dict(duration=50, redraw=True))]),
                        dict(label="Pause", method="animate", args=[[None], dict(frame=dict(duration=0, redraw=False))])
                    ]
                )
            ]
        ),
        frames=frames
    )
    
    fig.show()

monte_carlo_heat_eq_plotly(num_particles=1000, num_steps=100, domain_length=1.0, dt=0.001, dx=0.02, alpha=0.01)
