def implicit_method_heat_equation(L, T, alpha, x_points, t_points, u0, u_left, u_right):
    dx = L / (x_points - 1)
    dt = T / (t_points - 1)
    r = alpha * dt / dx**2

    x = np.linspace(0, L, x_points)
    t = np.linspace(0, T, t_points)
    u = np.zeros((t_points, x_points))

    u[0, :] = u0(x)
    u[:, 0] = u_left
    u[:, -1] = u_right

    A = np.zeros((x_points - 2, x_points - 2))
    np.fill_diagonal(A, 1 + 2 * r)
    np.fill_diagonal(A[:-1, 1:], -r)
    np.fill_diagonal(A[1:, :-1], -r)

    for n in range(0, t_points - 1):
        b = u[n, 1:-1]
        u_next = np.linalg.solve(A, b)
        u[n + 1, 1:-1] = u_next

    return x, t, u

x, t, u = implicit_method_heat_equation(L, T, alpha, x_points, t_points, u0, u_left, u_right)

plt.imshow(u, extent=[0, L, 0, T], origin='lower', aspect='auto', cmap='hot')
plt.colorbar(label="Temperature")
plt.title("Implicit Method Solution")
plt.xlabel("Position (x)")
plt.ylabel("Time (t)")
plt.show()
