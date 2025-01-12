# the code here determines in an input function is a solution to the heat equation

import sympy as sp

def is_solution_to_heat_equation(u_expr, x, t, alpha):
    """
    Test if a given function u(x, t) satisfies the heat equation u_t = alpha * u_xx.

    Parameters:
        u_expr (sympy expression): The candidate solution u(x, t).
        x (sympy.Symbol): The spatial variable.
        t (sympy.Symbol): The time variable.
        alpha (float or sympy.Symbol): Thermal diffusivity.

    Returns:
        bool: True if the function satisfies the heat equation, False otherwise.
        sympy expression: The residual of the equation (u_t - alpha * u_xx).
    """
    # Compute partial derivatives
    u_t = sp.diff(u_expr, t)         # Time derivative
    u_xx = sp.diff(u_expr, x, 2)    # Second spatial derivative

    # Heat equation: u_t = alpha * u_xx
    heat_eq_residual = u_t - alpha * u_xx

    # Simplify the residual to check if it is identically zero
    heat_eq_residual = sp.simplify(heat_eq_residual)

    # Return whether the residual is zero and the residual expression itself
    return heat_eq_residual == 0, heat_eq_residual


# Define symbols
x, t = sp.symbols('x t')
alpha = sp.Symbol('alpha', positive=True, real=True)

# Define candidate solution (example: u(x, t) = exp(-alpha * pi^2 * t) * sin(pi * x))
u_candidate = sp.exp(-alpha * sp.pi**2 * t) * sp.sin(sp.pi * x)

# Test if the candidate solution satisfies the heat equation
is_solution, residual = is_solution_to_heat_equation(u_candidate, x, t, alpha)

# Output results
print(f"Candidate Solution: {u_candidate}")
print(f"Satisfies Heat Equation: {is_solution}")
print(f"Residual: {residual}")
