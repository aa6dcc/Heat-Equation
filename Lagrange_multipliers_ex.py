#1: give the point (x,y) that minimizes the function f(x,y) = x^2 + y^2 subject to g(x,y) = x + y - 1 = 0 and the value of the Lagrange multiplier λ.

import sympy as sp

# Define variables
x, y, λ = sp.symbols('x y λ')

# Objective function
f = x**2 + y**2

# Constraint
g = x + y - 1

# Lagrangian
L = f + λ * g

# Take partial derivatives
grad_L = [sp.diff(L, var) for var in (x, y, λ)]

# Solve the system of equations
solution = sp.solve(grad_L, (x, y, λ))
print(solution)

#2:give the approximate solution for (x,y) and the minimized value of f(x,y).

from scipy.optimize import minimize # supports constraints and can solve optimization problems numerically

def objective(vars):
    x, y = vars
    return x**2 + y**2

def constraint(vars):
    x, y = vars
    return x + y - 1  # Equality constraint

# Initial guess
x0 = [0.5, 0.5]

constraints = {'type': 'eq', 'fun': constraint}

result = minimize(objective, x0, constraints=constraints)

print("Optimal solution:", result.x)
print("Optimal value:", result.fun)

#3: Minimizing an energy functional J(u) subject to Lu=0 (e.g., −Δu=f):

u, v, λ = sp.symbols('u v λ')
f = sp.Function('f')(u)
constraint = sp.Eq(sp.diff(u, v), f)

L = f + λ * constraint.rhs  # Lagrangian

