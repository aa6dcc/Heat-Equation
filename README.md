# Solving the heat equation

The heat equation is a partial differential equation.

It is most commonly seen under its standard form in one dimension:

![Heat Equation](https://latex.codecogs.com/svg.latex?\frac{\partial%20u}{\partial%20t}%20=%20\alpha%20\frac{\partial^2%20u}{\partial%20x^2})

But it in higher dimensions it can be expressed as:

![Heat Equation](https://latex.codecogs.com/svg.latex?\frac{\partial%20u}{\partial%20t}%20=%20\alpha%20\nabla^2%20u)

Where:
- `u(x, t)`: Temperature at position `x` and time `t`
- `α`: Thermal diffusivity of the material
- `∂u/∂t`: Rate of change of temperature with respect to time
- `∂²u/∂x²`: Second spatial derivative
- `∇²`: Laplace operator

It was first developed by French mathematician and physicist Joseph Fourier in 1822 for modelling how a quantity such as heat diffuses in a given region. 
