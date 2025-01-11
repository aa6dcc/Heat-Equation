# Solving the heat equation

The heat equation is a partial differential equation.
It tells us that the rate at which the temperature changes at a given point over time depends on the second derivative of that temperature at that point with respect to space. 

It is most commonly seen under its standard form in one dimension:

<p align="center">
    <img src="https://latex.codecogs.com/svg.latex?\frac{\partial%20u}{\partial%20t}%20=%20\alpha%20\frac{\partial^2%20u}{\partial%20x^2}" alt="1D Heat Equation">
    <br>
</p>

But it in higher dimensions it can be expressed as:

<p align="center">
    <img src="https://latex.codecogs.com/svg.latex?\frac{\partial%20u}{\partial%20t}%20=%20\alpha%20\nabla^2%20u" alt="Higher Dimensional Heat Equation">
    <br>
</p>

Where:
- `u(x, t)`: Temperature at position `x` and time `t`
- `α`: Thermal diffusivity of the material
- `∂u/∂t`: Rate of change of temperature with respect to time
- `∂²u/∂x²`: Second spatial derivative
- `∇²`: Laplace operator

It was first developed by French mathematician and physicist Joseph Fourier in 1822 for modelling how a quantity such as heat diffuses in a given region, but many great minds (such as Joule, Carnot, Kelvin... or even Benjamin Thompson, who measured heat produced during the process of boring a cannon) had worked on the concept of heat diffusion, from the Ancient Greeks to modern scientists. 

The heat equation doesn't fully describe the behaviour of heat, as we also need initial conditions and boundary conditions to fully analyze it. 

| Summary   |
|-----------|
| [Introduction](#Introduction) |
| [Newton-Raphson estimation method](#Newton-Raphson-estimation-method) |
| [Finding solution bounds](#Finding-solution-bounds) |
| [Deriving the heat equation in 1D](#Deriving-the-heat-equation-in-1-D)|
| [Generalizing the solution technique](#Generalizing-the-solution-technique) |
| [Computational methods](#Computational-methods) |
| [Optimization](#Optimization) |

## Introduction

As mentioned above, the heat equation is a partial differential equation which arises in problems of heat conduction.

At steady state (a state in which the system doesn't change through time), the heat equation simplifies down to Laplace's equation, a second order partial differential equation:

<p align="center">
    <img src="https://latex.codecogs.com/svg.latex?\nabla^2u%20=%200" alt="Laplace's Equation">
</p>

Where:
- `u`: Scalar potential function (e.g., temperature, electric potential, etc.)
- `∇²`: Laplace operator, defined as `∂²/∂x² + ∂²/∂y²` in 2D or `∂²/∂x² + ∂²/∂y² + ∂²/∂z²` in 3D

In effect, the heat equation can be thought of as an extension of Laplace's equation to include the effect of time-dependent heat flow. 

## Newton-Raphson estimation method

## Finding solution bounds

## Deriving the heat equation in 1D

### Heat flux 
Consider the temperature distribution u(x,t) in a thin metal rod of length L. 

<p align="center">
    <img src="https://github.com/user-attachments/assets/1f24faf7-305e-40d1-85f6-eb92d40e4cb2" alt="Metal Rod">
    <br>
</p>

The rate of change of heat energy in time is the sum of the heat flux through the boundary to its neighbours and the heat energy generated at any point in space/time. 
This equation is a general statement regarding the conservation of thermal energy. 

<!DOCTYPE html>
<html lang="en">
<body>
    <p>The formula for heat energy is given by:</p>
    <p>
        $$c(x) \rho(x) \frac{\partial u(x, t)}{\partial t} = -\frac{\partial q}{\partial x} + Q(x, t)$$
    </p>
    <p>
        <small>Where:</small>
        <ul>
            <li><code>c(x)</code>: Specific heat</li>
            <li><code>ρ(x)</code>: Density</li>
            <li><code>u(x, t)</code>: Internal energy</li>
            <li><code>q</code>: Heat flux</li>
            <li><code>Q(x, t)</code>: Heat source</li>
        </ul>
    </p>
</body>
</html>

### Fourier's law of heat conduction

Fourier worked on the heat equation (in fact he invented the Fourier transform as a tool to solve the heat equation) and derived Fourier's law of heat conduction.
Whilst he was analyzing the heat equation, he made the following three observations: 
    - There is no heat flux when the temperature is constant
    - Heat energy flows from high temperature to low temperature
    - There is more heat flux when there is a larger temperature difference

<!DOCTYPE html>
<html lang="en">
<body>
    <p>The heat flux is defined as:</p>
    <p>
        $$q(x, t) = -k \frac{\partial u}{\partial x}$$
    </p>
    <p>The rate of change of heat flux tells us how the heat varies from x to (x + Δx) and is given by:</p>
    <p>
        $$-\frac{\partial q}{\partial x} = \lim_{\Delta x \to 0} \frac{q(x, t) - q(x + \Delta x, t)}{\Delta x}$$
    </p>
    <p>
        <small>Where:</small>
        <ul>
            <li><code>q(x, t)</code>: Heat flux (thermal energy per unit time)</li>
            <li><code>k</code>: Thermal conductivity</li>
            <li><code>u(x, t)</code>: Temperature</li>
            <li><code>Δx</code>: Small spatial increment</li>
        </ul>
    </p>
</body>
</html>

### Overall 1D solution

When we combine the information above, we get:

<!DOCTYPE html>
<html lang="en">
<body>
    <p>The heat conduction equation is given by:</p>
    <p>
        $$c(x) \rho(x) \frac{\partial u}{\partial t} = k \frac{\partial^2 u}{\partial x^2} + Q$$
    </p>
    <p>Assuming c, ρ, and k are constant in space, the equation simplifies to:</p>
    <p>
        $$\frac{\partial u}{\partial t} = \frac{k}{c \rho} \frac{\partial^2 u}{\partial x^2} + \frac{1}{c \rho} Q$$
    </p>
    <p>since we have that</p>
    <p>
        $$- \frac{\partial}{\partial x} \left( -K \frac{\partial u}{\partial x} \right) = +K \frac{\partial^2 u}{\partial x^2}$$
    </p>
    <p>
        <ul>
            <li><code>c(x)</code>: Specific heat</li>
            <li><code>ρ(x)</code>: Density</li>
            <li><code>u(x, t)</code>: Temperature</li>
            <li><code>k</code>: Thermal conductivity</li>
            <li><code>Q</code>: Heat source</li>
            <li><code>α^2 = K/(cρ)</code>: Thermal diffusivity</li>
        </ul>
    </p>
</body>
</html>


## Generalizing the solution technique

An idealized solution to the heat equation 

<p align="center">
    <img src="https://latex.codecogs.com/svg.latex?\frac{\partial%20u}{\partial%20t}%20=%20\alpha%20\frac{\partial^2%20u}{\partial%20x^2}" alt="1D Heat Equation">
    <br>
</p>

would be plugging in 

<!DOCTYPE html>
<html lang="en">
<head>
    <script src="https://polyfill.io/v3/polyfill.min.js?features=es6"></script>
    <script type="text/javascript" id="MathJax-script" async
        src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js">
    </script>
</head>
<body>
    <h1>Solution to the Heat Equation</h1>
    <p>The heat equation in one dimension is:</p>
    <p>
        $$\frac{\partial u}{\partial t} = \alpha \frac{\partial^2 u}{\partial x^2}$$
    </p>
    <p>Let's test the proposed solution:</p>
    <p>
        $$u(x, t) = \sin(x) e^{-\alpha t}$$
    </p>
    <h2>Step 1: Compute \(\frac{\partial u}{\partial t}\)</h2>
    <p>
        $$\frac{\partial u}{\partial t} = \sin(x) \cdot (-\alpha) e^{-\alpha t}$$
    </p>
    <h2>Step 2: Compute \(\frac{\partial^2 u}{\partial x^2}\)</h2>
    <p>First, find \(\frac{\partial u}{\partial x}\):</p>
    <p>
        $$\frac{\partial u}{\partial x} = \cos(x) e^{-\alpha t}$$
    </p>
    <p>Then, find \(\frac{\partial^2 u}{\partial x^2}\):</p>
    <p>
        $$\frac{\partial^2 u}{\partial x^2} = -\sin(x) e^{-\alpha t}$$
    </p>
    <h2>Step 3: Substitute into the heat equation</h2>
    <p>Substitute \( \frac{\partial u}{\partial t} \) and \( \frac{\partial^2 u}{\partial x^2} \) into the heat equation:</p>
    <p>
        $$\sin(x)(-\alpha)e^{-\alpha t} = \alpha(-\sin(x)e^{-\alpha t})$$
    </p>
    <p>Simplify:</p>
    <p>
        $$-\alpha = -\alpha$$
    </p>
    <h2>Hence we can see a product of a sine and exponential function is a solution. </h2>
</body>
</html>


### Fourier series

## Computational methods

## Optimization


