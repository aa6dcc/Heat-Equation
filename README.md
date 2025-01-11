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

The heat equation doesn't fully describe the behaviour of heat, as we also need initial conditions and Dirichlet boundary conditions to fully analyze it. 

| Summary   |
|-----------|
| [Introduction](#Introduction) |
| [Estimation methods](#Estimation-methods) |
| [Finding solution bounds](#Finding-solution-bounds) |
| [Deriving the heat equation in 1D](#Deriving-the-heat-equation-in-1-D)|
| [Generalizing the solution technique](#Generalizing-the-solution-technique) |
| [Computational methods](#Computational-methods) |
| [Optimization](#Optimization) |

## Introduction

As mentioned above, the heat equation is a partial differential equation which arises in problems of heat conduction.

At steady state (a state in which the system doesn't change through time), the heat equation simplifies down to Laplace's equation, a second order partial differential equation.
This means the solution to Laplace's equation can often serve as the long-term behavior of the heat equation.

<p align="center">
    <img src="https://latex.codecogs.com/svg.latex?\nabla^2u%20=%200" alt="Laplace's Equation">
</p>

Where:
- `u`: Scalar potential function (e.g., temperature, electric potential, etc.)
- `∇²`: Laplace operator, defined as `∂²/∂x² + ∂²/∂y²` in 2D or `∂²/∂x² + ∂²/∂y² + ∂²/∂z²` in 3D

In effect, the heat equation can be thought of as an extension of Laplace's equation to include the effect of time-dependent heat flow. 

## Estimation methods

Before directly solving the heat equation, one could also choose to estimate its solutions.
There are several estimation methods for this, but these only deal with linear systems. 
As the heat equation is a linear PDE, we have to restrict the estimation methods to linear ones (ie. not the Newton-Raphson method for example). 

### Fourier series

A Fourier series is an expansion of a periodic function into an infinite sum of trigonometric functions.

<html>
<h2 style="color:darkblue;">Initial Conditions in the Heat Equation</h2>

<p>The <strong>heat equation</strong> in one dimension is:</p>
<p style="text-align: center;">
  <code>∂u(x,t)/∂t = α ∂²u(x,t)/∂x²</code>
</p>
<p>where:</p>
<ul>
  <li><code>u(x, t)</code> is the temperature at position <code>x</code> and time <code>t</code>.</li>
  <li><code>α</code> is the thermal diffusivity (a constant).</li>
</ul>

<h3 style="color:darkgreen;">Initial Conditions</h3>
<p>To solve this equation, we need to define the temperature distribution at the initial time <code>t = 0</code>. This is called the <em>initial condition</em>:</p>
<p style="text-align: center;">
  <code>u(x, 0) = f(x)</code>
</p>
<p>where <code>f(x)</code> is a function that specifies the initial temperature along the rod.</p>

<h3 style="color:darkgreen;">Why Sine Functions Are Common in Fourier Series Solutions</h3>
<p>The solution to the heat equation often uses a <strong>Fourier series</strong> expansion:</p>
<p style="text-align: center;">
  <code>u(x, t) = Σ<sub>n=1</sub>∞ b<sub>n</sub> sin(nπx/L) exp(-n²π²αt/L²)</code>
</p>

<p>Here:</p>
<ul>
  <li><code>b<sub>n</sub></code> are the Fourier coefficients determined from the initial condition <code>f(x)</code>.</li>
  <li><code>sin(nπx/L)</code> naturally satisfies the boundary conditions <code>u(0, t) = u(L, t) = 0</code> (Dirichlet conditions).</li>
</ul>

<p>The coefficients <code>a<sub>n</sub></code> are zero because sine functions are <strong>odd</strong>, and their integral over a symmetric interval cancels out.</p>


<h3 style="color:darkgreen;">Generalizing the Initial Condition</h3>
<p>Although a sine function is commonly used for <code>f(x)</code>, the method is general and works for any valid <code>f(x)</code>. The Fourier coefficients <code>b<sub>n</sub></code> are computed using:</p>
<p style="text-align: center;">
  <code>b<sub>n</sub> = (2/L) ∫<sub>0</sub><sup>L</sup> f(x) sin(nπx/L) dx</code>
</p>
<p>This allows for initial conditions like:</p>
<ul>
  <li><code>f(x) = sin(πx)</code>: A sine wave initial temperature.</li>
  <li><code>f(x) = x(1 - x)</code>: A parabolic initial temperature.</li>
  <li><code>f(x) = 1</code>: A uniform temperature along the rod.</li>
</ul>
</html>

### Laplace transforms

The Laplace transform is an integral transform that converts a variable in the time domain to a variable in the frequency domain. 

### Green functions

A Green function is the impulse response of linear inhomogeneous PDE differential operators. It's a mathematical tool which describes how a cause an initial time affects an effect at a later time. 

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

    1) There is no heat flux when the temperature is constant
    
    2) Heat energy flows from high temperature to low temperature
    
    3) There is more heat flux when there is a larger temperature difference


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

An idealized solution to the heat equation would be plugging in a product of a sine and exponential function.

<!DOCTYPE html>
<html lang="en">
<body>
    <p>The heat equation in one dimension is:</p>
    <p>
        $$\frac{\partial u}{\partial t} = \alpha \frac{\partial^2 u}{\partial x^2}$$
    </p>
    <p>Let's test the proposed solution:</p>
    <p>
        $$u(x, t) = \sin(x) e^{-\alpha t}$$
    </p>
    <p>Step 1: Compute ∂u/∂t</p>
    <p>
        $$\frac{\partial u}{\partial t} = \sin(x) \cdot (-\alpha) e^{-\alpha t}$$
    </p>
    <p>Step 2: Compute ∂²u/∂x²</p>
    <p>First, find ∂u/∂x::</p>
    <p>
        $$\frac{\partial u}{\partial x} = \cos(x) e^{-\alpha t}$$
    </p>
    <p>Then, find ∂²u/∂x²:</p>
    <p>
        $$\frac{\partial^2 u}{\partial x^2} = -\sin(x) e^{-\alpha t}$$
    </p>
    <p>Step 3: Substitute into the heat equation</p>
    <p>Substitute ∂u/∂t and ∂²u/∂x² into the heat equation:</p>
    <p>
        $$\sin(x)(-\alpha)e^{-\alpha t} = \alpha(-\sin(x)e^{-\alpha t})$$
    </p>
    <p>Simplify:</p>
    <p>
        $$-\alpha = -\alpha$$
    </p>
    <p>Hence we can see a product of a sine and exponential function is a solution. </p>
</body>
</html>

## Computational methods

## Optimization


