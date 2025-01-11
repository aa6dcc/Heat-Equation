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

Please find the Python code used for estimations here: [Estimation Methods in Python](https://github.com/aa6dcc/Heat-Equation/blob/main/estimation-methods/Estimation_methods.ipynb)

### Fourier series

A Fourier series is an expansion of a periodic function into an infinite sum of trigonometric functions.
In our case, Fourier series decomposes the initial temperature distribution f(x) into a sum of sine and cosine functions. These basis functions satisfy the boundary conditions and allow the heat equation to be solved as a series of time-evolving terms. Each term represents a mode of heat distribution, with higher modes decaying faster due to diffusion.

<html>
<h4 style="color:darkblue;">Initial Conditions in the Heat Equation</h4>

<p>The <strong>heat equation</strong> in one dimension is:</p>
<p style="text-align: center;">
  <code>∂u(x,t)/∂t = α ∂²u(x,t)/∂x²</code>
</p>
<p>where:</p>
<ul>
  <li><code>u(x, t)</code> is the temperature at position <code>x</code> and time <code>t</code>.</li>
  <li><code>α</code> is the thermal diffusivity (a constant).</li>
</ul>

<h5 style="color:darkgreen;">Initial Conditions</h5>
<p>To solve this equation, we need to define the temperature distribution at the initial time <code>t = 0</code>. This is called the <em>initial condition</em>:</p>
<p style="text-align: center;">
  <code>u(x, 0) = f(x)</code>
</p>
<p>where <code>f(x)</code> is a function that specifies the initial temperature along the rod.</p>

<h5 style="color:darkgreen;">Why Sine Functions Are Common in Fourier Series Solutions</h5>
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


<h5 style="color:darkgreen;">Generalizing the Initial Condition</h5>
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

![image](https://github.com/user-attachments/assets/59b737a1-9264-4e87-8f05-73be0aeaaca9)

Above we plotted the Fourier series for a sine wave as an initial condition.

### Laplace transforms

The Laplace transform is an integral transform that converts a variable in the time domain to a variable in the frequency domain. 
The Laplace transform simplifies the heat equation by converting it from a partial differential equation in time and space to an algebraic equation in the Laplace domain. This approach is particularly useful for problems with time-dependent boundary or initial conditions, as it allows systematic analysis. Once solved, the inverse transform retrieves the time-domain solution, capturing the evolution of heat distribution.

<html lang="en">
  <body>
    <p>The Laplace transform of <code>u(x, t)</code> is defined as:</p>
    <div class="equation">
      <p>U(x, s) = ∫<sub>0</sub><sup>∞</sup> u(x, t) e<sup>-st</sup> dt
      </p>
    </div>
    <h4>The Heat Equation</h4>
    <p>The one-dimensional heat equation is given by:</p>
    <div class="equation">
      <p>
        ∂u(x, t)/∂t = α ∂<sup>2</sup>u(x, t)/∂x<sup>2</sup>
      </p>
    </div>
    <p>with initial condition:</p>
    <div class="equation">
      <p>u(x, 0) = f(x)</p>
    </div>
    <h4>Applying the Laplace Transform</h4>
    <p>
      Taking the Laplace transform of both sides of the heat equation with
      respect to <code>t</code>, we get:
    </p>
    <div class="equation">
      <p>
        sU(x, s) - f(x) = α ∂<sup>2</sup>U(x, s)/∂x<sup>2</sup>
      </p>
    </div>
    <p>Rearranging gives a second-order ordinary differential equation in <code>x</code>:</p>
    <div class="equation">
      <p>
        ∂<sup>2</sup>U(x, s)/∂x<sup>2</sup> - (s/α) U(x, s) = -f(x)/α
      </p>
    </div>
    <h4>Solving the Transformed Equation</h4>
    <p>
      The general solution of this ODE depends on the boundary conditions.
      Once the solution <code>U(x, s)</code> is obtained, the inverse Laplace
      transform is applied to recover <code>u(x, t)</code> in the time domain.
    </p>
  </body>
</html>

![image](https://github.com/user-attachments/assets/49fd437d-6058-4b0b-8b15-dcae363ec8ef)

Above we plotted the Fourier series for a sine wave as an initial condition.

### Green functions

A Green function is the impulse response of linear inhomogeneous PDE differential operators. It's a mathematical tool which describes how a cause an initial time affects an effect at a later time. 
Green's functions offer a way to represent the solution of the heat equation as a convolution with the initial condition, describing how heat propagates over time. They naturally handle boundary conditions and source terms, making them suitable for complex geometries or varying heat sources. Additionally, the method provides a fundamental solution that illustrates how heat spreads from a point source.

<h4 style="color:darkblue;">Green's Function Approach for the Heat Equation</h4>

<p>The Green's function method provides a powerful way to solve the 1D heat equation:</p>

<pre>
∂u(x, t)/∂t = α ∂²u(x, t)/∂x²,
</pre>

<p>with arbitrary initial conditions and Dirichlet boundary conditions.</p>

<h5>Key Ideas:</h5>
<ul>
  <li>
    <strong>Green's Function Definition:</strong>
    The Green's function <code>G(x, ξ, t)</code> represents the temperature at position <code>x</code>
    and time <code>t</code>, due to an initial impulse at position <code>ξ</code>.
  </li>
  <li>
    <strong>Heat Equation Green's Function:</strong>
    For an infinite domain, the Green's function is given by:
    <pre>
    G(x, ξ, t) = (1 / √(4π α t)) * exp(-(x - ξ)² / (4 α t)).
    </pre>
  </li>
  <li>
    <strong>Solution:</strong>
    The temperature <code>u(x, t)</code> is found by integrating the initial condition <code>f(ξ)</code>
    weighted by the Green's function:
    <pre>
    u(x, t) = ∫₀ᴸ f(ξ) G(x, ξ, t) dξ,
    </pre>
    where <code>L</code> is the domain length.
  </li>
</ul>

<h5>Why Use Green's Functions?</h5>
<ul>
  <li>
    Allows solving the heat equation for <strong>arbitrary initial conditions</strong>.
  </li>
  <li>
    Handles time evolution naturally by weighting contributions from the initial temperature distribution.
  </li>
  <li>
    Provides insight into the diffusion process via the exponential decay in the Green's function.
  </li>
</ul>

<h5>Steps to Solve Using Green's Functions:</h5>
<ol>
  <li>Define the Green's function <code>G(x, ξ, t)</code> for the problem.</li>
  <li>Express the solution as <code>u(x, t) = ∫ f(ξ) G(x, ξ, t) dξ</code>.</li>
  <li>Compute the integral for all points <code>x</code> and times <code>t</code>.</li>
</ol>

![image](https://github.com/user-attachments/assets/0d963e9e-eabd-457a-9504-0858174f16a5)

Above we plotted the Fourier series for a sine wave as an initial condition.

## Finding solution bounds

Finding the solution bounds for the heat equation involves estimating the maximum and minimum values of the temperature u(x,t) over the spatial domain and time. 
Instead of using traditional statistical bounds derived from probability inequalities such as Chebyshev's inequality, we will try and use more precise and accurate physical methods. 

### Maximum Principle for Parabolic equations

The heat equation satisfies the maximum principle, which states:
    <li>The maximum and minimum values of the solution u(x,t) occur either at the initial time (t=0) or on the boundary (x=0 or x=L)</li>
    <li>This implies that if u(x,0)=f(x), the bounds of u(x,t) are determined by the bounds of f(x) and any boundary conditions.</li>

### Bounding Using Initial and Boundary Conditions
    <ol> 
    <li>If the initial condition f(x) is known, then u(x,t) is bounded by the maximum and minimum values of f(x).</li>
    <li>For Dirichlet boundary conditions (u(0,t)=u(L,t)=0), u(x,t) remains non-negative if f(x)≥0.</li>
    <li>For Neumann boundary conditions (∂u/∂x=0 at boundaries), u(x,t) may preserve symmetry or maintain constant total heat.</li>
    </ol>

### Energy methods
    <ol>
    <li>The energy (integral of u^2*(x,t)) of the solution decreases over time due to diffusion</li>
    <li>The norm ∥u∥∞ decreases monotonically, providing a practical way to estimate bounds as heat diffuses.</li>
    </ol>
    
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


