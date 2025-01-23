# Solving the heat equation

The heat equation is a partial differential equation.
It tells us that the rate at which the temperature changes at a given point over time depends on the second derivative of that temperature at that point with respect to space. 

Aside from analysing heat diffusion, it has a variety of applications, such as particle diffusion, Brownian motion, image analysis, Riemannian analysis along with financial mathematics for the modelling of options (Black-Scholes can be seen as a variant of the heat equation).

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
| [Bibliography](#Bibliography) |

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

Please find the Python code used for solution bounds here: [Solution bounds in Python](https://github.com/aa6dcc/Heat-Equation/tree/main/solution-bounds)

### Maximum Principle for Parabolic equations

The heat equation satisfies the maximum principle, which states:
    <li>The maximum and minimum values of the solution u(x,t) occur either at the initial time (t=0) or on the boundary (x=0 or x=L)</li>
    <li>This implies that if u(x,0)=f(x), the bounds of u(x,t) are determined by the bounds of f(x) and any boundary conditions.</li>

<p>
    <strong>The Maximum Principle</strong> for parabolic equations is a fundamental result that ensures the solution of 
    the heat equation attains its maximum and minimum values either at the initial time or on the boundaries of the 
    spatial domain. This principle reflects the physical intuition that heat naturally flows from regions of higher 
    temperature to lower temperature, preventing the creation of new extrema within the domain over time. The Maximum 
    Principle is crucial for proving uniqueness, stability, and boundedness of solutions, providing a theoretical 
    foundation for error analysis in numerical approximations and guaranteeing that solutions adhere to physical constraints.
</p>


![image](https://github.com/user-attachments/assets/1093a5cd-9c58-4a89-a32e-cc571cf568fa)

Above is an example of a plot bounded using the maximum principle. 

### Bounding Using Initial and Boundary Conditions

If the initial condition f(x) is known, then u(x,t) is bounded by the maximum and minimum values of f(x):
    <li>For Dirichlet boundary conditions (u(0,t)=u(L,t)=0), u(x,t) remains non-negative if f(x)≥0</li>
    <li>For Neumann boundary conditions (∂u/∂x=0 at boundaries), u(x,t) may preserve symmetry or maintain constant total heat</li>

<p>
    <strong>Bounding the solution</strong> to the heat equation using initial and boundary conditions ensures 
    that the temperature distribution remains physically meaningful and adheres to the problem's constraints. 
    Initial conditions provide the starting energy and structure of the solution, while boundary conditions define 
    the interaction of the system with its surroundings, such as fixed or insulated boundaries. These bounds guarantee 
    that the solution does not exceed the maximum or minimum values dictated by the initial and boundary inputs, 
    reflecting the maximum principle for parabolic equations. This approach is essential for stability analysis, 
    validating numerical approximations, and ensuring the consistency of the solution with physical laws.
</p>


![image](https://github.com/user-attachments/assets/088fea1f-cda8-4732-b12b-3da97f80eacd)

Above is an example of a plot bounded using initial and boundary conditions.

### Energy methods

We can also use the energy of the signal:
    <li>The energy (integral of u^2*(x,t)) of the solution decreases over time due to diffusion</li>
    <li>The norm ∥u∥ of ∞ decreases monotonically, providing a practical way to estimate bounds as heat diffuses</li>

<p>
    <strong>Energy</strong> plays a central role in understanding the behavior of solutions to the heat equation, 
    as it quantifies the total "strength" or magnitude of the temperature distribution over the spatial domain. 
    The energy, defined as <math>
        <msub><mrow>&#x2225;</mrow><mrow>u(x, t)</mrow></msub><sup>2</sup> = <mo>&#x222b;</mo> u(x, t)<sup>2</sup> dx
    </math>, decreases over time due to the dissipative nature of the heat equation, reflecting the physical principle 
    of heat flow from high to low temperatures. By analyzing the energy, we can establish bounds on the solution, 
    ensuring it remains well-behaved and converges as time progresses. Energy estimates are also critical for proving 
    stability, error bounds in numerical methods, and demonstrating the uniqueness of solutions.
</p>

![image](https://github.com/user-attachments/assets/aaae09bc-15aa-4304-8b6f-f1d0b46ae540)

Above is an example of a plot bounded using the energy method

In fact, we can also use Cauchy-Schwarz's inequality.
<p>The <strong>Cauchy-Schwarz inequality</strong> plays a crucial role in analyzing the energy dynamics of the heat equation. 
    The energy, defined as <math>
        <msub><mrow>&#x2225;</mrow><mrow>u(x, t)</mrow></msub><sup>2</sup> = <mo>&#x222b;</mo> u(x, t)<sup>2</sup> dx
    </math>, quantifies the total "magnitude" of the solution over the spatial domain. During analysis, terms like 
    <code>|⟨u, v⟩|</code> appear, representing inner products of functions. Using Cauchy-Schwarz, these terms are bounded as 
    <code>|⟨u, v⟩| ≤ ||u|| ||v||</code>, ensuring numerical stability and aiding in deriving bounds for the solution. This 
    inequality also supports demonstrating exponential decay of energy, highlighting the dissipative nature of the heat equation 
    and ensuring the solution stabilizes over time.</p>

![image](https://github.com/user-attachments/assets/4ad6b432-9cec-44e7-ade6-85d8de180d77)

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

<!DOCTYPE html>
<html lang="en">
<body>
    <h4>Generalizing the Solution Technique</h4>
    <p>
        The solution technique used above can be extended to many other types of equations. 
        The idea is that the operator <span class="math">u<sub>xx</sub></span> with the zero boundary conditions 
        can be represented in terms of its eigenfunctions. This leads naturally to one of the basic ideas of the spectral theory 
        of linear self-adjoint operators.
    </p>
    <p>
        Consider the linear operator <span class="math">Δu = u<sub>xx</sub></span>. The infinite sequence of functions:
    </p>
    <p class="math">
        e<sub>n</sub>(x) = √(2/L) sin(nπx / L), &nbsp; for &nbsp; n ≥ 1
    </p>
    <p>
        are eigenfunctions of <span class="math">Δ</span>. Indeed:
    </p>
    <p class="math">
        Δe<sub>n</sub> = −(n²π² / L²)e<sub>n</sub>.
    </p>
    <p>
        Moreover, any eigenfunction <span class="math">f</span> of <span class="math">Δ</span> with the boundary conditions 
        <span class="math">f(0) = f(L) = 0</span> is of the form <span class="math">e<sub>n</sub></span> for some <span class="math">n ≥ 1</span>. 
        The functions <span class="math">e<sub>n</sub></span> for <span class="math">n ≥ 1</span> form an orthonormal sequence 
        with respect to a certain inner product on the space of real-valued functions on <span class="math">[0, L]</span>. 
        This means:
    </p>
    <p class="math">
        ⟨e<sub>n</sub>, e<sub>m</sub>⟩ = ∫<sub>0</sub><sup>L</sup> e<sub>n</sub>(x)e<sub>m</sub>* (x) dx = δ<sub>mn</sub>
    </p>
    <p>
        Finally, the sequence {e<sub>n</sub>}<sub>n ∈ N</sub> spans a dense linear subspace of 
        <span class="math">L²((0, L))</span>. This shows that in effect we have diagonalized the operator 
        <span class="math">Δ</span>.
    </p>
</body>
</html>

<!DOCTYPE html>
<html lang="en">
<body>
    <h4>General Solution to the Heat Equation</h4>
    <p>The general solution to the one-dimensional heat equation:</p>
    <div class="equation">
        u<sub>t</sub> = α u<sub>xx</sub>
    </div>
    <p>with zero boundary conditions and initial condition <code>u(x, 0) = f(x)</code> is given by:</p>
    <div class="equation">
        u(x, t) = ∑<sub>n=1</sub><sup>∞</sup> b<sub>n</sub> e<sub>n</sub>(x) e<sup>−λ<sub>n</sub>t</sup>
    </div>
    <p>where:</p>
    <ul>
        <li><b>Eigenfunctions:</b> e<sub>n</sub>(x) = √(2/L) sin(nπx/L)</li>
        <li><b>Eigenvalues:</b> λ<sub>n</sub> = α (nπ/L)<sup>2</sup></li>
        <li><b>Fourier Coefficients:</b> b<sub>n</sub> = ∫<sub>0</sub><sup>L</sup> f(x) e<sub>n</sub>(x) dx</li>
    </ul>
    <h5>Steps to Derive the General Solution</h5>
    <ol>
        <li><b>Eigenfunction Expansion:</b> Represent <code>u(x, t)</code> as a sum of eigenfunctions:
            <div class="equation">
                u(x, t) = ∑<sub>n=1</sub><sup>∞</sup> c<sub>n</sub>(t) e<sub>n</sub>(x)
            </div>
        </li>
        <li><b>Substitute into the Heat Equation:</b> Substitute this expansion into <code>u<sub>t</sub> = α u<sub>xx</sub></code>:
            <div class="equation">
                ∑<sub>n=1</sub><sup>∞</sup> c<sub>n</sub>'(t) e<sub>n</sub>(x) = α ∑<sub>n=1</sub><sup>∞</sup> c<sub>n</sub>(t) λ<sub>n</sub> e<sub>n</sub>(x)
            </div>
        </li>
        <li><b>Orthogonality of Eigenfunctions:</b> Use the orthogonality property of <code>e<sub>n</sub>(x)</code> to isolate equations for each mode.</li>
        <li><b>Solve for Time Coefficients:</b> The solution for <code>c<sub>n</sub>(t)</code> is:
            <div class="equation">
                c<sub>n</sub>(t) = c<sub>n</sub>(0) e<sup>−αλ<sub>n</sub>t</sup>
            </div>
        </li>
        <li><b>Determine Initial Coefficients:</b> Use the initial condition <code>u(x, 0) = f(x)</code> to find:
            <div class="equation">
                b<sub>n</sub> = ∫<sub>0</sub><sup>L</sup> f(x) e<sub>n</sub>(x) dx
            </div>
        </li>
        <li><b>Combine Results:</b> Substitute the coefficients into the solution:
            <div class="equation">
                u(x, t) = ∑<sub>n=1</sub><sup>∞</sup> b<sub>n</sub> e<sub>n</sub>(x) e<sup>−αλ<sub>n</sub>t</sup>
            </div>
        </li>
    </ol>
    <h5>Key Properties of the Solution</h5>
    <ul>
        <li><b>Decay of Modes:</b> Higher modes decay faster due to the exponential term <code>e<sup>−λ<sub>n</sub>t</sup></code>.</li>
        <li><b>Orthogonality:</b> The eigenfunctions form an orthonormal basis, ensuring convergence of the solution.</li>
        <li><b>Spectral Representation:</b> The solution decomposes the initial condition into modes of the operator <code>u<sub>xx</sub></code>.</li>
    </ul>
</body>
</html>

## Computational methods

Computational methods are essential for solving the heat equation, especially when analytical solutions are unavailable due to complex geometries, non-linearities, or boundary conditions. These methods discretize the equation into manageable components for numerical approximation.

Tool to verify if a given function is a solution to the heat equation using Sympy: [Code](https://github.com/aa6dcc/Heat-Equation/blob/main/verify_solution_heat_equation.py)

Python scripts running the different computational methods: [Relevant folder](https://github.com/aa6dcc/Heat-Equation/tree/main/computation-methods)

SciPy ODE solver: [Function](https://github.com/aa6dcc/Heat-Equation/blob/main/SciPy_numerical_solver.py)

### Finite Difference Methods

Finite methods are numerical techniques used to approximate the solutions of differential equations. They involve breaking down a continuous domain (such as a region in space or time) into a discrete set of points, enabling the equations to be solved computationally. 

![image](https://github.com/user-attachments/assets/5dfb589c-ac57-4f73-b52f-d2779c9ece68)

#### Explicit Method

The Explicit Method is a straightforward approach that uses a forward time difference and a central spatial difference to approximate the heat equation. It is computationally simple and easy to implement, making it suitable for small problems. However, it is conditionally stable, meaning that the time step size must satisfy a specific stability criterion (e.g., r=αΔt/(Δx)^2 ≤0.5). Larger time steps can result in instability and inaccuracies, limiting its applicability for complex or long-duration simulations.

#### Implicit Method

The Implicit Method employs a backward time difference and a central spatial difference, making it unconditionally stable. This method is suitable for stiff problems and longer simulations, as stability is guaranteed regardless of the time step size. However, it requires solving a linear system of equations at every time step, which increases computational complexity and runtime compared to the explicit method. It is ideal for scenarios where stability is a higher priority than computational simplicity.

#### Crank-Nicolson Method

The Crank-Nicolson Method combines features of both explicit and implicit methods by averaging the spatial derivative at the current and next time steps. This results in a method that is second-order accurate in both space and time, offering improved precision for smaller time steps. It is both stable and accurate, making it an excellent choice for practical applications. Like the implicit method, it requires solving a system of equations at each step, which can be computationally demanding. However, its balance of stability and accuracy often makes it the preferred method for solving the heat equation in real-world scenarios.

### Runge-Kutta method

This method reduces the PDE to a system of ordinary differential equations (ODEs), which can then be solved. The Runge-Kutta method approximates solutions for differential equations by rewriting the higher order equation into a system of first order equations.

<!DOCTYPE html>
<html lang="en">
<body>
    <h4>Explanation of Runge-Kutta Method for the Heat Equation</h4>
    <h5>Discretization</h5>
    <p>
        Spatial derivatives 
        <math>
            <msup>
                <mfrac>
                    <mrow>&#x2202;u</mrow>
                    <mrow>&#x2202;x</mrow>
                </mfrac>
                <mn>2</mn>
            </msup>
        </math>
        are approximated using a second-order central difference method.<br>
        The partial differential equation (PDE) becomes a system of ordinary differential equations (ODEs) with respect to time.
    </p>
    <h5>Runge-Kutta (RK2)</h5>
    <p>The RK2 method steps are as follows:</p>
    <ul>
        <li>
            <math>
                <msub>k<mn>1</mn></msub>
                = &#x2206;t &sdot; f(&#x1D45D;<sub>n</sub>)
            </math>, 
            where 
            <math>
                f(&#x1D45D;) = &#x1D6FC; 
                <mfrac>
                    <msup>&#x2202;u<mn>2</mn></msup>
                    <msup>&#x2202;x</msup>
                </mfrac>
            </math>.
        </li>
        <li>
            <math>
                <msub>k<mn>2</mn></msub>
                = &#x2206;t &sdot; f(&#x1D45D;<sub>n</sub> + 0.5 &sdot; k<sub>1</sub>)
            </math>
        </li>
        <li>
            <math>
                &#x1D45D;<sub>n+1</sub> = &#x1D45D;<sub>n</sub> + k<sub>2</sub>
            </math>
        </li>
    </ul>
    <h5>Summary</h5>
    <p>
        In this approach, the heat equation is discretized spatially to reduce it to a system of ODEs. The Runge-Kutta method 
        (RK2) is then applied to advance the solution in time, providing a stable and accurate method for solving the heat equation.
    </p>
</body>
</html>

![image](https://github.com/user-attachments/assets/60ec96e9-0d29-4a82-8dc0-c02b28e2805f)

### Monte Carlo methods

Monte Carlo simulation is a mathematical technique that uses random sampling to predict the range of possible outcomes for an uncertain event. 

Here, the Monte Carlo method for solving the heat equation simulates the diffusion of heat by modeling the random movement of particles in a domain. Each particle represents a unit of heat energy, and its random walk approximates the diffusion process described by the heat equation. Particles are allowed to jump to neighboring spatial positions with a probability determined by the thermal diffusivity α, the time step Δt, and the spatial step Δx, ensuring that the simulation satisfies stability conditions. The concentration of particles at each spatial point over time represents the temperature distribution. By averaging the positions and movements of a large number of particles, the Monte Carlo method provides an approximation of the solution to the heat equation. 

![image](https://github.com/user-attachments/assets/72eb9c45-04b4-499d-babf-d80e96c9454a)

![image](https://github.com/user-attachments/assets/3390a8cd-db1b-4e63-8997-52d9b92dd21a)

## Optimization

There are many ways to approach the optimization section. In general, this section relates to how we could enhance the efficiency, accuracy and scalability of the methods used above.
Optimization in the context of the heat equation focuses on improving system performance or achieving specific objectives, such as minimizing energy usage or designing temperature profiles. By incorporating constraints like boundary conditions, energy conservation, or material properties, optimization techniques provide practical solutions for controlling heat flow, solving inverse problems, and enhancing thermal management in various applications. 

### Numerical stability and convergence

We want to make sure  r=αΔt/(Δx)^2 ≤0.5.
This ensures that the numerical solution converges to the correct result without instability by preventing numerical divergence. 

### Computational efficiency

Another approach would be to use more efficient code, more powerful Python libraries (software), or possibly making for larger calculations, notably for the Monte Carlo simulations, using multi-core CPUs or GPUs (hardware). We could also focus more on certain aspects of the heat equation to concentrate the computational power of the computer on something more precise yet more detailed. 

### Energy and error analysis

This reduces errors and improves accuracy, as it monitors the energy decay of the system to verify physical correctness.
By analyzing the discretization error (O(Δx^2,Δt^2)) and adjusting parameters accordingly, we can use extrapolation techniques to reduce numerical errors.
This ensures solutions are both physically meaningful and numerically accurate.

### Lagrange multipliers

<!DOCTYPE html>
<html lang="en">
<body>
    <h4>Lagrange Multipliers in the Heat Equation</h4>
    <h5>1. Constrained Optimization Problems</h5>
    <p>Lagrange multipliers are useful for solving the heat equation while ensuring certain constraints are satisfied. For example:</p>
    <ul>
        <li><strong>Energy Conservation:</strong> Enforce that the total energy, such as the integral of the temperature over the domain, remains constant.</li>
        <li><strong>Boundary Constraints:</strong> Impose specific temperature values or heat fluxes at the boundaries beyond standard Dirichlet or Neumann conditions.</li>
    </ul>
    <p>The optimization problem is often written as: minimize J(u), subject to g(u) = 0, where:</p>
    <ul>
        <li>J(u): The cost or objective function, for example, the integral of temperature differences over time.</li>
        <li>g(u): A constraint function, such as energy conservation or fixed boundary flux.</li>
    </ul>
    <p>Using Lagrange multipliers, the problem is reformulated as:</p>
    <p>
        $$L(u, \lambda) = J(u) + \lambda \cdot g(u)$$
    </p>
    <p>The solution satisfies:</p>
    <p>
        $$\frac{\partial L}{\partial u} = 0$$
    </p>
    <p>and</p>
    <p>
        $$\frac{\partial L}{\partial \lambda} = 0$$
    </p>
    <h5>2. Inverse Problems</h5>
    <p>Lagrange multipliers can help determine unknown parameters, such as thermal diffusivity, or initial conditions that lead to a desired temperature distribution. For example:</p>
    <p>The constraint ensures that the solution satisfies the heat equation:</p>
    <p>
        $$L(u, \lambda) = \int (u_t - \alpha u_{xx})^2 \, dx + \int \lambda \cdot (u_t - \alpha u_{xx}) \, dx$$
    </p>
    <p>Minimizing L ensures that u satisfies the heat equation while fitting observed data.</p>
    <h5>3. Optimization of Control Inputs</h5>
    <p>In practical applications, you might control heat sources to achieve a specific temperature profile. For example:</p>
    <p>Minimize the energy used to maintain a target temperature distribution:</p>
    <p>
        $$J(u, q) = \int (u(x, t) - u_{\text{target}}(x, t))^2 \, dx + \int q(x, t)^2 \, dx$$
    </p>
    <p>Subject to:</p>
    <p>
        $$u_t = \alpha u_{xx} + q(x, t)$$
    </p>
    <p>Here, q(x, t) represents the heat input, and u_target is the desired temperature profile. Lagrange multipliers ensure the solution respects the heat equation dynamics.</p>
    <h5>Applications of Lagrange Multipliers</h5>
    <ul>
        <li><strong>Thermal Management:</strong> Design heat source placement or control inputs for devices.</li>
        <li><strong>Material Optimization:</strong> Adjust properties like thermal conductivity or diffusivity to achieve desired temperature behavior.</li>
        <li><strong>Boundary Condition Design:</strong> Solve for boundary fluxes or interface conditions to meet operational requirements.</li>
    </ul>
    <h5>Advantages of Using Lagrange Multipliers</h5>
    <ul>
        <li>Provide a systematic way to handle constraints during optimization.</li>
        <li>Help incorporate physical laws, such as conservation of energy, directly into the problem formulation.</li>
        <li>Enable solving inverse or control problems in thermal systems.</li>
    </ul>
    <p>By combining Lagrange multipliers with computational methods, such as finite difference or finite element methods, you can address complex, constrained problems in heat transfer efficiently.</p>
</body>
</html>

Here are examples of basic Python code using Lagrange multipliers: [Examples](https://github.com/aa6dcc/Heat-Equation/tree/main) 

This project examined the heat equation, different ways of approaching it, solving it, estimating it and visualizing it.
Its purpose is also to give the reader a better understanding of the topic - I hope it was successful in this regard!

## Bibliography

• *Heat equation* [Wikipedia]: [https://en.wikipedia.org/wiki/Heat_equation](https://en.wikipedia.org/wiki/Heat_equation)

• *Deriving the heat equation: A Parabolic Partial Differential Equation for Heat Energy Conservation* by Steve Brunton [Youtube]: [https://www.youtube.com/watch?v=9d8PwnKVA-U&ab_channel=SteveBrunton](https://www.youtube.com/watch?v=9d8PwnKVA-U&ab_channel=SteveBrunton)

• *Solving the heat equation* by 3Blue1Brown [Youtube]: [https://www.youtube.com/watch?v=ToIXSwZ1pJU&ab_channel=3Blue1Brown](https://www.youtube.com/watch?v=ToIXSwZ1pJU&ab_channel=3Blue1Brown)

• *Green's function for the heat equation* [Mathematics Stack Exchange]: [https://math.stackexchange.com/questions/2588931/green-s-function-for-the-heat-equation](https://math.stackexchange.com/questions/2588931/green-s-function-for-the-heat-equation)

• *The tangled history of the second law of thermodynamics* [Stephen Wolfram]: [https://writings.stephenwolfram.com/2023/01/how-did-we-get-here-the-tangled-history-of-the-second-law-of-thermodynamics/](https://writings.stephenwolfram.com/2023/01/how-did-we-get-here-the-tangled-history-of-the-second-law-of-thermodynamics/)

• *Mathematics for Engineers 1 (ELEC40012) (chapter 11, p.50)* [Dan Nucinkis, Imperial College London]
