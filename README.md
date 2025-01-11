# Solving the heat equation

The heat equation is a partial differential equation.

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

At steady state, the heat equation simplifies down to Laplace's equation, a second order partial differential equation:

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

From there, we have:

<svg viewBox="0 0 400 200" xmlns="http://www.w3.org/2000/svg">
  <!-- Left cylinder -->
  <circle cx="100" cy="100" r="30" stroke="white" fill="none" stroke-width="2"/>
  <!-- Right cylinder -->
  <circle cx="200" cy="100" r="30" stroke="white" fill="none" stroke-width="2" stroke-dasharray="5,5"/>
  
  <!-- Arrow -->
  <path d="M130 100 L170 100" stroke="white" fill="none" stroke-width="2"/>
  <path d="M170 100 L160 95 L160 105 Z" fill="white"/>
  
  <!-- Text labels -->
  <text x="70" y="100" fill="white" font-family="Arial" font-size="14">q(x,t)</text>
  <text x="210" y="100" fill="white" font-family="Arial" font-size="14">q(x+Δx,t)</text>
</svg>

<div class="equation">
  <p>q(x,t) = heat flux from left to right (thermal energy/area)</p>
  
  <p>
    <sup>∂q</sup>⁄<sub>∂x</sub> = lim<sub>Δx→0</sub>
    <span class="fraction">
      q(x,t) - q(x+Δx,t)
      <hr>
      Δx
    </span>
  </p>

  <p>q(x,t) = -k<span class="fraction">∂u<hr>∂x</span></p>
</div>

.equation {
  font-family: "Times New Roman", serif;
  font-size: 16px;
}

.fraction {
  display: inline-block;
  vertical-align: middle;
  text-align: center;
}

.fraction hr {
  border: none;
  border-top: 1px solid black;
  margin: 2px 0;
}

Fourier worked on the heat equation (in fact he invented the Fourier transform as a tool to solve the heat equation) and derived Fourier's law of heat conduction, which we 

## Generalizing the solution technique

## Computational methods

## Optimization


