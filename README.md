# Pendulums.jl

A Julia simulation and visualization project for single and double pendulums, illustrating their classical and chaotic dynamics.

<p align="center">
  <img src="/images/double_pendulum_initial_condition.gif" alt="Double pendulums" width="400"/>
  <img src="/images/double_pendulum.linear_approximation1.png" alt="Comparison with linear approximation" width="400"/>
</p>


---

## Overview

This repository provides a collection of numerical simulations and animations for:

- Simple pendulum (undamped and damped)
- Forced oscillation and resonance
- Double pendulum in chaotic and linear regimes
- Comparison between numerical and analytical solutions (normal modes)

The core solvers are implemented from scratch, including a 4th order Runge-Kutta and a symplectic Euler integrator.

---

## Project Structure

```
Pendulums.jl/
├── src/
│   ├── ODESolvers.jl           # Runge-Kutta and symplectic integrator
│   ├── Pendulum.jl             # Dynamics of single and double pendulums
│   ├── examples1.jl            # Basic single pendulum animation
│   ├── examples2.jl            # Damping and forced oscillation analysis
│   ├── examples3.jl            # Double pendulum and chaos visualization
│   ├── examples4.jl            # Comparison: numerical vs analytical solution
├── images/                    # GIFs and figures
└── README.md                  # This file
```

---

## Features & Demos

- GIF animations with Plots.jl
- Custom integrator framework with method dispatch
- Chaotic trajectory tracing for the double pendulum
- Phase space & FFT spectrum for harmonic analysis
- Analytical comparison of small-angle normal modes

---

## Getting Started

```julia
using Pkg
Pkg.add("Plots")
Pkg.add("FFMPEG")  # For GIF generation

include("ODESolvers.jl")
include("Pendulum.jl")
include("examples1.jl")  # or any other example
```

---

## Physics Concepts Illustrated

- Harmonic motion
- Energy dissipation and resonance
- Nonlinear dynamics and chaos
- Hamiltonian systems and symplectic integration
- Normal mode analysis

---

## Background of Double Pendulum

The Euler--Lagrange equations for the double pendulum are given as

$$(m_1+m_2){L_1}^2\ddot \theta_1+m_2L_1L_2\left[\ddot\theta_2\cos(\theta_1-\theta_2)+{\dot\theta_2}^2\sin(\theta_1-\theta_2)\right]=-(m_1+m_2)gL_2 \sin\theta_1,$$
$$m_2{L_2}^2\ddot\theta_2+m_2L_1L_2\left[\ddot\theta_1\cos(\theta_1-\theta_2)-{\dot\theta_1}^2\sin(\theta_1-\theta_2)\right]=-m_2gL_2\sin\theta_2.$$

Under the conditions of $m_1=m_2=m$ and $L_1=L_2=L$, the equations are simplified as follows:

$$\ddot\theta_1=\cfrac{-{\dot\theta_1}^2\sin\Delta\cos\Delta-{\dot\theta_2}^2\sin\Delta+\cfrac{g}{L}(\sin\theta_2\cos\Delta-2\sin\theta_1)}{2-\cos\Delta},$$
$$\ddot\theta_2=\cfrac{{\dot\theta_2}^2\cos\Delta\sin\Delta+2{\dot\theta}^2\sin\Delta+\cfrac{2g}{L}(\sin\theta_1\cos\Delta-\sin\theta_2)
}{2-\cos\Delta},$$

where $\Delta:=\theta_1-\theta_2$.

<p align="center">
  <img src="/images/double_pendulum_trace.gif" alt="Double pendulums with trace" width="400"/>
</p>

In particular when the amplitudes $\theta_1$ and $\theta_2$ are small, they are linearized:

$$\ddot\theta_1=-\dfrac{2g}{L}\theta_1+\dfrac{g}{L}\theta_2,$$
$$\ddot\theta_2=\dfrac{2g}{L}\theta_1-\dfrac{2g}{L}\theta_2.$$

For simplicity, let us consider the initial condition of $\dot\theta_1(t=0)=\dot\theta_2(t=0)=0$. Then the linearized equation explicitly solved and the solution is given by

$$\theta_1(t)=\dfrac{1}{2}\left(\theta_1(t=0)-\dfrac{\theta_2(t=0)}{\sqrt2}\right)\cos\Omega_+ t+\dfrac{1}{2}\left(\theta_1(t=0)+\dfrac{\theta_2(t=0)}{\sqrt2}\right)\cos\Omega_-t$$

---

## Notes

- Created by Ushihara as a side project for learning and exploration.
- A part of this README is generated by using ChatGPT, based on my original code implementation.
