1D Jominy End-Quench Cooling Simulation and Visualization
Overview

This project implements a simplified one-dimensional heat-transfer model of a Jominy end-quench test and provides both numerical cooling curves and a spatial–temporal visualization of the temperature field along the specimen. The workflow is designed to illustrate the fundamental thermal behavior of a quenched steel bar under idealized conditions and to generate data suitable for qualitative analysis and educational purposes.

The model solves transient heat conduction in a cylindrical bar with convection applied at the quenched end and an insulated boundary at the opposite end. Material properties and heat-transfer coefficients are assumed constant, representing an effective, time-averaged cooling condition.

Physical Model

The governing equation is one-dimensional transient heat conduction:

𝜌
𝑐
𝑝
∂
𝑇
∂
𝑡
=
𝑘
∂
2
𝑇
∂
𝑥
2
ρc
p
	​

∂t
∂T
	​

=k
∂x
2
∂
2
T
	​


where 
𝑇
(
𝑥
,
𝑡
)
T(x,t) is temperature, 
𝜌
ρ is density, 
𝑐
𝑝
c
p
	​

 specific heat capacity, and 
𝑘
k thermal conductivity.

Boundary conditions

Quenched end (x = 0):
Convection to a coolant:

−
𝑘
∂
𝑇
∂
𝑥
=
ℎ
(
𝑇
−
𝑇
∞
)
−k
∂x
∂T
	​

=h(T−T
∞
	​

)

Far end (x = L):
Insulated boundary:

∂
𝑇
∂
𝑥
=
0
∂x
∂T
	​

=0

The initial condition is a uniform high temperature representing austenitization.

Numerical Method

Spatial discretization: second-order central finite differences

Time integration: implicit Backward Euler scheme (unconditionally stable)

Linear system: tridiagonal, solved using the Thomas algorithm

The implicit formulation allows stable simulation of rapid cooling without restrictive time-step constraints.

Project Structure and Workflow
1. T_Hist files.py

Runs the full transient heat-transfer simulation and stores the complete temperature history:

Generates spatial grid (x.npy)

Generates time vector (t.npy)

Stores full temperature evolution (T_hist.npy)

These files serve as input for post-processing and animation.

2. Cooling_Curves.py

Computes and plots cooling curves at selected distances from the quenched end:

Extracts temperature histories at fixed spatial intervals

Produces logarithmic time-scale cooling curves

Highlights representative positions along the bar

This script is used for quantitative inspection of cooling behavior versus distance.

3. Cooling animation.py

Creates a spatial–temporal visualization of the cooling process:

Loads x.npy, t.npy, and T_hist.npy

Displays temperature as a 1D heatmap along the bar

Saves an MP4 animation showing the progression of cooling over time

The animation is intended for visualization and presentation rather than numerical analysis.

Assumptions and Limitations

One-dimensional heat flow

Constant material properties

Constant effective convection coefficient

No phase transformation or latent heat effects

As such, the model does not reproduce real Jominy hardness curves quantitatively, but captures the dominant thermal trends and spatial gradients.

Intended Use

This project is suitable for:

Educational demonstrations of transient heat conduction

Qualitative analysis of quenching behavior

Visualization of thermal gradients in end-quenched specimens

Baseline studies prior to more advanced thermo-metallurgical modeling

Requirements

Python 3.x

NumPy

Matplotlib

FFmpeg (for MP4 animation export)

How to Run

Run T_Hist files.py to generate simulation data

Run Cooling_Curves.py to plot cooling curves

Run Cooling animation.py to generate and save the animation