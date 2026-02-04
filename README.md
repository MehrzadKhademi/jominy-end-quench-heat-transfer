1D Jominy End-Quench Heat-Transfer Simulation and Visualization
Overview

This project implements a simplified one-dimensional transient heat-transfer model of a Jominy end-quench test. The workflow computes cooling curves along a steel bar and produces a spatial–temporal animation of the temperature field during quenching. The model is intended for qualitative analysis, visualization, and educational use rather than quantitative prediction of hardness.

Transient heat conduction is solved in a bar with convection applied at the quenched end and an insulated boundary at the opposite end. Material properties and the heat-transfer coefficient are assumed constant, representing an effective time-averaged cooling condition.

Physical Model

Heat transfer is governed by the one-dimensional transient conduction equation:

rho * c_p * dT/dt = k * d²T/dx²

where T(x,t) is temperature, rho is density, c_p is specific heat capacity, and k is thermal conductivity.

Boundary and Initial Conditions

At the quenched end (x = 0), heat is removed by convection to a coolant:

k * dT/dx = h * (T − T_inf)

where h is the effective convection coefficient and T_inf is the coolant temperature.

At the far end of the bar (x = L), the boundary is insulated:

dT/dx = 0

The initial condition is a uniform temperature along the bar:

T(x,0) = T_init

representing an austenitized specimen prior to quenching.

Numerical Method

The governing equation is discretized using second-order central finite differences in space and an implicit Backward Euler scheme in time. The implicit formulation is unconditionally stable and allows simulation of rapid cooling without restrictive time-step limitations.

At each time step, the discretized equations form a tridiagonal linear system, which is solved efficiently using the Thomas algorithm.

Project Workflow

The project consists of three Python scripts executed in sequence.

1. T_Hist files.py

This script performs the full transient heat-transfer simulation. It defines the spatial grid and time discretization, applies boundary and initial conditions, and computes the temperature evolution along the bar. The complete temperature history is saved to disk as NumPy arrays: spatial coordinates, time vector, and temperature history.

2. Cooling_Curves.py

This script post-processes the simulated temperature history to extract cooling curves at selected distances from the quenched end. Temperature–time curves are plotted on a logarithmic time scale, and representative locations along the bar are highlighted to illustrate the spatial variation in cooling behavior.

3. Cooling animation.py

This script visualizes the transient cooling process. It loads the saved temperature history and generates an animation showing temperature as a function of position and time along the bar. The animation is exported as an MP4 file for presentation and demonstration purposes.

Assumptions and Limitations

The model is intentionally simplified and includes the following assumptions:

One-dimensional heat flow

Constant material properties

Constant effective convection coefficient

No phase transformations or latent heat effects

As a result, the model does not reproduce experimental Jominy hardness curves quantitatively but captures the dominant thermal trends and spatial temperature gradients.

Intended Use

This project is suitable for:

Educational demonstrations of transient heat conduction

Qualitative analysis of end-quench cooling behavior

Visualization of thermal gradients in steel specimens

Baseline studies prior to more advanced thermo-metallurgical modeling

Requirements

Python 3.x

NumPy

Matplotlib

FFmpeg (required for saving MP4 animations)

How to Run

Run T_Hist files.py to generate the temperature history

Run Cooling_Curves.py to plot cooling curves

Run Cooling animation.py to generate and save the animation