# README

This repository contains MATLAB code to reproduce and explore the models presented in the article A Theory for Cyclic Shifts between Alternative States in Shallow Lakes (Van Nes, Rip & Scheffer). The project investigates the two alternative states typical of shallow lakes — clear water with macrophytes and turbid water dominated by phytoplankton — as well as the mechanisms (fast and slow feedbacks) that may trigger regime shifts or multi-year cycles between these states.

Repository objectives:

Implementation of two simple dynamic models inspired by the article:

Mechanism 1 — Phosphorus retention (effect of macrophytes on the partitioning of phosphorus between water and sediments).

Mechanism 2 — Anoxic phosphorus release linked to the accumulation of organic matter / SOD.

For both mechanisms, the study focuses on the model proposed in the article, in particular:

- Solving the equations numerically and reproducing the temporal evolution of the system

- Identifying and analysing the equilibrium points and studying their stability

- Performing sensitivity analyses and, if desired, bifurcation analyses (Hopf / infinite-period)

## Model 1

- modelisation.m: Computation of the non-zero equilibrium, plotting of the nullclines and local simulation of the dynamics using Euler
- jacobienne.m: computation and plotting of the determinant and trace of the Jacobian for several values of Pw

## Model 2

- constants.m: script containing the parameters used in the simulations

The second mechanism allows the reduction of the three-dimensional system (V, P, SOD) to a two-dimensional system (V, P). The scripts are therefore divided into two folders corresponding to each formulation.

### 3D System

- three_ode.m: implementation of the system of three time-dependent differential equations describing the evolution of the lake
- three_ode_res.m: numerical solution of the 3D system and plotting of the results (V and P versus time)

### 2D System

- two_ode.m: implementation of the system of two time-dependent differential equations describing the evolution of the reduced system
- two_ode_res.m: numerical solution of the 2D system and plotting of the results (V and P versus time)
- nullclines_and_eq.m: graphical determination and visualisation of the equilibrium points
- stability_eq.m: stability analysis of the equilibria and computation of the Jacobian