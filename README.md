# Numerical Methods for Nonlinear PDEs

This repository includes MATLAB scripts for solving various types of nonlinear partial differential equations (PDEs). The focus is on conservation laws, Hamilton-Jacobi equations, and optimal control problems.

## Table of Contents

1. [Overview](#overview)
2. [Conservation Laws](#conservation-laws)
3. [Hamilton-Jacobi Equations](#hamilton-jacobi-equations)
4. [Control Theory](#control-theory)

## Overview

The repository contains code for several numerical methods to address different types of PDEs. Each type has its unique set of methods implemented, providing a versatile toolkit for researchers and practitioners.

## Conservation Laws
Impements schemes for the transport equation.
All conservation laws are demonstrated in `app_conservation_laws_exported.m`, which serves as the main script for testing these schemes. The following schemes are included:
- Upwind
- Upwind correction Enquit-Auscher for convex and concave flux case 
- Lax-Friedrichs 
- MacCormack 
- Godunov for convex and concave case
- Lax-Wendroff.

Also, the script `entropy.m` describes the entropy solution of the transport equation when there is a shock.
## Hamilton-Jacobi Equations

### Semilagrangian Approach
- **Stationary Case:** Implemented in `semilagrangian_stat.m`.
- **Evolutive Case:** Found in `semilagrangian_evo.m`.
The algorithms are tested on `Burger's` and `Eikonal` equations.

### Parallelism with convervation laws
Parallel algorithms for Hamilton-Jacobi equations are provided in `hj_evo.m` and `hj_ptfix` for evolutive and stationary case respectively.

## Control Theory
The approach to solve control problems is derived from Bellman **Dynami Programming Principle (DPP)** discretized in Semi-Lagrangian schemes applied to:
- Finite horizion problem
- Minimun time problem
- Infinite horizion problem 
Optimal control problems include scripts to derive and run simulations. Finite horizon examples are handled with `finite_horizon.m` and infinite horizon scenarios with `infinite_horizon.m`.
A couple of examples are aviable. To run the simulation you need to generate the V-function matrix data with the relative `...Data.m` and the run the .m file.
