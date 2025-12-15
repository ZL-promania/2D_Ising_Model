# 2D Ising Model

Monte Carlo simulations of the two-dimensional Ising model with Metropolis and Wolff algorithms.

---

## Overview

This repository contains numerical simulations of the two-dimensional ferromagnetic Ising model, carried out as a group project for a computational physics course.

We study the equilibrium and dynamical properties of the model across the thermal phase transition by implementing and comparing two Monte Carlo update schemes:

- **Metropolis single-spin updates**
- **Wolff cluster updates**

The project focuses on physical observables, autocorrelation effects, critical slowing down, and finite-size behavior near the critical temperature.

---

## Physics Background

The two-dimensional Ising model consists of spins \( s_i = \pm 1 \) on a square lattice with nearest-neighbor interactions,

\[
H = -J \sum_{\langle i,j \rangle} s_i s_j ,
\]

and exhibits a continuous phase transition at the exact critical temperature

\[
k_B T_c = \frac{2J}{\ln(1+\sqrt{2})}.
\]

Because of its exact solution and well-understood critical behavior, the 2D Ising model serves as a benchmark system for testing Monte Carlo algorithms and studying critical phenomena.

---

## Monte Carlo Methods

### Metropolis Algorithm

The Metropolis algorithm performs local single-spin flips that satisfy detailed balance with respect to the Boltzmann distribution.  
While it correctly reproduces equilibrium thermodynamic quantities, it suffers from **critical slowing down** near the phase transition, leading to long autocorrelation times for global observables.

---

### Wolff Cluster Algorithm

The Wolff algorithm performs non-local cluster updates based on the Fortuin–Kasteleyn representation.  
By flipping correlated spin clusters in a single update, it significantly reduces autocorrelation times near the critical temperature and improves sampling efficiency.

---

## Measured Observables

From Monte Carlo time series generated at equilibrium, we compute:

- Energy per spin  
- Magnetization (absolute value)  
- Susceptibility  
- Binder cumulant  
- Integrated autocorrelation time  

Autocorrelation times are used as the primary diagnostic for sampling efficiency and critical slowing down.

---

## Finite-Size Effects

Simulations are performed on square lattices of various linear sizes \( L \) with periodic boundary conditions.

Finite-size effects are analyzed by comparing temperature-dependent magnetization curves for different \( L \), as well as through finite-size scaling and data collapse using the exact critical exponents of the 2D Ising universality class.

---

## Linear Annealing

In addition to equilibrium simulations at fixed temperature, a linear annealing protocol is implemented using the Wolff cluster algorithm.

The inverse temperature is varied linearly during the simulation to study non-equilibrium relaxation and the effect of cooling rate on energy and magnetization trajectories.

---

## Repository Structure

```text
2D_Ising_Model/
├── wolff_ising.py        # Wolff cluster Monte Carlo implementation
├── metropolis_ising.py  # Metropolis Monte Carlo implementation
├── analysis/            # Autocorrelation, finite-size scaling, plotting scripts
├── README.md
