# 2D Ising Model

Monte Carlo simulations of the two-dimensional Ising model using Metropolis and Wolff algorithms.

---

## Overview

This repository contains a group project for a computational physics course, focusing on numerical simulations of the two-dimensional ferromagnetic Ising model.

We study equilibrium properties, autocorrelation effects, critical slowing down, finite-size behavior, and non-equilibrium dynamics across the thermal phase transition by implementing and comparing two Monte Carlo algorithms:

- Metropolis single-spin update
- Wolff cluster update

---

## Core Simulation Code

### Metropolis Algorithm

- [`Metropolis.py`](./Metropolis.py)  
  Core Metropolis Monte Carlo implementation.

- [`Metropolis.ipynb`](./Metropolis.ipynb)  
  Notebook version of Metropolis simulations.

- [`drawing.ipynb`](./drawing.ipynb)  
  Analysis of Metropolis results: autocorrelation time, temperature dependence, finite-size effects.

---

### Wolff Cluster Algorithm

- [`wolff_ising.py`](./wolff_ising.py)  
  Wolff cluster algorithm implementation, including equilibrium simulations and linear annealing.

- [`drawing_wolff.py`](./drawing_wolff.py)  
  Analysis and visualization for Wolff simulations.

---

## Report and LaTeX Files

All report-related files are stored in the [`Report Latex/`](./Report%20Latex) directory.

### LaTeX Source and PDF

- [`final1.0.tex`](./Report%20Latex/final1.0.tex)
- [`final1.0.pdf`](./Report%20Latex/final1.0.pdf)
- [`final1.0Notes.bib`](./Report%20Latex/final1.0Notes.bib)
- [`reference.bib`](./Report%20Latex/reference.bib)
- [`sample.bib`](./Report%20Latex/sample.bib)

---

### Figures Used in the Report

- [`fig1.png`](./Report%20Latex/fig1.png)
- [`fig2.png`](./Report%20Latex/fig2.png)
- [`fig3.png`](./Report%20Latex/fig3.png)
- [`fig4.png`](./Report%20Latex/fig4.png)
- [`fig5.png`](./Report%20Latex/fig5.png)
- [`fig6.png`](./Report%20Latex/fig6.png)
- [`fig7.png`](./Report%20Latex/fig7.png)
- [`fig8.png`](./Report%20Latex/fig8.png)
- [`fig9.png`](./Report%20Latex/fig9.png)
- [`fig10.png`](./Report%20Latex/fig10.png)

---

### Data and Analysis Plots

- [`E_theory.csv`](./Report%20Latex/E_theory.csv)
- [`tau_vs_T.png`](./Report%20Latex/tau_vs_T.png)
- [`tau_vs_L.png`](./Report%20Latex/tau_vs_L.png)
- [`autt_vs_tem.jpg`](./Report%20Latex/autt_vs_tem.jpg)
- [`autt_vs_tem_2.png`](./Report%20Latex/autt_vs_tem_2.png)
- [`autt_vs_leng.jpg`](./Report%20Latex/autt_vs_leng.jpg)
- [`autt_vs_leng_2.png`](./Report%20Latex/autt_vs_leng_2.png)
- [`collapse_M.png`](./Report%20Latex/collapse_M.png)
- [`collapse_M_single.png`](./Report%20Latex/collapse_M_single.png)
- [`raw_M_vs_T.png`](./Report%20Latex/raw_M_vs_T.png)
- [`raw_M_vs_T_single.png`](./Report%20Latex/raw_M_vs_T_single.png)
- [`fint_sca.jpg`](./Report%20Latex/fint_sca.jpg)
- [`wolff.png`](./Report%20Latex/wolff.png)

---

## Notes

- All files listed above can be opened directly by clicking the links.
- For source files (`.py`, `.tex`), use **“Raw”** on GitHub to download.
- For figures and PDFs, clicking the link opens or downloads the file directly.

