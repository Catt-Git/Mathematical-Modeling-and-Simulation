# Mathematical Modeling of HIV Infection Dynamics

**Authors:** Alberto Catalano, Anna Salsano  
**Course:** Mathematical Modeling and Simulation (2025–2026)

## 📌 Overview

This repository contains MATLAB implementations of **four mathematical models** describing the **within-host dynamics of HIV infection**.  
The models are based on systems of **ordinary differential equations (ODEs)** and progressively increase in biological complexity, capturing viral replication, immune responses, latency, and treatment effects.

## 📂 Repository Structure

| File | Description |
|------|------------|
| `basic_hiv_model.m` | Classical Nowak & May model with uninfected cells (T), infected cells (I), and virus (V). |
| `cell_to_cell_hiv_model.m` | Adds latent reservoirs (L) and cell-to-cell transmission. |
| `Immune_system_hiv_model_def.m` | Includes immune responses: CTLs (E) and antibodies (Z). |
| `complex_hiv_model_def.m` | High-fidelity model (16 ODEs) with thymic output, macrophages, and viral tropism switching. |
| `complex_hiv_model_def_ART.m` | Complex model configured to simulate antiretroviral therapy (ART). |
| `Report_Catalano_Salsano.pdf` | Report of the project with the full explenation of everything |
| `Presentation_Catalano_Salsano.ppt` | Powerpoint presentation of the project |

## 🔬 Models Summary

- **Basic model:** Reproduces acute infection and chronic viral set-point.
- **Latency model:** Explains viral persistence via long-lived latent reservoirs.
- **Immune model:** Demonstrates immune control through cellular and humoral responses.
- **Complex model:** Captures triphasic progression (acute → latency → AIDS) and ART-driven immune reconstitution.

## 🚀 Usage

### Requirements
- MATLAB (R2020a or later recommended)
- No additional toolboxes required (`ode45`, `ode15s`)

### Run a Simulation
```matlab
complex_hiv_model_def_ART
