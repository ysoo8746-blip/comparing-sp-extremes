# Comparing Spatial Extremes (JCGS Reproducibility Repository)

Code and data to reproduce the analyses in:

**Multiple Testing for Spatial Extremes with Application to Reanalysis Data Evaluation**  
Sooin Yun, Bo Li, Xianyang Zhang (2026)

## Repository structure

- `Rcodes/functions/`  
  Reusable functions used across simulation and data analysis.

- `Rcodes/Simulation/`  
  Scripts to generate simulated data, run competing methods, summarize results, and produce simulation figures/tables.

- `Rcodes/Climate_model_comparison/`  
  Scripts for the real-data analysis comparing reanalysis and observational products (e.g., winter precipitation extremes), including hypothesis testing and plots.

- `Data/Simdata/`  
  Saved simulation datasets/objects used to generate figures/tables.

- `Data/Realdata/`  
  Analysis-ready real data used in the manuscript (or instructions to obtain it if redistribution is restricted).

- `output/`  
  Generated figures/tables (may be re-created by running the pipeline).

## Requirements

- R (>= 4.x recommended)
- Suggested: `renv` for package reproducibility

If using renv:
```r
install.packages("renv")
renv::restore()
