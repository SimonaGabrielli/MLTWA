# MLTWALab: Multi-Layer Multiple Time Window Analysis & Bias Correction

MLTWALab is a MATLAB-based toolbox designed to model seismic wave propagation through complex multi-layered structures
using Monte Carlo simulations. It includes a dedicated calibration and inversion framework to correct the geometric
and structural biases that occur when interpreting multi-layer structural data through a simplified
homogeneous half-space model.

This methodology aligns with the classical multiple time window analysis (MLTWA) framework and incorporates advanced correction techniques
inspired by Del Pezzo & Bianco (2010a) and Akinci et al. (2020).

* Project Overview

When evaluating seismic attenuation, scattering (B0), and inverse total leakage lengths (L_e^-1),
assuming a completely homogeneous half-space can introduce systematic errors if the true medium contains a multi-layered velocity or scattering profile (e.g., a crust-mantle boundary).

This project provides a complete laboratory environment to:
* Calculate B0 and Le^-1 for an half-space model, obtaining also values of Qtot, Qi and Qs.
* Simulate by running millions of particles through multi-layered scattering structures via Monte Carlo methods.
* Calibrate by mapping the relationship between the true input parameters and the retrieved apparent parameters.

------------

* Workflow & Theory

The toolbox utilizes a two-stage operational framework:

1. A discrete meshgrid of true structural parameters is defined.
For each point, a Monte Carlo simulation generates raw energy envelopes.
A traditional Grid Search algorithm finds the best-fitting apparent parameters by minimizing the residuals between the simulated envelopes and the theoretical half-space model curves (`computeEnergy3`).

2. Using the forward mapping dataset, the script builds a data kernel matrix G based on the apparent formulation.
A least-squares inversion yields 5 unique polynomial coefficients, allowing direct correction of any experimental observation.

---

* Repository Structure

1. MLTWA_Code.mlx: Implementation of the Multiple Time Window Analysis window extraction routines.
2. Scattering_simulation.mlx: Live Script handling the baseline setup and execution of the core layer-scattering Monte Carlo tracker.
3. Correction.mlx: runs the forward grid calibration loop, executes the least-squares inversion, applies the bias correction to the target station.

---

* Prerequisites
MATLAB (R2020a or later recommended)
Parallel Computing Toolbox (for accelerating the Monte Carlo `parfor` loops) - Not fundamental but highly suggested
The user can select to not use parallelization, or the number of cores for it, by selecting 0 (no parallelization), N (number of cores), auto (uses all available cores (maximum speed).

In cases where a local MATLAB license is not available, the package can also be executed using the MATLAB Online platform.
In this case, consider that the Parallel Computing Toolbox is not available, so you need to set the number of workers to 0.

* Running the Calibration and Correction

1. Open MATLAB and navigate to the repository directory.
2. Provide your experimental dataset as input: the code requires the integrated energy values computed within the three distinct time windows for each station-distance combination, plus the coda normalization window (more details in `MLTWA_Code.mlx`).
3. Run the main execution script `MLTWA_Code.mlx`.
4. Once obtained the values of B0 and Le, you can run the code `Scattering_simulation.mlx` to obtain the first energy corrections for a depth dependent simulations.
5. In case you have the MLTWA results, you can also run the code `Correction.mlx`, to obtain a further corrections with polynomials.
6. Now, you can run again the `Scattering_simulation.mlx` code to visualize the differences in energies with the polynomial corrections.

As a default, the values for the smoothing in the MLTWA_Code.mlx are 20 (windows size) and 30 (number of bins). In case you want to change these values, you can modify them
in the function `run_MLTWA_automation.m`, line 25 (`prepare_smoothed_data`).
---

* Outputs & Visualizations

The script automatically exports publication-ready graphics illustrating the spatial dynamics of the geometric bias.

---

* References
Akinci, et al. (2020) *Intrinsic and scattering seismic wave attenuation in the Central Apennines (Italy). Physics of the Earth and Planetary Interiors, 303, 106498.*
Del Pezzo, E., & Bianco, F. (2010a) *MathLTWA: Multiple lapse time window analysis using Wolfram Mathematica 7. Computers & geosciences, 36(10), 1388-1392.*
Del Pezzo, E., & Bianco, F. (2010b) *Two-layer earth model corrections to the MLTWA estimates of intrinsic-and scattering-attenuation obtained in a uniform half-space. Geophysical Journal International, 182(2), 949-955*
