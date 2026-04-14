# 15O Lifetime Analysis — AGATA @ LNL

Data analysis pipeline for measuring nuclear state lifetimes in **¹⁵O** and **¹⁸F** using the Doppler-shift attenuation method (DSAM) with the AGATA γ-ray tracking array at Laboratori Nazionali di Legnaro (LNL).

---

## Experiment Overview

A **¹⁶O beam at 50 MeV** was directed onto a **³He-implanted gold target**. Two reaction channels are simultaneously active:

| Channel | Reaction | Recoil |
|---------|----------|--------|
| 1 | ¹⁶O(³He, p) | ¹⁸F |
| 2 | ¹⁶O(³He, α) | ¹⁵O |

Recoiling nuclei are produced in excited states inside the gold stopper. As they slow down, they emit γ-rays detected by **AGATA**. Light ejectiles (protons or α particles) are detected by **SAURON**, which provides the reaction kinematics event-by-event.

The **DSAM observable** is:

```
Δβ = β_reaction − β_emission
```

where `β_reaction` is the recoil velocity reconstructed from the ejectile kinematics and `β_emission` is the velocity at the moment of γ emission inferred from the Doppler shift. A longer lifetime means more slowing-down in the stopper, resulting in a larger Δβ. The linear relationship between Δβ and lifetime is calibrated using GEANT4 simulations and literature values of known ¹⁸F state lifetimes.

---

## Repository Structure

```
.
├── data/
│   └── 15O_Lifetime_AGATA.npy       # Full dataset (event-by-event, keyed by γ energy)
│
├── results/
│   ├── delta_beta.txt               # Measured Δβ values per transition (input to Bayesian fit)
│   ├── simulation.txt               # Simulated (τ, mean Δβ) calibration curve
│   └── lifetime_samples.txt         # MCMC posterior samples for each measured lifetime
│
├── plots/                           # All output figures (one set per γ transition)
│
├── pykinematic/
│   ├── kinematics.py                # Relativistic kinematics library
│   ├── __init__.py
│   └── masstable.txt                # AME nuclear mass table
│
├── doppler.py                       # Relativistic Doppler shift utilities
│
├── Analyze_<E>.ipynb                # Per-transition analysis notebooks (see below)
├── Bayesian_Fit.ipynb               # Global Bayesian fit → lifetimes
├── Plot_Doppler.ipynb               # Diagnostic Doppler-shift plots
└── Plot_Simulation.ipynb            # Simulation calibration curve plots
```

---

## Analysis Notebooks

### `Analyze_<E>.ipynb` — Per-transition analysis

One notebook exists for each γ-ray transition energy `E` (in keV). The analysed transitions are:

| Energy (keV) | Nucleus | Ex (MeV) |
|-------------|---------|----------|
| 777  | ¹⁸F | 3.839 |
| 1041 | ¹⁸F | 3.724 |
| 1298 | ¹⁸F | 4.360 |
| 1619 | ¹⁵O | 6.859 |
| 2682 | ¹⁸F | 3.724 |
| 3061 | ¹⁸F | 3.839 |
| 3531 | ¹⁸F | 4.652 |
| 3839 | ¹⁸F | 3.839 |
| 4963 | ¹⁸F | 4.963 |
| 5180 | ¹⁵O | 5.180 |
| 6171 | ¹⁵O | 6.176 |
| 6792 | ¹⁵O | 6.793 |

Each notebook follows the same pipeline:

1. **Event selection** — Load the dataset for the given γ energy; select events where the measured AGATA energy lies within tolerance of the Doppler curve computed with `calculate_doppler_energy`.
2. **Kinematics correction** — Use SAURON ejectile (proton) angle and energy to reconstruct the ¹⁸F/¹⁵O recoil via `recoil_kinematics`; apply energy gain and angle offset corrections.
3. **Excitation energy gate** — Reconstruct the recoil excitation energy with `excitation_energy` and apply a narrow window around the level of interest to suppress contaminants.
4. **Δβ calculation** — Compute `Δβ = β_reaction − β_emission` for each selected event.
5. **Sliding-window Gaussian fit** — Bin events in overlapping θ windows; fit a Gaussian with bootstrap to the Δβ distribution in each window to extract its mean and get the covariance matrix due to data overlap.
6. **MCMC Doppler fit** — Fit the θ-dependent Δβ profile with the relativistic Doppler model using `emcee` considering the covariance matrix obtained previously; free parameters are the energy offset ΔE, angle offset Δθ (optional), and the integrated Δβ. Corner and fit plots are saved to `plots/`.

### `Bayesian_Fit.ipynb` — Global lifetime extraction

Combines all per-transition Δβ values with GEANT4 simulation output and literature ¹⁸F lifetimes in a single Bayesian model:

- **Calibration**: fits a linear relation `τ = slope × Δβ + intercept` using the simulation curve as a Gaussian prior on the slope (uncertainty dominated by stopping-power and beam-energy systematics).
- **MCMC** (`emcee`, 32 walkers, 1000 burn-in + 3000 production steps): jointly infers the calibration parameters and individual lifetimes for the newly measured states.
- Posterior samples are saved to `results/lifetime_samples.txt`.

### `Plot_Doppler.ipynb` / `Plot_Simulation.ipynb`

Diagnostic notebooks for visualising the Doppler-shift curves and the simulation-based Δβ–lifetime calibration.

---

## Module Reference

### `doppler.py`

| Function | Description |
|----------|-------------|
| `calculate_doppler_energy(Eg0, beta, theta_deg)` | Relativistic Doppler-shifted γ energy: `Eg = Eg0 √(1−β²) / (1 − β cos θ)` |
| `calculate_doppler_beta(Eg0, energy_offset, angle_offset, fixed_beta)` | Inverse: solve for β given a measured energy curve; returns both solution branches |

### `pykinematic/kinematics.py`

| Function | Description |
|----------|-------------|
| `get_mass(particle)` | Look up atomic mass from the AME table |
| `qvalue(beam, target, ejectile, recoil, ex)` | Q-value for a two-body reaction |
| `calculate_kinematics(beam, target, ejectile, recoil, beam_energy, ex)` | Full kinematic curves (ejectile and recoil) over all lab angles |
| `recoil_kinematics(…, ejectile_energy, angle)` | Reconstruct recoil energy and angle from measured ejectile kinematics |
| `excitation_energy(…, ejectile_energy, angle)` | Reconstruct recoil excitation energy from 4-momentum conservation |
| `energy_to_beta(T, mass)` | Convert kinetic energy to β = v/c |

---

## Data Format

`data/15O_Lifetime_AGATA_{gamma}.npy` (where "gamma" is the gamma energy under analysis) is a NumPy dictionary (`.item()`) keyed by γ energy (integer keV). Each entry contains per-event arrays:

| Key | Description |
|-----|-------------|
| `energy_agata` | AGATA γ energies (list of arrays, one per event) |
| `energy_agata_dc` | Doppler-corrected AGATA energies |
| `theta_agata_sauron` | AGATA detection angles relative to SAURON (rad) |
| `beta_agata` | β at γ emission from Doppler correction |
| `energy_reac` | Ejectile (proton) kinetic energy from SAURON (MeV) |
| `theta_reac_sauron` | Ejectile lab angle from SAURON (rad) |
| `beta_reaction` | Recoil β reconstructed from kinematics |
| `energy_binary` | Reconstructed recoil energy (MeV) |
| `theta_binary` | Reconstructed recoil lab angle (rad) |
| `energy_ex` | Reconstructed excitation energy (MeV) |

---

## Results Format

`results/delta_beta.txt` columns:

```
E_gamma(keV)   delta_beta   delta_beta_err   beta_mean   beta_mean_err
```

`results/simulation.txt` columns (used as calibration prior):

```
tau(fs)   mean_delta_beta   sigma_delta_beta
```

---

## R-Matrix Analysis (`r-matrix/`)

The `r-matrix/` folder contains a Bayesian R-matrix analysis of the astrophysically important **¹⁴N(p,γ)¹⁵O** reaction, performed with the AZURE2 R-matrix code. The measured lifetime of the 6.793 MeV state in ¹⁵O (from this experiment) is fed back in as a prior on the partial radiative width, directly linking the DSAM measurement to the astrophysical S-factor.

### Folder structure

```
r-matrix/
├── 15O.azr                  # AZURE2 input file (level scheme, channel radii, data segments)
├── mcmc.py                  # Bayesian MCMC script (emcee + AZURE2 via pyazr)
├── lifetime_samples.txt     # MCMC posterior samples imported from results/ (6.793 MeV width prior)
│
├── data/                    # Literature cross-section datasets
│   ├── Chen_0.dat / Chen_6793.dat
│   ├── Formicola_F0810_0.dat
│   ├── Frentz_C2889_0.dat / Frentz_C2889_6793.dat
│   ├── Imbriani_O1195_0.dat / Imbriani_O1195_6793.dat
│   ├── Lambert_O0885.dat
│   ├── Li_C3100_0_{1,2}.dat / Li_C3100_6794_{1,2}.dat
│   ├── Marta_O1912_0.dat / Marta_O1912_6793.dat
│   ├── Runkle_C1340_0.dat / Runkle_C1340_6793.dat
│   ├── Schroeder_A0651_0.dat
│   ├── Wagner_O2439_0.dat / Wagner_O2439_6793.dat
│   └── deBoer_C2160_30deg.dat
│
├── output/                  # AZURE2 output files
│   ├── AZUREOut_aa=1_R=1.out        # Cross sections over data energy range
│   ├── AZUREOut_aa=1_R=2.out        # Cross sections over extrapolation range
│   ├── AZUREOut_aa=1_R=2.extrap     # S-factor extrapolation to stellar energies
│   ├── intEC.dat / intEC.extrap     # Internal/external capture contributions
│   ├── chiSquared.out               # Per-segment χ² values
│   ├── normalizations.out           # Dataset normalization factors
│   ├── param.fit / param.par / param.sav   # R-matrix parameter files
│   └── parameters.out               # Full parameter table
│
└── pyazr/                   # Python interface to AZURE2
    ├── azure2.py            # Main wrapper class (spawns AZURE2 servers, manages calculations)
    ├── client.py            # TCP client for communicating with AZURE2 processes
    └── server.py            # Launches AZURE2 as a subprocess on a TCP port
```

### `mcmc.py` — Bayesian R-matrix fit

Runs a joint Bayesian fit of all literature ¹⁴N(p,γ)¹⁵O datasets (ground-state and 6.793 MeV capture) using `emcee`:

- **R-matrix parameters**: reduced-width amplitudes (RWA) and ANCs for all levels, with uninformative flat priors except where constrained.
- **6.793 MeV state width prior**: the lifetime posterior from `lifetime_samples.txt` is converted to a partial width distribution (Γ = ℏ/τ) and used as an empirical KDE prior on the corresponding RWA, propagating the DSAM measurement uncertainty directly into the S-factor.
- **Normalization nuisance parameters**: each dataset is assigned an independent normalization factor with a half-Cauchy prior (scale set per-experiment to reflect published systematic uncertainties).
- **Parallelisation**: `nprocs` AZURE2 instances are spawned on consecutive TCP ports; `emcee` walkers are distributed across them via `multiprocess.Pool`.
- Results are saved to `results/samples-complete.h5` (HDF5 backend, resumable).

### `pyazr/` — AZURE2 Python wrapper

| Class / file | Description |
|---|---|
| `azure2` | Spawns `nprocs` AZURE2 servers, connects TCP clients, exposes `calculate()`, `calculate_sfactor()`, `calculate_chi2()`, parameter I/O, and mode switching (data / extrapolation). |
| `client` | Sends binary-packed command+data buffers over TCP and receives results as NumPy arrays. |
| `server` | Launches a single AZURE2 process listening on a given port. |

### Data files

Each `.dat` file contains ¹⁴N(p,γ)¹⁵O cross-section or S-factor measurements from the literature. Files with suffix `_0` correspond to capture to the ground state; files with suffix `_6793` (or `_6794`) correspond to capture to the 6.793 MeV state. Columns are `E_cm (MeV)`, `angle (deg, 0 for angle-integrated)`, `cross section or S-factor`, `uncertainty`.

### Connection to the main analysis

`lifetime_samples.txt` is generated by `Bayesian_Fit.ipynb` (in the parent directory) and copied here. The script `mcmc.py` reads the first column (lifetime in fs for the 6.793 MeV state), converts to partial width Γ = ℏ/τ, builds a KDE, and uses it as the prior for the corresponding R-matrix parameter — closing the loop between the DSAM lifetime measurement and the astrophysical rate.

---

## Dependencies

```
numpy
scipy
matplotlib
lmfit
emcee
corner
```

Install with:

```bash
pip install numpy scipy matplotlib lmfit emcee corner
```

---

## Usage

Run the notebooks in order:

1. **`Analyze_<E>.ipynb`** for each transition to produce per-transition Δβ values (already stored in `results/delta_beta.txt`).
2. **`Bayesian_Fit.ipynb`** to extract lifetimes from the combined dataset.
3. **`Plot_*.ipynb`** for publication-quality figures.
