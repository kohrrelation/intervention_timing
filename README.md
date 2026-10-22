# Code and data for
# *Reactive management timing shapes bias from omitted temporal risk in ecological intervention studies*

This archive contains the R code and data required to reproduce the simulation study and empirical analysis presented in the manuscript:

> **Reactive management timing shapes bias from omitted temporal risk in ecological intervention studies**

Author: Jonathan Koh

## Contents

### R code

- `0_functions.R`
  - Functions for simulating the managed self-exciting point process.
  - Functions for constructing ecological-event and intervention-history processes.
  - Point process likelihoods and model-fitting functions.
  - Functions used in the simulation experiments.

- `1_simulation_study.R`
  - Reproduces the simulation experiments in the manuscript.
  - Generates Figures 2–5.
  - Examines:
    1. early versus quiet reactive-management protocols;
    2. intervention timing relative to the persistence timescale of event risk;
    3. the relationship between empirical bias and the theoretical covariance-over-variance expression;
    4. variation in bias with self-excitation strength.

- `2_data_analysis.R`
  - Reproduces the empirical analysis of livestock depredations and wolf removals in France.
  - Generates Figure 6.
  - Includes the matched pseudo-removal analysis.

### Data

The `Data/` directory contains the data used in the empirical illustration:

- `wolves_processed.RData`

These data originate from the wolf-management study described in:

Grente, O., Opitz, T., Duchamp, C., Drouet-Hoguet, N., Chamaillé-Jammes, S. & Gimenez, O. (2025).  
*Spatio-temporal dynamics of attacks around deaths of wolves: A statistical assessment of lethal control efficiency in France*.  
*Peer Community Journal*, 5, e128.

The empirical data retain the terms and conditions of the original data source. The MIT licence included in this archive applies to the R code only and does not alter the licence, ownership, or reuse conditions of the empirical data.

## Software requirements

The analyses were conducted using:

- **R version 4.6.1 (2026-06-24)**

The empirical analysis uses the following R packages:

- `dplyr`
- `tidyr`
- `ggplot2`
- `lubridate`

Package versions for the final analysis can be recorded using:

```r
sessionInfo()
```

## Directory structure

Run the scripts from the root directory of the archive. The expected structure is:

```text
.
├── 0_functions.R
├── 1_simulation_study.R
├── 2_data_analysis.R
├── LICENSE.txt
├── README.md
├── Data/
│   ├── wolves_processed.RData
└── Plots/
```

If the `Plots/` directory is not already present, create it in R using:

```r
dir.create("Plots", showWarnings = FALSE)
```

## Reproducing the simulation study

Run:

```r
source("1_simulation_study.R")
```

The script sources `0_functions.R` and reproduces the simulation analyses and Figures 2–5.

The main baseline parameter values are:

- observation period: `T = 500`
- background event rate: `mu = 0.2`
- self-excitation decay rate: `beta_N = 1`
- intervention-effect decay rate: `beta_R = 1`
- baseline branching ratio: `eta = 0.6`
- management trigger: at least `K = 3` events within a window of length `W = 4`
- baseline early-management delay: `D = 0.1`
- quiet-management period: `Q = 3`

Following an intervention, the triggering mechanism is reset and no new intervention can be triggered for three time units.

The fitted point process likelihood evaluates the compensator numerically using a midpoint grid of 4,000 time points. The decay rates `beta_N` and `beta_R` are fixed at their data-generating values in the fitted models.

Fixed random-number seeds are used in the simulation analyses for reproducibility.

## Reproducing the empirical wolf-management analysis

Run:

```r
source("2_data_analysis.R")
```

The script reads the attack and removal data from the `Data/` directory and generates the three panels of Figure 6.

The empirical analysis uses:

- a spatial radius of 10 km;
- a 60-day lookback period;
- 500 Monte Carlo pseudo-removal replicates.

### Matched pseudo-removal analysis

For each observed wolf removal, matched pseudo-removals:

1. retain the observed removal location;
2. retain the observed calendar month and day;
3. retain the same broad management period (2011–2014, 2015–2018, or 2019–2020);
4. assign a year sampled uniformly from the other eligible years within that management period;
5. use only years for which the complete preceding 60-day observation window is available.

This preserves location, seasonality, and broad temporal changes in attack occurrence across the study period, while breaking the association between the observed removal and local attack activity in its actual removal year.

For Figure 6b, local attack histories are summarized in successive five-day intervals over the 60 days preceding each observed or pseudo-removal. The observed mean curve is compared with the Monte Carlo mean and pointwise 95% simulation envelope.

For Figure 6c, for each threshold `q = 1, ..., 30` days, the analysis records whether at least one attack occurred within 10 km during the preceding `q` days. The cumulative proportion therefore represents the proportion of all removals occurring within `q` days of a local attack. Removals with no attack during the preceding 60 days remain in the denominator and contribute zero at all thresholds shown.

The manuscript additionally reports the proportions occurring within three and seven days, with Monte Carlo p-values calculated using the plus-one correction:

```text
p = (1 + number of simulated statistics at least as extreme as observed) /
    (number of Monte Carlo replicates + 1)
```

## Output

The scripts write figure files to the `Plots/` directory.

The principal outputs are:

- Figure 2: bias under early and quiet management;
- Figure 3: bias as a function of intervention timing;
- Figure 4: empirical bias versus the covariance-over-variance quantity;
- Figure 5: bias across intervention timing and self-excitation strength;
- Figure 6: empirical wolf attack and removal analyses.

## Reproducibility notes

The simulation study may take some time to reproduce because several hundred independent point process realizations are generated and fitted for each scenario.

Small numerical differences may occur across operating systems or package versions because the analyses involve numerical optimization and Monte Carlo simulation. Fixed random-number seeds are used where appropriate.

## Licence

The R code in this archive is released under the MIT License; see `LICENSE.txt`.

The empirical data are not covered by the MIT License and remain subject to the terms of the original data source.

## Contact

Jonathan Koh, jonathan.koh@stat.math.ethz.ch
