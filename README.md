# Bayesian Instrumental Variable Analysis

Research code for Bayesian instrumental variable (IV) analysis with censored
time-to-event outcomes. The repository contains implementations and analysis
scripts for parametric Bayesian IV models and semiparametric Dirichlet process
mixture IV (DPMIV) models.

This repository supports work on IV analysis for partly interval-censored
time-to-event outcomes, including simulation studies and applications to
biomedical cohort data.

## Related Work

- Parametric Bayesian IV method for survival outcomes:
  [Bayesian instrumental variable analysis with censored time-to-event outcomes](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4314427/).
- Semiparametric DPMIV extension for partly interval-censored outcomes:
  *A Semiparametric Bayesian Method for Instrumental Variable Analysis with
  Partly Interval-Censored Time-to-Event Outcomes*.
- Related dissertation material:
  [Bayesian Instrumental Variable Analysis for Survival Outcomes](https://escholarship.org/uc/item/8223z6fp).

## Repository Structure

| Path | Contents |
| --- | --- |
| `PBIV/` | R code for parametric Bayesian IV models, including right-censored and interval-censored variants. |
| `DPMIV_Simulation/` | C source code for simulation studies under several bivariate error distributions. |
| `DPMIV_UKB/` | Code and study-specific analysis scripts for the UK Biobank example. |
| `DPMIV_ARIC/` | Code and study-specific analysis scripts for the ARIC example. |
| `DPMIV_Analysis/` | Downstream scripts for posterior summaries, comparisons, diagnostics, and figures. |

The scripts preserve the analysis workflow used for the manuscript. They are
research code rather than a packaged R library, so paths, input files, and
runtime settings may need to be adapted for a new computing environment.

## Data and Permissions

The UK Biobank and Atherosclerosis Risk in Communities (ARIC) Study data used
in the manuscript are controlled-access cohort data. Researchers must obtain
access through the respective study data-access procedures and comply with the
applicable data-use agreements.

The authors are not permitted to redistribute individual-level UK Biobank or
ARIC data as open data. Before reusing any study-specific input files, users are
responsible for confirming that they have the required permissions. Simulation
code can be run independently of these cohort data.

## Requirements

The code uses a mixture of R scripts and C programs.

Common R package dependencies include:

- `reshape2`
- `tidyverse`
- `bayesplot`
- `icenReg`
- `bayesSurv`
- `coda`
- `viridis`
- `jcolors`
- `MASS`

The DPMIV simulation and data-analysis programs require a C compiler such as
`gcc` or `clang`. Exact run times depend on the number of MCMC iterations,
number of chains, and thinning settings.

## Basic Use

1. Clone the repository.
2. Install the R package dependencies needed for the relevant script.
3. Compile the C program in the relevant DPMIV directory when using the
   semiparametric sampler.
4. Run the simulation or analysis scripts from their corresponding directory,
   checking file paths and MCMC settings before launch.

For manuscript reproduction, start with the directory matching the target
analysis:

- `PBIV/` for parametric Bayesian IV comparisons.
- `DPMIV_Simulation/` for simulation scenarios.
- `DPMIV_UKB/` and `DPMIV_ARIC/` for the real-data examples.
- `DPMIV_Analysis/` for posterior summaries and figures.

## Authors

- [Gang Li](https://gangli.faculty.biostat.ucla.edu/), Department of
  Biostatistics, UCLA.
- Xuyang Lu, Department of Biostatistics, UCLA.
- [Eliuvish Han Cui](https://elviscuihan.github.io/), Department of
  Biostatistics, UCLA.

## License

This repository is released under the MIT License. See [`LICENSE`](LICENSE) for
details. Data-use restrictions from UK Biobank, ARIC, and other third-party data
providers are separate from the software license and must be followed
independently.
