# R code from *A response-adaptive multi-arm design for continuous endpoints based on a weighted information measure* by G. Caruso and P. Mozgunov (2024+)

This repository contains R functions to implement the procedures presented in the original paper *A response-adaptive multi-arm design for continuous endpoints based on a weighted information measure* by G. Caruso and P. Mozgunov (2024+).

**Status:** this repository is under active development. Core components of the methodology presented in the paper are available, while remaining parts are being polished, commented and added progressively.


## Repository Structure

- `info-gain-RAD-normal/`
  - `WE_main_functions_univ.R` – All core functions to run procedures described in the paper (univariate case).
  - `run_trials_univ.R` – Code to run single-trial and multi-trials simulations under all the considered designs (univariate endpoint).
  - `WE_main_functions_multiv.R` – All core functions to run procedures described in the paper (multivariate case)
  - `run_trials_multiv.R` – Code to run single-trial and multi-trials simulations under all the considered designs (multivariate endpoint).
  - `gittinsTableFull.csv` - Values for the normalised Gittins Index from the table 8.1 in Gittins et al. (2011, p.261-262); missing rows in the table have been estimated via linear interpolation.

## Main Functions (both univariate and multivariate)

1. **Weight function**: computes the symmetric weight function with Gaussian kernel introduced in the paper.

2. **Shannon entropies**: computes the unweighted and weighted Shannon entropy of a normal distribution. 

3. **Information gain**: computes the information gain as the difference between the unweighted Shannon entropy and the weighted Shannon entropy.

4. **Single trial simulation**: runs a single trial under the specified design.

5. **Summaries of the single trial**: computes several summaries for a single trial.

7. **Bayesian hypothesis testing**: performs Bayesian decision-making for a single trial realization (`Monte Carlo procedure A` - details in the paper).

8. **Multiple trial simulation**: runs multiple trial replicates under the specified design.

9. **Operating characteristics over many replicates**: compute operating characteristics over many trial replicates. The probability of rejecting the null hypothesis is computed following `Monte Carlo procedure B` (details in the paper).

10. **Derive critical values**: calculates critical values to achieve strong or mean type-I error control across a grid of null scenarios.

11. **Robust kappa selection**: selects a robust value of kappa over a grid of candidate values.




