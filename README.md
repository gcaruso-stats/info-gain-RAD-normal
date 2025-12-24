# R code from *A response-adaptive multi-arm design for continuous endpoints based on a weighted information measure* by G. Caruso and P. Mozgunov (2024+)

This repository contains R functions to implement the procedures presented in the original paper *A response-adaptive multi-arm design for continuous endpoints based on a weighted information measure* by G. Caruso and P. Mozgunov (2024+).


## Repository Structure

- `info-gain-RAD-normal/`
  - `WE_main_functions_univ.R` – All core functions to run procedures described in the paper (univariate case)
  - `WE_main_functions_multiv.R` – All core functions to run procedures described in the paper (multivariate case)
  - `run_trials_univ.R` – Code to run single-trial and multi-trials simulations under all the considered designs (univariate endpoint).
  - `run_trials_multiv.R` – Code to run single-trial and multi-trials simulations under all the considered designs (multivariate endpoint).

## Main Functions (both univariate and multivariate)

1. **Weight function**  
   Computes the symmetric weight function with Gaussian kernel introduced in the paper.

2. **Shannon entropies**  
   Computes the unweighted and weighted Shannon entropy of a normal distribution. 

3. **Information gain**  
   Computes the information gain as the difference between the unweighted Shannon entropy and the weighted Shannon entropy.

4. **Single trial simulation**  
   Runs a single trial under the specified design.

5. **Bayesian hypothesis testing**  
   Performs Bayesian decision-making for a single trial realization (`Monte Carlo procedure A` - details in the paper).

6. **Frequentist probability of rejection**  
   Computes the probability of rejecting the null hypothesis (`Monte Carlo procedure B` - details in the paper).

7. **Derive critical values**  
   Calculates critical values to achieve strong or mean type-I error control across a grid of null scenarios.

8. **Robust kappa selection**  
   Selects a robust value of kappa over a grid of candidate values.




