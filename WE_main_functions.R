## Symmetric weight function for univariate endpoint ----
#
# Inputs:
#   mu     - true mean response
#   n      - sample size
#   xBar   - sample mean
#   sigma0 - true standard deviation 
#   gamma  - clinical target mean (pre-specified in advance)
#   kappa  - exponent for the sample size (penalization parameter)
#   p      - exponent for the standard deviation
#   log    - logical; if TRUE, returns log-weight
# Output:
#   weight (or log-weight) value at mu

normalWeightFunction = function(mu, n, xBar, sigma0, gamma, kappa, p, log = FALSE) {
  
  # Compute weight standard deviation
  sigmaW = sqrt(sigma0^p / n^kappa)
  
  # Posterior-like combination
  sigmaN = sigma0 / sqrt(n)
  muTilde = (gamma * sigmaN^2 + xBar * sigmaW^2) / (sigmaN^2 + sigmaW^2)
  sigmaTilde = sqrt((sigmaN^2 * sigmaW^2) / (sigmaN^2 + sigmaW^2))
  
  # Standardized ratios
  ratioTilde = muTilde / sigmaTilde
  ratioW = gamma / sigmaW
  ratioN = xBar / sigmaN
  
  # Compute weight (or log-weight)
  if (log) {
    logC = log(sigmaN) - log(sigmaTilde) - 0.5 * (ratioTilde^2 - ratioW^2 - ratioN^2)
    return(logC - 0.5 * (mu - gamma)^2 / sigmaW^2)
  } else {
    C = sigmaN / sigmaTilde * exp(-0.5 * (ratioTilde^2 - ratioW^2 - ratioN^2))
    return(C * exp(-0.5 * (mu - gamma)^2 / sigmaW^2))
  }
  
}

# Vectorizing the previous function
normalWeightFunctionWithMu = Vectorize(normalWeightFunction, vectorize.args = "mu")

## Symmetric weight function for multivariate endpoint ----
#
# Inputs:
#   mu     - true mean vector
#   n      - sample size
#   xBar   - sample mean vector
#   Omega0 - true precision matrix (inverse of covariance matrix)
#   gamma  - clinical target mean vector (pre-specified)
#   OmegaW - weight precision matrix; if NULL, defaults to n^(kappa-1)*OmegaN
#   kappa  - exponent controlling weight scaling
#   log    - logical; if TRUE, returns log-weight
# Output:
#   weight (or log-weight) value at mu

multinormalWeightFunction = function(mu, n, xBar, Omega0, gamma, OmegaW = NULL, kappa = NULL, log = FALSE) {
  
  # Compute precision matrix for the sample
  OmegaN = Omega0 * n
  
  # Set weight precision matrix if not provided
  if (is.null(OmegaW)) {
    OmegaW = n^(kappa - 1) * OmegaN
  }
  
  # Posterior-like combination
  b = OmegaW %*% gamma + OmegaN %*% xBar
  SigmaTilde = solve(OmegaW + OmegaN)
  
  # Quadratic terms for normalization
  addTilde = t(b) %*% SigmaTilde %*% b
  addW = t(gamma) %*% OmegaW %*% gamma
  addN = t(xBar) %*% OmegaN %*% xBar
  
  # Compute weight (or log-weight)
  if (log) {
    logC = -0.5 * log(det(OmegaN %*% SigmaTilde)) - 0.5 * (addTilde - addW - addN)
    return(logC - 0.5 * t(mu - gamma) %*% OmegaW %*% (mu - gamma))
  } else {
    C = (det(OmegaN %*% SigmaTilde))^-0.5 * exp(-0.5 * (addTilde - addW - addN))
    return(C * exp(-0.5 * t(mu - gamma) %*% OmegaW %*% (mu - gamma)))
  }
  
}



## Unweighted Shannon entropy (univariate) ----
#
# Inputs:
#   n      - sample size
#   sigma0 - standard deviation 
# Output:
#   Shannon entropy value
shannonNorm = function(n, sigma0) {
  return(0.5 * (log(2 * pi * exp(1) * sigma0^2) - log(n)))
}

## Weighted Shannon entropy (univariate) ----
#
# Inputs:
#   n      - sample size
#   xBar   - sample mean
#   sigma0 - standard deviation under the null
#   gamma  - clinical target mean
#   kappa  - penalization parameter
#   p      - exponent for standard deviation
# Output:
#   Shannon entropy value
shannonNormWeighted = function(n, xBar, sigma0, gamma, kappa = NULL, p = NULL) {
  
  # Compute weight standard deviation
  sigmaW = sqrt(sigma0^p / n^kappa)
  
  # Posterior-like combination
  sigmaN = sigma0 / sqrt(n)
  muTilde = (gamma * sigmaN^2 + xBar * sigmaW^2) / (sigmaN^2 + sigmaW^2)
  sigmaTilde = sqrt((sigmaN^2 * sigmaW^2) / (sigmaN^2 + sigmaW^2))
  
  # Shannon entropy formula for weighted normal
  H = 0.5 * (log(2 * pi * exp(1) * sigmaN^2)) + (sigmaTilde^2 + (muTilde - xBar)^2) / (2 * sigmaN^2)
  
  return(H)
  
}


## Unweighted Shannon entropy (multivariate) ----
#
# Inputs:
#   n      - sample size
#   Omega0 - precision matrix 
# Output:
#   Shannon entropy value
shannonMultinorm = function(n, Omega0) {
  
  q = ncol(Omega0)
  OmegaN = Omega0 * n
  return(q / 2 * log(2 * pi) - 0.5 * log(det(OmegaN)) + q / 2)
  
}


## Weighted Shannon entropy (multivariate) ----
#
# Inputs:
#   n      - sample size
#   xBar   - sample mean vector
#   Omega0 - precision matrix 
#   gamma  - clinical target mean vector
#   OmegaW - weight precision matrix (optional)
#   kappa  - penalization parameter (optional if OmegaW is NULL)
# Output:
#   Weighted Shannon entropy value
shannonMultinormWeighted = function(n, xBar, Omega0, gamma, OmegaW = NULL, kappa = NULL) {
  
  if (!is.null(kappa) & !is.null(OmegaW)) {
    stop("You cannot specify both OmegaW and kappa at the same time.")
  }
  
  q = ncol(Omega0)
  OmegaN = Omega0 * n
  
  if (is.null(OmegaW)) {
    OmegaW = n^(kappa - 1) * OmegaN
  }
  
  SigmaTilde = solve(OmegaW + OmegaN)
  
  traceMat = sum(diag(OmegaN %*% SigmaTilde))
  weightDistTerm = t(gamma - xBar) %*% OmegaW %*% SigmaTilde %*% OmegaN %*% SigmaTilde %*% OmegaW %*% (gamma - xBar)
  
  return(q / 2 * log(2 * pi) - 0.5 * log(det(OmegaN)) + 0.5 * traceMat + 0.5 * weightDistTerm)
  
}



## Information gain (univariate, symmetric weight) ----
#
# Inputs:
#   n      - sample size
#   xBar   - sample mean
#   sigma0 - standard deviation
#   gamma  - clinical target mean
#   kappa  - penalization parameter
#   p      - SD exponent 
# Output:
#   Information gain value

informationGain = function(n, xBar, sigma0, gamma, kappa, p) {
  
  # Compute penalization factor for the weight
  nPenal = sigma0^(2 - p) * n^kappa / (sigma0^(2 - p) * n^kappa + n)
  
  # Standard deviation under the null
  sigmaN = sigma0 / sqrt(n)
  
  # Information gain formula for symmetric weight
  IG = 0.5 * nPenal - 0.5 * ((gamma - xBar) / sigmaN)^2 * nPenal^2
  
  return(IG)
}

# Vectorized versions 
informationGainWithXBar = Vectorize(informationGain, vectorize.args = "xBar")
informationGainWithSigma0 = Vectorize(informationGain, vectorize.args = "sigma0")
informationGainWithN = Vectorize(informationGain, vectorize.args = "n")
informationGainWithKappa = Vectorize(informationGain, vectorize.args = "kappa")
informationGainWithP = Vectorize(informationGain, vectorize.args = "p")

# Example of multi-argument vectorization
informationGainWithXBarAndN = Vectorize(informationGain, vectorize.args = c("n", "xBar"))
informationGainWithKappaAndN = Vectorize(informationGain, vectorize.args = c("kappa", "n"))
informationGainWithKappaAndNAndXBar = Vectorize(informationGain, vectorize.args = c("kappa", "n", "xBar"))

## Information gain (multivariate, symmetric weight) ----
#
# Inputs:
#   n      - sample size
#   xBar   - sample mean vector
#   Omega0 - precision matrix
#   gamma  - clinical target mean vector
#   kappa  - penalization parameter
# Output:
#   Information gain value

informationGainMultinormal = function(n, xBar, Omega0, gamma, kappa) {
  
  # Number of endpoints
  q = length(xBar)
  
  # Precision matrix for the sample
  OmegaN = Omega0 * n
  
  # Penalization factor for the weight
  nPenal = n^kappa / (n^kappa + n)
  
  # Multivariate information gain formula for symmetric weight
  IG = q / 2 * nPenal - 0.5 * (t(gamma - xBar) %*% OmegaN %*% (gamma - xBar)) * nPenal^2
  
  return(IG)
}


## Posterior probability piHat (Monte Carlo procedure A)  ----
#
# Computes the posterior probability that the estimated best arm is closer
# to the target gamma than the estimated second-best arm for a given trial realization.
#
# Inputs:
#   xBar   - vector of sample means
#   gamma  - vector of clinical target means
#   sigma0 - vector of standard deviations 
#   n      - vector of sample sizes for each arm
#   MCiter - number of Monte Carlo iterations (default: 2e4)
# Output:
#   Posterior probability piHat for a given trial realization

prob.1bestVS2best = function(xBar, gamma, sigma0, n, MCiter = 2e4) {
  
  # Compute absolute distances from gamma
  distFromGamma = abs(gamma - xBar)
  
  # Find the indices of the best and second-best arms
  idx1best = which(distFromGamma == sort(distFromGamma, decreasing = FALSE)[1])
  idx2best = which(distFromGamma == sort(distFromGamma, decreasing = FALSE)[2])
  
  # Monte Carlo draws for posterior distances
  postDistance1best = abs(
    rnorm(MCiter, mean = xBar[idx1best], sd = sigma0[idx1best] / sqrt(n[idx1best])) - gamma[idx1best]
  )
  
  postDistance2best = abs(
    rnorm(MCiter, mean = xBar[idx2best], sd = sigma0[idx2best] / sqrt(n[idx2best])) - gamma[idx2best]
  )
  
  # Posterior probability that best arm is closer to gamma than second-best
  prob = mean(postDistance1best < postDistance2best)
  
  return(prob)
}
