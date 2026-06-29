## ===================================================================================
## Main functions from "A response-adaptive multi-arm design for continuous endpoints 
## based on a weighted information measure" by G. Caruso and P. Mozgunov (2024+)
## ===================================================================================

## Symmetric weight function ----
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

## Unweighted Shannon entropy  ----
#
# Inputs:
#   n      - sample size
#   sigma0 - standard deviation 
# Output:
#   Shannon entropy value
shannonNorm = function(n, sigma0) {
  return(0.5 * (log(2 * pi * exp(1) * sigma0^2) - log(n)))
}

## Weighted Shannon entropy  ----
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


## Information gain (symmetric weight) ----
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

## Run a single trial ----
#
# This function runs a single simulated trial under a specified allocation
# design (details below). It performs patients' treatment assignments, generates
# responses and stores all quantities needed for post-trial summaries.
#
# INPUTS:
#
# - N: total sample size
# - K: number of treatment arms
# - mu: vector of true mean responses (length K)
# - sigma0: vector of true standard deviations (length K)
# - gamma: vector of target values (length K)
#
# - areSigma0Known: logical; if TRUE, sigma0 is treated as known
# - initialSizePerArm: number of burn-in observations per arm
#
# - design: character string specifying the allocation rule.
#   Allowed values:
#     * "FR"  : Fixed and equal randomization
#     * "RAD" : Response-adaptive design based on information gain derived from symmetric weight function (WE in the paper)
#     * "CB"  : Current Belief design
#     * "GI"  : Gittins Index-based designs (symmetric and targeted)
#     * "TS"  : Thompson Sampling-type design (adjusted version)
#
# - parDesign: list containing design-specific parameters.
#   Its required fields depend on the chosen design:
#
#   * FR:
#       - parDesign = NULL
#
#   * RAD:
#       - kappa      : penalization parameter
#       - p          : exponent for the standard deviation
#
#   * CB:
#       - parDesign = NULL
#
#   * GI:
#       - GItype : "symm" or "targ" 
#       - d      : discount factor for the Gittins index
#
#   * TS:
#       - MCiter : number of Monte Carlo samples used to approximate posterior
#
# OUTPUT:
#
# A list with the following components:
#
# - trialDetails: list with N, K, mu, sigma0, gamma, true best and second best arms
# - design, parDesign, areSigma0Known, initialSizePerArm
# - a: vector of treatment assignments (length N)
# - x: N x K matrix of observed responses (NA when patient i is not assigned to arm j)
# - n: vector with the number of patients assigned to each arm (length K)


runTrial = function(N, K, mu, sigma0, gamma,
                    areSigma0Known = TRUE,
                    initialSizePerArm = 2,
                    design,
                    parDesign = NULL) {
  
  if(!any(design == c("RAD","FR","CB","GI","TS"))){
    stop("Invalid design name.")
  }
  
  trueBestArm = which.min(abs(gamma - mu))
  trueSecondBestArm = which(abs(gamma - mu) == sort(abs(gamma - mu))[2])[1]
  
  outList = list(
    trialDetails = list(
      N = N,
      K = K,
      mu = mu,
      sigma0 = sigma0,
      gamma = gamma,
      trueBestArm = trueBestArm,
      trueSecondBestArm = trueSecondBestArm
    ),
    areSigma0Known = areSigma0Known,
    initialSizePerArm = initialSizePerArm,
    design = design,
    parDesign = parDesign
  )
  
  
  ## ------------------------
  ## START ALLOCATION
  ## ------------------------
  
  # Burn-in phase (if the design is response-adaptive)
  
  if(design != "FR"){
    
    a = rep(NA, N)
    x = matrix(NA, nrow = N, ncol = K)
    
    a[1:(K * initialSizePerArm)] = rep(1:K, initialSizePerArm)
    
    for(t in 1:(K * initialSizePerArm)){
      x[t, a[t]] = rnorm(1, mu[a[t]], sigma0[a[t]])
    }
    
    aveResponseByArm = colMeans(x, na.rm = TRUE)
    sdResponseByArm = apply(x, 2, sd, na.rm = TRUE)
  }
  
  ## ---------- FR ----------
  
  if(design == "FR"){
    
    n = NULL
    while(length(n) < K){
      a = sample.int(K, size = N, replace = TRUE)
      n = as.numeric(table(a))
    }
    
    x = matrix(NA, nrow = N, ncol = K)
    for(t in 1:N){
      x[t, a[t]] = rnorm(1, mu[a[t]], sigma0[a[t]])
    }
    
    aveResponseByArm = colMeans(x, na.rm = TRUE)
    sdResponseByArm = apply(x, 2, sd, na.rm = TRUE)
  }
  
  ## ---------- RAD (WE design) ----------
  
  if(design == "RAD"){
    
    delta = matrix(NA, nrow = N, ncol = K)
    
    for(t in (K * initialSizePerArm + 1):N){
      for(j in 1:K){
        delta[t, j] = informationGain(
          n = length(which(a==j)),
          xBar = aveResponseByArm[j],
          sigma0 = ifelse(areSigma0Known, sigma0[j], sdResponseByArm[j]),
          gamma = gamma[j],
          kappa = parDesign$kappa,
          p = parDesign$p
        )
      }
      
      a[t] = which.max(delta[t, ])
      x[t, a[t]] = rnorm(1, mu[a[t]], sigma0[a[t]])
      
      aveResponseByArm = colMeans(x, na.rm = TRUE)
      sdResponseByArm = apply(x, 2, sd, na.rm = TRUE)
    }
    
    n = as.numeric(table(a))
    outList$parDesign$out$delta = delta
  }
  
  ## ---------- CB ----------
  
  if(design=="CB"){
    
    for(t in (K*initialSizePerArm+1):N){
      
      #assign the patient to the arm with closest mean response to the target
      a[t] = which.min(abs(aveResponseByArm-gamma))
      
      #observe the response for that patient
      x[t,a[t]] = rnorm(1,mu[a[t]],sigma0[a[t]])
      
      aveResponseByArm = colMeans(x, na.rm = T)
      sdResponseByArm = apply(x, 2, sd, na.rm=T)
      
    }
    
    n = as.numeric(table(a)) 
    
  }
  
  ## ---------- GI ----------
  
  if(design=="GI"){
    
    GItabfull = read.csv("gittinsTableFull.csv",header = T) #load a table with normalised Gittins Index (d=0.5,0.6,0.7,0.8,0.9,0.95,0.99,0.995 and n=1,...,500)
    GInormfull = GItabfull[,paste0("d",gsub("\\.","",parDesign$d))] #select the required column (according to d)
    
    GI = matrix(rep(NA, N*K), ncol=K) #N x K matrix with values of the Gittins Index at each step for each arm
    
    for(t in (K*initialSizePerArm+1):N){
      
      s=ifelse(rep(areSigma0Known==T,K),
               sigma0,
               sdResponseByArm)
      
      if(parDesign$GItype=="symm"){
        GI[t,] = GInormfull[table(a)]*s - abs(aveResponseByArm-gamma) 
        
        #assign the patient to the arm with largest symmetric Gittins index
        a[t] = which.max(GI[t,])
        
      }
      if(parDesign$GItype=="targ"){
        GI[t,] = aveResponseByArm + ifelse(aveResponseByArm>gamma,-1,1)*GInormfull[table(a)]*s
        
        #assign the patient to the arm with Gittins index closest to the target
        a[t] = which.min(abs(GI[t,]-gamma))
        
      }
      
      #observe the response for that patient
      x[t,a[t]] = rnorm(1,mu[a[t]],sigma0[a[t]])
      
      aveResponseByArm = colMeans(x, na.rm = T) 
      sdResponseByArm = apply(x, 2, sd, na.rm=T) 
      
    }
    
    n = as.numeric(table(a)) 
    
    outList$parDesign$out$GInormfull = GInormfull
    outList$parDesign$out$GI = GI
    
  }
  
  ## ---------- TS ----------
  
  if(design=="TS"){
    
    require(Rfast)
    
    zSim = matrix(rnorm(K*parDesign$MCiter),ncol=K) 
    postDistance = sapply(1:K, function(j) abs( aveResponseByArm[j]-gamma[j] + ifelse(areSigma0Known==T,sigma0[j],sdResponseByArm[j])/sqrt(initialSizePerArm) * zSim[,j] )) #samples from posterior of |mu_j-gamma|, for all j 
    probBestByArm = matrix(rep(NA, N*K), ncol=K) #N x K matrix with probabilities that a particular arm is the best one for each arm (columns) and at for each iteration (rows)
    probBestByArm_adj = matrix(rep(NA, N*K), ncol=K)
    
    for(t in (K*initialSizePerArm+1):N){
      
      probBestByArm[t,] = table(factor(rowMins(postDistance), levels = 1:K))/nrow(postDistance) #N.B. TIME-CONSUMING LINE
      
      probBestByArm_adj[t,] = probBestByArm[t,]^(t/(2*N))/sum(probBestByArm[t,]^(t/(2*N))) 
      
      a[t] = which.max(probBestByArm_adj[t,]) #assign patient to the arm with highest posterior "adjusted" probability of being the best
      
      
      #observe the response for that patient
      x[t,a[t]] = rnorm(1,mu[a[t]],sigma0[a[t]])
      
      aveResponseByArm = colMeans(x, na.rm = T)
      sdResponseByArm = apply(x, 2, sd, na.rm=T)
      
      postDistance[,a[t]] = abs( aveResponseByArm[a[t]]-gamma[a[t]] + ifelse(areSigma0Known==T,sigma0[a[t]],sdResponseByArm[a[t]])/sqrt(table(a)[a[t]]) * zSim[,a[t]] ) #update the posterior sample of |mu_j-gamma|, where j is the last played arm 
      
    }
    
    n = as.numeric(table(a)) 
    
    outList$parDesign$out$probEstByArm = probBestByArm
    outList$parDesign$out$probEstByArm_adj = probBestByArm_adj
    
  }
  
  
  ## ------------------------
  ## END ALLOCATION
  ## ------------------------
  
  outList$a = a
  outList$x = x
  outList$n = as.numeric(table(a))
  
  return(outList)
}


## Compute summaries for a single trial ----
#
# This function computes summaries at the end of a single trial
#
# INPUT:
#   trialOut : output of runTrial() (list)
#
# OUTPUT:
#   A list with the following elements:
#
#   - N, K
#   - selectedBestArm           : arm whose sample mean is closest to gamma
#   - selectedSecondBestArm     : estimated second closest arm to gamma
#   - trueBestArm               : arm whose true mean is closest to gamma
#   - trueSecondBestArm         : true second closest arm to gamma
#   - correctSelectionBest      : is best arm correctly identified? yes (1), no (0)
#   - correctSelectionSecondBest: is second best arm correctly identified? yes (1), no (0)
#   - correctSelectionTop2      : are two best arms correctly identified? yes (1), no (0)
#   - nByArm                    : number of patients allocated to each arm
#   - allocationProportions     : proportion of patients allocated to each arm
#   - aveResponseByArm          : sample mean by arm
#   - sdResponseByArm           : sample standard deviation by arm
#   - trueDistanceToTarget      : |mu_j - gamma_j| for each arm
#   - estimatedDistanceToTarget : |xBar_j - gamma_j| for each arm
#   - postProb_1bestVS2best     : posterior probability that |mu1best - gamma| < |mu2best - gamma|

computeTrialSummaries = function(trialOut) {
  
  a = trialOut$a
  x = trialOut$x
  n = trialOut$n
  
  mu = trialOut$trialDetails$mu
  gamma = trialOut$trialDetails$gamma
  sigma0 = trialOut$trialDetails$sigma0
  
  areSigma0Known =  trialOut$areSigma0Known
  
  trueBestArm = trialOut$trialDetails$trueBestArm
  trueSecondBestArm = trialOut$trialDetails$trueSecondBestArm
  
  K = length(mu)
  N = sum(n)
  
  # Sample means and sd's
  
  aveResponseByArm = colMeans(x, na.rm = TRUE)
  sdResponseByArm = apply(x, 2, sd, na.rm = TRUE)
  
  # final arm recommendations
  
  selectedBestArm = which.min(abs(aveResponseByArm - gamma))
  selectedSecondBestArm = which(abs(aveResponseByArm - gamma) == sort(abs(aveResponseByArm - gamma))[2])[1]
  
  # Correct selection
  
  correctSelectionBest = as.integer(selectedBestArm == trueBestArm)
  correctSelectionSecondBest = as.integer(selectedSecondBestArm == trueSecondBestArm)
  correctSelectionTop2 = correctSelectionBest * correctSelectionSecondBest #1 iff top2 correctly identified; 0, otherwise
  
  # Distances from target
  
  trueDistanceToTarget = abs(mu - gamma)
  estimatedDistanceToTarget = abs(aveResponseByArm - gamma)
  
  # Bayesian test
  
  postProb_1bestVS2best = prob.1bestVS2best(xBar = aveResponseByArm,
                                            gamma = gamma,
                                            sigma0 = ifelse(rep(areSigma0Known, K), 
                                                            sigma0, 
                                                            sdResponseByArm),
                                            n = n)
  
  # Output
  
  outSummary = list(
    
    N = N,
    K = K,
    
    selectedBestArm = selectedBestArm,
    selectedSecondBestArm = selectedSecondBestArm,
    
    trueBestArm = trueBestArm,
    trueSecondBestArm = trueSecondBestArm,
    
    correctSelectionBest = correctSelectionBest,
    correctSelectionSecondBest = correctSelectionSecondBest,
    correctSelectionTop2 = correctSelectionTop2,
    
    nByArm = n,
    allocationProportions = n / N,
    
    aveResponseByArm = aveResponseByArm,
    sdResponseByArm = sdResponseByArm,
    
    trueDistanceToTarget = trueDistanceToTarget,
    estimatedDistanceToTarget = estimatedDistanceToTarget,
    
    postProb_1bestVS2best = postProb_1bestVS2best
  )
  
  return(outSummary)
}

## Run multiple trials ----
#
# This function runs multiple (M) trial replicates. 
#
# Inputs:
#   M                         - number of simulated trial replicates
#   N, K, mu, sigma0, gamma   - trial parameters
#   design, parDesign         - design specification
#   seed                      - RNG seed for reproducibility
#
# Output:
#   A list with:
#     - trials: list of M trial objects (output of runTrial)
#     - elapsedTime: total simulation time

runMultipleTrials = function(M = 5e2,
                             N,
                             K,
                             mu,
                             sigma0,
                             gamma,
                             areSigma0Known = TRUE,
                             initialSizePerArm = 2,
                             design,
                             parDesign = NULL,
                             seed = 123,
                             printTimeElaps = TRUE) {
  
  set.seed(seed)
  startTime = Sys.time()
  
  trials = replicate(
    M,
    {
    trial = runTrial(N = N,
             K = K,
             mu = mu,
             sigma0 = sigma0,
             gamma = gamma,
             areSigma0Known = areSigma0Known,
             initialSizePerArm = initialSizePerArm,
             design = design,
             parDesign = parDesign)
     list(trial = trial, summary = computeTrialSummaries(trial))
     },
    simplify = FALSE)
  
  
  elapsedTime = Sys.time() - startTime
  
  if (printTimeElaps) {
    message("Time elapsed for simulations: ",round(elapsedTime, 2), " ", attr(elapsedTime, "units"))
  }
  
  return(list(trials = trials, elapsedTime = elapsedTime))
}

## Operating characteristics over many replicates (Monte Carlo procedure B embedded here) ----
#
# This function computes operating characteristics over many trial replicates. 
# It embeds Monte Carlo procedure B to compute frequentist rejection rates.
#
# Inputs:
#   trials      - list of trial objects (output of runTrial)
#   cutoffProb  - vector of Bayesian posterior probability cutoffs
#
# Output:
#   List of operating characteristics:
#     - allocation summaries
#     - selection probabilities
#     - response and target-distance metrics
#     - Bayesian power (naive, two-comp, conditional)

computeOCs = function(trials, cutoffProb = NULL) {
  
  M = length(trials)
  K = trials[[1]]$summary$K
  N = trials[[1]]$summary$N
  
  # Allocation summaries
  
  nPercentMat = sapply(trials, function(x) x$summary$nByArm)/N
  
  armPerc_average = rowMeans(nPercentMat)
  armPerc_SE = apply(nPercentMat, 1, sd) / sqrt(M)
  
  # Final recommendation
  
  armFinal = sapply(trials, function(x) x$summary$selectedBestArm)
  armFinalPerc_average = table(factor(armFinal, levels = 1:K)) / M * 100
  armFinalModal = as.numeric(which.max(armFinalPerc_average))
  
  # Average responses by arm
  
  aveResponse = rowMeans(sapply(trials, function(x) x$summary$aveResponseByArm))
  response_SE = apply(sapply(trials, function(x) x$summary$aveResponseByArm), 1, sd) / sqrt(M) # Monte Carlo Standard Error
  
  # Estimation performance (taking arm-specific average response as estimator)
  
  mu = trials[[1]]$trial$trialDetails$mu # underlying truth
  aveResponse_bias = aveResponse - mu
  aveResponse_var = apply(sapply(trials, function(x) x$summary$aveResponseByArm), 1, var)*(M-1)/M
  aveResponse_mse = aveResponse_bias^2 + aveResponse_var
  # mse = rowMeans((sapply(trials, function(x) x$summary$aveResponseByArm) - mu)^2) #alternative way
  
  # Average distance from the target
  
  aveTargetDistance = rowMeans(sapply(trials, function(x) x$summary$estimatedDistanceToTarget))
  targetDistance_SE = apply(sapply(trials, function(x) x$summary$estimatedDistanceToTarget), 1, sd) / sqrt(M)
  
  # Correct identification best arm(s)
  
  selectionBest = sapply(trials, function(x) x$summary$correctSelectionBest) #1 if best arm correctly identified; 0, otherwise
  selectionSecondBest = sapply(trials, function(x) x$summary$correctSelectionSecondBest) #1 if second best arm correctly identified; 0, otherwise
  selectionTop2 = sapply(trials, function(x) x$summary$correctSelectionTop2) #1 if two best arms correctly identified; 0, otherwise
  
  percSelectionBest = mean(selectionBest)
  percSelectionSecondBest = mean(selectionSecondBest)
  percSelectionTop2 = mean(selectionTop2)
  
  # Frequentist rejection rates (power or type-I error)
  
  TS_bayes = sapply(trials, function(x) x$summary$postProb_1bestVS2best)
  TScond_bayes = ifelse(selectionTop2, TS_bayes, NA)
  
  powerNaive = sapply(cutoffProb, function(c) mean(TS_bayes > c)) # NOT CONSIDERED IN THE PAPER
  powerTwoComp = sapply(cutoffProb, function(c) mean(selectionTop2 & TS_bayes > c))
  powerCond = sapply(cutoffProb, function(c) mean(TScond_bayes > c, na.rm = TRUE))
  
  return(list(
    armPerc_average = armPerc_average,
    armPerc_SE = armPerc_SE,
    armFinalPerc_average = armFinalPerc_average,
    armFinalModal = armFinalModal,
    aveResponse = aveResponse,
    response_SE = response_SE,
    aveResponse_bias = aveResponse_bias,
    aveResponse_var = aveResponse_var,
    aveResponse_mse = aveResponse_mse,
    aveTargetDistance = aveTargetDistance,
    targetDistance_SE = targetDistance_SE,
    percSelectionBest = percSelectionBest,
    percSelectionSecondBest = percSelectionSecondBest,
    percSelectionTop2 = percSelectionTop2,
    cutoffProb = cutoffProb,
    powerNaive = powerNaive,
    powerTwoComp = powerTwoComp,
    powerCond = powerCond
  ))
}


## Cut-off probability search for Bayesian hypothesis testing ----
#
# This function calculates critical values to achieve strong or 
# mean type-I error control across a grid of null scenarios when  
# Bayesian hypothesis testing (prob.1bestVS2best) is performed.
#
# Inputs:
#   M                 - number of simulated trials per null scenario
#   K                 - number of arms
#   N                 - total sample size
#   nullScenarios     - list of null scenarios; each element is a list 
#                       containing a vector of true means 'mu' (length K) and 
#                       a vector of true standard deviations 'sigma' (length K)
#   gamma             - vector of target values (length K)
#   areSigma0Known    - logical; whether variances are known or not
#   initialSizePerArm - burn-in sample size per arm
#   design            - trial design ("FR","CB","RAD","GI","TS") - more details 
#                       can be found in 'Run a single trial' section
#   parDesign         - design-specific parameters - more details 
#                       can be found in 'Run a single trial' section
#   alphaTarg         - maximum or average (across the set of null scenarios) 
#                       tolerated type-I-error rate (TIE)
#   alphaMaxSingle    - maximum tolerated TIE on a single scenario; 
#                       this cap may avoid extreme type-I error rates on single
#                       scenarios, when cut-off values are derived based on average TIE
#                       (default is 1, i.e. no cap)
#   ncores            - number of cores for parallel computing
#   seed              - RNG seed for reproducibility
#
# Output:
#   A list containing:
#     - draws from posterior probabilities (prob.1bestVS2best) under the null scenarios 
#     - cut-off probabilities ensuring strong and mean type-I error control
#     - scenario-specific type-I error rates
#     - elapsed computing time

criticalValueSearchBayes = function(M, K, N,
                                    nullScenarios, gamma,
                                    areSigma0Known = FALSE, initialSizePerArm = 2,
                                    design, parDesign = NULL, 
                                    alphaTarg = 0.05, alphaMaxSingle = 1,
                                    ncores = NULL,
                                    seed = 123) {
  
  set.seed(seed)
  startTime = Sys.time()
  
  # Parallel setup
  
  require(doParallel)
  require(doRNG)
  
  if (is.null(ncores)) {
    ncores = detectCores() - 1
  }
  
  cl = makeCluster(ncores, outfile = "")
  registerDoParallel(cl)
  
  # Simulate null scenarios
  
  bayesNull = foreach(i = seq_along(nullScenarios),
                      .export = c("computeTrialSummaries", "runTrial", "prob.1bestVS2best", "informationGain"),
                      .inorder = TRUE) %dorng% {
    
    trials = replicate(M,
                       runTrial(N = N, K = K,
                                mu = nullScenarios[[i]]$mu, 
                                sigma0 = nullScenarios[[i]]$sigma,
                                gamma = gamma,
                                areSigma0Known = areSigma0Known,
                                initialSizePerArm = initialSizePerArm,
                                design = design,
                                parDesign = parDesign),
                       simplify = FALSE)
    
    summaries = sapply(trials, 
                       function(s) computeTrialSummaries(s)$postProb_1bestVS2best)

  }
  
  stopCluster(cl)
  
  TSbayes_Null = do.call(rbind, bayesNull)  # rows are scenarios, cols are trial replicates
  
  # Strong type-I error control
  
  eta = apply(TSbayes_Null, 1, function(x) quantile(x, probs = 1 - alphaTarg))
  etaMax = max(eta)
  cvStrongControl = list(eta = eta,
                         etaMax = etaMax,
                         alphaSingles = rowMeans(TSbayes_Null > etaMax))
  
  # Mean type-I error control
  
  etaGrid = seq(from = max(apply(TSbayes_Null, 1,
                                 function(x) quantile(x, probs = 1 - alphaMaxSingle))),
                to = max(TSbayes_Null),
                length.out = 2000 #this allows to have a fine grid
                )
  
  alphaMeanGrid = sapply(etaGrid,
                         function(e) mean(rowMeans(TSbayes_Null > e)))
  
  idx = which.min(abs(alphaMeanGrid - alphaTarg))
  
  cvMeanControl = list(
    etaMean = etaGrid[idx],
    alphaMean = alphaMeanGrid[idx],
    alphaSingles = rowMeans(TSbayes_Null > etaGrid[idx])
  )
  
  elapsedTime = Sys.time() - startTime
  
  # Output
  
  return(list(trialDetails = list(N = N, 
                                  K = K, 
                                  nullScenarios = nullScenarios, 
                                  gamma = gamma),
              areSigma0Known = areSigma0Known,
              initialSizePerArm = initialSizePerArm,
              design = design,
              parDesign = parDesign,
              TSbayes_Null = TSbayes_Null,
              alphaTarg = alphaTarg,
              alphaMaxSingle = alphaMaxSingle,
              cvStrongControl = cvStrongControl,
              cvMeanControl = cvMeanControl,
              elapsedTime = elapsedTime))
}


## Robust selection of the penalization parameter kappa for the WE design ----
#
# This function optimizes values of kappa that optimize either the patient benefit
# criterion or the power criterion (Section 4 of the paper, for more details).
# This robust selection is performed across a wide set of alternative scenarios (specified in muMatrix). 
#
# Inputs:
#   kappaGrid         - vector of kappa parameters to evaluate
#   muMatrix          - S x K matrix where each row represents a scenario's true arm means
#   M                 - number of simulated trial replicates per configuration
#   N, K              - trial sample size and number of arms
#   sigma0, gamma     - vectors of standard deviations and clinical targets (length K)
#   p                 - exponent for standard deviation
#   cutoffProb_FR     - single cutoff probability for the Fixed Randomization design
#   cutoffProb_RAD    - vector of cutoff probabilities (must match length of kappaGrid)
#   xi                - target proportion of scenarios where power rule must hold (default = 0.9)
#   areSigma0Known    - logical; if TRUE, sigma0 is treated as known
#   initialSizePerArm - burn-in observations per arm
#   seed              - RNG seed for reproducibility
#   ncores            - number of CPU cores to use (defaults to available cores - 1)
#   printProgress     - logical; if TRUE, prints progress logs to the console
#
# Output:
#   A list containing the optimal kappa parameters per each criterion, the raw performance metrics, 
#   and few other additional information.

selectKappa = function(kappaGrid,
                       muMatrix,
                       M = 5e2,
                       N,
                       K,
                       sigma0,
                       gamma,
                       p,
                       cutoffProb_FR,
                       cutoffProb_RAD,
                       xi = 0.9,
                       areSigma0Known = TRUE,
                       initialSizePerArm = 2,
                       ncores = NULL,
                       seed = 123,
                       printProgress = TRUE) {
  
  if(length(cutoffProb_RAD)!=length(kappaGrid)){
    stop("You need to provide a cut-off probability for each value of kappa in the grid.")
  }
  
  overallStartTime = Sys.time() 
  
  S = nrow(muMatrix)
  nKappa = length(kappaGrid)
  
  # Define consistent naming for dimensions
  scenarioNames = paste0("Scenario_", 1:S)
  kappaNames = paste0("kappa_", kappaGrid)
  
  ## =================================================================================
  ## STEP 1: SIMULATION AND OPERATING CHARACTERISTICS FOR FR (BENCHMARK) AND WE DESIGN
  ## =================================================================================
  
  # Parallel setup
  
  require(doParallel)
  require(doRNG)
  
  if (is.null(ncores)) {
    ncores = detectCores() - 1
  }
  
  cl = makeCluster(ncores, outfile = "")
  registerDoParallel(cl)
  
  if (printProgress) {
    message("Running simulations across ", S, " scenarios using ", ncores, " cores...")
  }
  
  # Computing operating characteristics
  
  out = foreach(s = 1:S, 
                .export = c("computeTrialSummaries", "runTrial", "prob.1bestVS2best", "informationGain", "runMultipleTrials", "computeOCs"),
                .inorder = TRUE) %dorng% {
    
    mu_s = muMatrix[s, ]
    trueBestArm = which.min(abs(mu_s - gamma))
    
    # Run Fixed Randomization (FR) benchmark for this scenario
    FR_trials = runMultipleTrials(
      M = M, N = N, K = K, mu = mu_s, sigma0 = sigma0, gamma = gamma,
      areSigma0Known = areSigma0Known, initialSizePerArm = initialSizePerArm,
      design = "FR", parDesign = NULL, seed = seed, printTimeElaps = FALSE)
    
    FR_powerBenchmark = computeOCs(trials = FR_trials$trials, cutoffProb = cutoffProb_FR)$powerTwoComp
    
    # Run RAD across the grid of values for kappa for this scenario
    patientBenefit = numeric(nKappa)
    power = numeric(nKappa)
    
    for (i in seq_along(kappaGrid)) {
      RAD_trials = runMultipleTrials(
        M = M, N = N, K = K, mu = mu_s, sigma0 = sigma0, gamma = gamma,
        areSigma0Known = areSigma0Known, initialSizePerArm = initialSizePerArm,
        design = "RAD", parDesign = list(kappa = kappaGrid[i], p = p), 
        seed = seed, printTimeElaps = FALSE)
      
      RAD_OC = computeOCs(trials = RAD_trials$trials, cutoffProb = cutoffProb_RAD[i])
      
      # Extract metrics 
      patientBenefit[i] = RAD_OC$armPerc_average[trueBestArm]
      power[i] = RAD_OC$powerTwoComp
    }
    
    return(list(FR_powerBenchmark = FR_powerBenchmark,
                patientBenefit = patientBenefit,
                power = power))
    
  }
  
  # Re-assemble matrices
  
  FR_powerBenchmarkVec = sapply(out, function(x) x$FR_powerBenchmark)
  patientBenefitMat = t(sapply(out, function(x) x$patientBenefit))
  powerMat = t(sapply(out, function(x) x$power))
  
  ## ==========================================
  ## STEP 2: COMPUTE OBJECTIVE FUNCTIONS
  ## ==========================================
  
  if (printProgress) {
    message("\nSimulations complete. Optimizing objective functions...")
  }
  
  # Patient Benefit 
  xBest = apply(patientBenefitMat, 1, max)
  g_percExp = colMeans( sweep(patientBenefitMat, 1, xBest, "-")^2 )
  
  # Power metric
  probRADvs80FR = colMeans(sweep(powerMat, 1, 0.8 * FR_powerBenchmarkVec, ">"))
  
  # Record the results in a table
  resultsTable = data.frame(kappa = kappaGrid,
                            g_percExp = g_percExp,
                            probRADvs80FR = probRADvs80FR)
  
  ## ==========================================
  ## STEP 3: SELECT BEST KAPPA PER CRITERION
  ## ==========================================
  
  # Patient benefit criterion: kappa minimizing g_percExp
  
  optimalKappa_PB = kappaGrid[which.min(g_percExp)]
  
  # Power criterion: minimum kappa where RAD > 0.8*FR in at least (xi*100)% of scenarios
  
  validKappas = resultsTable$kappa[resultsTable$probRADvs80FR >= xi]
  optimalKappa_power = if(length(validKappas) > 0) min(validKappas) else NA
  
  totalTime = Sys.time() - overallStartTime
  if (printProgress) {
    message("Execution finished. Total time: ", round(totalTime, 2), " ", attr(totalTime, "units"), "\n")
  }
  
  return(list(
    optimalKappa_PB = optimalKappa_PB,
    optimalKappa_power = optimalKappa_power,
    evaluationTable = resultsTable,
    FR_powerBenchmarkVec = FR_powerBenchmarkVec,
    patientBenefitMatrix = patientBenefitMat,
    powerMat = powerMat,
    totalTime = totalTime
  ))
}
