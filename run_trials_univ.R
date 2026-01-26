source("WE_main_functions_univ.R") #load main functions 

## ================================================================================
## Run one trial replicate for each design and analyze the output (univariate case)
## ================================================================================

## General trial settings

N = 100                 # total sample size
K = 4                   # number of arms

mu = c(1.13,-3.48,-3.57,0.34)   # true means
sigma0 = c(2,2,2,4)    # true standard deviatins
gamma = rep(0,4)  # target value (always the same for all arms)

idxBest = which.min(abs(mu-gamma)) # index of the best arm

areSigma0Known = TRUE
initialSizePerArm = 5 #burn-is size (can be changed)

## Design-specific settings

parDesign_RAD = list(  # parameters can be changed 
  kappa = 0.7,
  p = 2,
  out = list()
)

parDesign_GI_symm = list(
  GItype = "symm",
  d = 0.99,
  out = list()
)

parDesign_GI_targ = list(
  GItype = "targ",
  d = 0.99,
  out = list()
)

parDesign_TS = list(
  MCiter = 2e3, #HIGH VALUES ARE TIME CONSUMING - default in paper simulations is 2e4
  random = TRUE,
  out = list()
)

## Run single trials

results = list()
seed = 123

## ---- FR ----
set.seed(seed)

trial_FR = runTrial(
  N = N,
  K = K,
  mu = mu,
  sigma0 = sigma0,
  gamma = gamma,
  areSigma0Known = areSigma0Known,
  initialSizePerArm = initialSizePerArm,
  design = "FR"
)

summary_FR = computeTrialSummaries(trial_FR)

results$FR = list(
  trial = trial_FR,
  summary = summary_FR
)

## ---- CB ----
set.seed(seed)

trial_CB = runTrial(
  N = N,
  K = K,
  mu = mu,
  sigma0 = sigma0,
  gamma = gamma,
  areSigma0Known = areSigma0Known,
  initialSizePerArm = initialSizePerArm,
  design = "CB"
)

summary_CB = computeTrialSummaries(trial_CB)

results$CB = list(
  trial = trial_CB,
  summary = summary_CB
)

## ---- RAD ----
set.seed(seed)

trial_RAD = runTrial(
  N = N,
  K = K,
  mu = mu,
  sigma0 = sigma0,
  gamma = gamma,
  areSigma0Known = areSigma0Known,
  initialSizePerArm = initialSizePerArm,
  design = "RAD",
  parDesign = parDesign_RAD
)

summary_RAD = computeTrialSummaries(trial_RAD)

results$RAD = list(
  trial = trial_RAD,
  summary = summary_RAD
)

## ---- symm GI ----
set.seed(seed)

trial_GI_symm = runTrial(
  N = N,
  K = K,
  mu = mu,
  sigma0 = sigma0,
  gamma = gamma,
  areSigma0Known = areSigma0Known,
  initialSizePerArm = initialSizePerArm,
  design = "GI",
  parDesign = parDesign_GI_symm
)

summary_GI_symm = computeTrialSummaries(trial_GI_symm)

results$GI_symm = list(
  trial = trial_GI_symm,
  summary = summary_GI_symm
)

## ---- targ GI ----
set.seed(seed)

trial_GI_targ = runTrial(
  N = N,
  K = K,
  mu = mu,
  sigma0 = sigma0,
  gamma = gamma,
  areSigma0Known = areSigma0Known,
  initialSizePerArm = initialSizePerArm,
  design = "GI",
  parDesign = parDesign_GI_targ
)

summary_GI_targ = computeTrialSummaries(trial_GI_targ)

results$GI_targ = list(
  trial = trial_GI_targ,
  summary = summary_GI_targ
)


## ---- TS ----
set.seed(seed)

trial_TS = runTrial(
  N = N,
  K = K,
  mu = mu,
  sigma0 = sigma0,
  gamma = gamma,
  areSigma0Known = areSigma0Known,
  initialSizePerArm = initialSizePerArm,
  design = "TS",
  parDesign = parDesign_TS
)

summary_TS = computeTrialSummaries(trial_TS)

results$TS = list(
  trial = trial_TS,
  summary = summary_TS
)

## Quick comparison on single trials (table format)


sapply(results, function(res) {
    
    idxBest = res$trial$trialDetails$trueBestArm #index best arm
    nBest = res$summary$nByArm[idxBest] #size best arm
    postMeanBest = res$summary$aveResponseByArm[idxBest] 
    postSDBest = res$summary$sdResponseByArm[idxBest] / sqrt(nBest) 
    
    round(c(correctBest = res$summary$correctSelectionBest, #is best arm correctly selected?
            correctTop2 = res$summary$correctSelectionTop2, #are top-2 arms correctly selected?
            propAllocationBest = res$summary$allocationProportions[idxBest], #proportion of participants allocated to best arm
            postMeanBest = postMeanBest, #posterior mean best arm
            postSDBest = postSDBest, #posterior standard deviation best arm
            post95CIBest_lb = qnorm(0.025, mean = postMeanBest, sd = postSDBest), #95% posterior CI lower bound
            post95CIBest_ub = qnorm(0.975, mean = postMeanBest, sd = postSDBest) #95% posterior CI upper bound
            ), 2)
    })


## =====================================================================================================
## Run multiple trial replicates for each design and compute operating characteristics (univariate case)
## =====================================================================================================

M = 1e2 # increase for higher accuracy

parDesign_list = list(FR = NULL,
                      CB = NULL,
                      TS = parDesign_TS,
                      GI_symm = parDesign_GI_symm,
                      GI_targ = parDesign_GI_targ,
                      RAD = parDesign_RAD)

cutoffProb = c(FR = 0.917, 
               CB = 0.912, 
               TS = 0.912, 
               GI_symm = 0.912, 
               GI_targ = 0.913, 
               RAD = 0.914)

## Run multiple trial replicates

multiRuns = lapply(names(parDesign_list), function(des) {
    
  runMultipleTrials(M = M,
                    N = N,
                    K = K,
                    mu = mu, sigma0 = sigma0, gamma = gamma,
                    areSigma0Known = TRUE, initialSizePerArm = 5,
                    design = ifelse(des %in% c("GI_symm", "GI_targ"), "GI", des),
                    parDesign = parDesign_list[[des]],
                    seed = seed)
  })

names(multiRuns) = names(parDesign_list)

## Compute operating characteristics

OCs = lapply(names(multiRuns), 
             function(des) {
  
               computeOCs(trials = multiRuns[[des]]$trials,
                          cutoffProb = cutoffProb[des])
               })

names(OCs) = names(multiRuns)


OC_table = t(sapply(names(OCs), 
                    function(des) {
                      
                      foo = c(percAllocBest = OCs[[des]]$armPerc_average[idxBest]/100,
                              percAllocBest_SE = OCs[[des]]$armPerc_SE[idxBest],
                              correctBest = OCs[[des]]$percSelectionBest,
                              correctTop2 = OCs[[des]]$percSelectionTop2,
                              powerCond = as.numeric(OCs[[des]]$powerCond),
                              powerTwoComp = as.numeric(OCs[[des]]$powerTwoComp))
                      
                      round(foo, 3)
                      }
                    )
             )

print(OC_table)



## ===================================================================
## Compute critical values given the set of null scenarios of interest
## ===================================================================

library(doParallel)
library(doRNG)

# Simple trial setup
K = 4                     
N = 100                    
M = 1e3         # small value for quick run (this should be ideally increase)
gamma = rep(0, K)    
sigma0 = c(rep(2,K-1), 4)

# Define null configurations (all mu's equal)
upperLimit = 10 * max(sigma0)
grid = seq(0, sqrt(upperLimit), length = 15)^2  # spacing on squared scale
muNull = matrix(rep(grid, K), ncol = K, byrow = FALSE)

# Convert to list of null scenarios
nullScenarios = lapply(1:nrow(muNull), function(i) list(mu = muNull[i, ], sigma = sigma0))

ncores = min(length(nullScenarios), detectCores() - 1)

# Run critical value search (FR design)
critValue_FR = criticalValueSearchBayes(M = M, K = K, N = N,
                                        nullScenarios = nullScenarios, gamma = gamma,
                                        areSigma0Known = TRUE, initialSizePerArm = 5,
                                        design = "FR", parDesign = NULL,
                                        alphaTarg = 0.05,
                                        ncores = ncores,
                                        seed = 123
)


## Inspection of results

# cut-off probability values

critValue_FR$cvStrongControl$etaMax # strong control of TIE at 5%
critValue_FR$cvMeanControl$etaMean # mean control of TIE at 5%

# scenario-specific TIE rates with strong control at 5%
critValue_FR$cvStrongControl$alphaSingles   

# scenario-specific TIE rates with mean control at 5%
critValue_FR$cvMeanControl$alphaSingles   

# prob.1bestVS2best (piHat) under H0: mu_j-gamma=0, for all j
hist(critValue_FR$TSbayes_Null[1,])  
