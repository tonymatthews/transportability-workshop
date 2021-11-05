###############################################################################
# Author: Alex Keil
# Program: 02_transport_analyses.R
# Language: R
# Date: January 8, 2020
# Project: An introduction to transporting treatment effects from randomized
#  clinical trials to clinical practice
# Description: Running the transport analysis
# Released under the GNU General Public License:
#  http: /  / www.gnu.org / copyleft / gpl.html
###############################################################################
# structure of this program:
# 1. Preliminaries (path names, etc.)
# 2. Functions
# 3. Analysis

# NOTE: you need to modify code for your computer in locations denoted MODIFY

#### 1. Preliminaries ####
library("dplyr")
library("haven")
library("survival")
library("future") # parallel processing
library("future.apply") # parallel processing

operating_system <- Sys.info()[["sysname"]]

#MODIFY these pathnames
if (operating_system == "Windows") {
  # where all R / SAS code is stored for this project
  prgdir <- "C:/Users/antmat/OneDrive/ki/teaching/transportability/code/"
  # Naming the folder where the trial data are being kept
  cohort <- "C:/Users/antmat/OneDrive/ki/teaching/transportability/trialdata/"
  #Naming the folder where the target data are being kept
  targ <- "C:/Users/antmat/OneDrive/ki/teaching/transportability/targetdata/"
  #Naming the folder where the figures are output
  figdir <- "C:/Users/antmat/OneDrive/ki/teaching/transportability/visualizations/"
}
if (operating_system == "Darwin" || operating_system == "Linux") {
  # where all R / SAS code is stored for this project
  prgdir <- paste0(path.expand("~"), "/repo/HybridDesignVisualsWorkshop/code/")
  # Naming the folder where the trial data are being kept
  cohort <- paste0(path.expand("~"), "/repo/HybridDesignVisualsWorkshop/data/Trial data/")
  #Naming the folder where the target data are being kept
  targ <- paste0(path.expand("~"), "/repo/HybridDesignVisualsWorkshop/data/Target data/")
  #Naming the folder where the figures are output
  figdir <- paste0(path.expand("~"), "/repo/HybridDesignVisualsWorkshop/Visualizations/")
}

#### 2. Functions ####
## Note about R functions ##
# If you are unfamiliar with the creation / use of functions in R,
# the following functions do the heavy lifting in the analyses and
# can be seen as a rough equivalent to SAS macro. They help to keep
# the code relatively tidy and keep you from having to create many different
# intermediate data steps. You can often directly run the code within a function
# provided that the function arguments have values (e.g. if "trialdata" is
# a function argument, then be sure that your trial dataset is called "trialdata"
# to use the calculations inside the function). This can make the code a little
# more abstract if you are not used to reading it, so please ask questions if
# you don't understand something.


## weight calculation ##
# This function calculates inverse odds weights based on a logistic model
# that estimates the "pseudo-probability" of target membership given a
# linear parameterization with possible indicator functions for categorical
# variables: this function is a "wrapper" function that
getweights <- function(trialdata,       # (name) The name of the trial data set
                       targetdata,      # (name) The name of the target data set
                       modelparams = NULL, # (character vector) All the variables to be used in the model
                       weightbyarm = TRUE
) {
  # create joint dataset
  # IF weightbyarm is TRUE we need to split the trial data into two groups,
  # one for each arm, then combine each of them with the target data
  if (weightbyarm) {
    message("Fit weights separately by trial arm")
    trialdata$order <- seq_len(nrow(trialdata))
    concat1 <- .joindata(filter(trialdata, treatment == 1), targetdata)
    ret1 <- .calcwts(concat1, modelparams)
    concat0 <- .joindata(filter(trialdata, treatment == 0), targetdata)
    ret0 <- suppressMessages(.calcwts(concat0, modelparams))
    # Next, we combine the data sets;
    ret <- rbind(filter(ret1, intrial == 1), ret0)
    ret <- ret[order(ret$order), ]
  }else{
    #If weightbyarm is FALSE, we instead simply concatenate the trial and target populations,
    # then run similar logistic regressions
    message("Fit weights using entire trial population in a single model")
    concat <- .joindata(trialdata, targetdata)
    ret <- .calcwts(concat, modelparams)
  }
  return(ret)
}

npbootstrap <- function(iters = 5, trialdata = trialset, targetdata = targset, ..., verbose = TRUE) {
  outset <- getweights(trialdata = trialdata, targetdata = targetdata, ...)
  wfit <- coxph(Surv(pfsdycr, pfscr) ~ treatment,
                data = filter(outset, intrial == 1),
                weights = sioperating_systemw)
  ocoef <- coef(summary(wfit))[, c("coef", "robust se")]
  message(paste("Bootstrapping", iters, "iterations"))
  # single bootstrap iteration helper function
  .bootiter <- function(i, trialdata = NULL, targetdata = NULL, ..., 
                        verbose = TRUE) {
    if (verbose) cat(".") # progress bar won't work in parallel processing
    ntrial <- nrow(trialdata)
    ntarg <- nrow(targetdata)
    trialidx <- sample(seq_len(ntrial), ntrial, replace = TRUE)
    targidx <- sample(seq_len(ntarg), ntarg, replace = TRUE)
    suppressMessages(outset <- getweights(trialdata = trialdata[trialidx, ], targetdata = targetdata[targidx, ], ...))
    wfit <- coxph(Surv(pfsdycr, pfscr) ~ treatment,
                  data = filter(outset, intrial == 1),
                  weights = sioperating_systemw)
    res <- coef(summary(wfit))[, c("coef", "robust se")]
    as.numeric(res)
  }
  
  boots <- future.apply::future_lapply(seq_len(iters), .bootiter, trialdata = trialdata, 
                                       targetdata = targetdata, ..., verbose = verbose,
                                       future.seed = 56789)
  if (verbose) cat("\n")
  res <- list()
  rr <- do.call("rbind", boots)
  colnames(rr) <- c("coef", "robust se")
  res[["boots"]] <- as.data.frame(rr)
  res[["original"]] <- as.data.frame(ocoef)
  class(res) <- "transboot"
  res
}

print.transboot <- function(x, ...) {
  # this just formats how npbootstrap results are shown in output
  # it takes advantage of the "print" generic R function
  oest <- x$original
  best <- x$boots
  # print out original estimate
  cat("Estimate:\n   ")
  est <- oest[1]
  cat(est[1, 1])
  cat("\n")
  cat("Robust standard error:\n   ")
  cat(oest[2, 1])
  cat("\n")
  # bootstrap standard error
  cat("Bootstrap standard error:\n   ")
  booterr <- sd(best[, 1])
  cat(booterr)
  cat("\n")
  cat("Bootstrap Wald 95% CI:\n   ")
  cat(est[1, 1] + qnorm(c(.025, .975)) * booterr)
  cat("\n")
}


### helper functions (prefixed with ".")

# This function joints two datasets in a clean fashion that avoids some of the
# ugly messages if trying to join the two datasets naively
.joindata <- function(trialdata,       # The name of the trial data set
                      targetdata          # The name of the target data set
) {
  # remove sas labels
  for (v in names(trialdata)) {
    attr(trialdata[[deparse(as.name(v))]], "label") <- NULL
  }
  for (v in names(targetdata)) {
    attr(targetdata[[deparse(as.name(v))]], "label") <- NULL
  }
  suppressMessages(full_join(trialdata, targetdata))
}


# This function calculates inverse odds weights based on a logistic model
# that estimates the "pseudo-probability" of target membership given a
# linear parameterization with possible indicator functions for categorical
# variables
.calcwts <- function(concat,           # (name) The name of the trial + target concatenated data set
                     modelparams = NULL    # (character vector) All the variables to be used in the model
) {
  # build model formula in R syntax
  outcome <- "intrial"
  f_den_str <- paste0(outcome, " ~ ", paste(c("1", modelparams), collapse = " + "))
  #To obtain the numerator of stabilized odds weights, we have to run additional regressions
  # with an empty set of variables (intercept only model)
  f_num_str <- paste0(outcome, " ~ 1")
  f_den <- as.formula(f_den_str)
  f_num <- as.formula(f_num_str)
  # fit model
  m_den <- glm(formula = f_den, data = concat, family = binomial())
  m_num <- glm(formula = f_num, data = concat, family = binomial())
  # get predicted probabilities
  p_den <- predict(m_den, type = "response", newdata = concat)
  p_num <- predict(m_num, type = "response", newdata = concat)
  # calculate stabilized weights
  # Now, we simply calculate stabilized / unstabilized inverse odds weights
  message(paste0("Numerator model: ", f_num_str))
  message(paste0("Denominator model: ", f_den_str))
  mutate(.data = concat,
         denomodds = p_den / (1 - p_den),
         numodds = p_num / (1 - p_num),
         ioperating_systemw = intrial * 1 / denomodds,
         sioperating_systemw = intrial * numodds / denomodds
  )
}



##### 3. Analysis #####

# First, let's import the trial and target data
source(paste0(prgdir, "01_trial_and_target_creation.R"))
# create local copies
trialset <- analytictrialcohort
targset  <- analytictargcohort


# Now, we want to obtain weights so that the trial resembles the target using the "getweights" function from above
outset <- getweights(trialdata = trialset, targetdata = targset,
                     modelparams = c("age", "female", "colon", "livermets", "wildkras"),
                     weightbyarm = TRUE)

#Now that we have weights, we can use them in a weighted analysis to estimate a hazard ratio for the effect of treatment
wfit <- coxph(Surv(pfsdycr, pfscr) ~ treatment,
              data = filter(outset, intrial == 1),
              weights = sioperating_systemw)
summary(wfit)

# We can compare thoperating_systeme results to the unweighted estimate of the hazard ratio
uwfit <- coxph(Surv(pfsdycr, pfscr) ~ treatment,
               data = filter(outset, intrial == 1)
)
summary(uwfit)

# compare weighted results with non-parametric bootstrap results
# this will take a while to do 250 iterations (set to 5 by default in case you run it by accident)!
# set up parallel processing using all available cores (takes a few seconds)
future::plan(multisession)
rr <- npbootstrap(iter = 5, trialdata = trialset, targetdata = targset,
                  modelparams = c("age", "female", "colon", "livermets", "wildkras"),
                  weightbyarm = TRUE, verbose = TRUE)
rr

# you can access the individual bootstrap iterations to get percentile based CI
#  (but not really advisable unless you have 1000+ bootstrap iterations)
quantile(rr$boots$coef, probs = c(.025, 0.975))

