###############################################################################
# Author: Alex Keil
# Program: Skyscraper_plot_code.R
# Language: R
# Date: January 8, 2020
# Project: An introduction to transporting treatment effects from randomized
#  clinical trials to clinical practice
# Description: Visualization #4: A skyscraper plot of the weights resulting from the sampling model
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
library("ggplot2")

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

# Main plotting function: most of the calculations are done in the helper function .getweights
skyscraper_plot <- function(trialdata,       # (name) Name of trial / study population file
                            targetdata,      # (name) Name of target population file
                            modelparams = NULL, # (character vector) Set of variables to be plotted in the figure - should be continuous
                            figurename = NULL      # (not used in R) The name of the output data set that includes weights
) {
  withweights <- .getweights(trialset, targset,
                             modelparams = c("age", "female", "colon", "livermets", "wildkras"),
                             weightbyarm = FALSE)
  if (any(is.na(withweights$siosw))) warning("Some weights are missing, plotting function excludes these")
  withweights$forscramble <- runif(nrow(withweights))
  plotdata <- data.frame(index = seq_len(nrow(withweights)), siosw = withweights$siosw[order(withweights$forscramble)]) %>%
    filter(!is.na(siosw) & siosw > 0)
  plotdata$place <- .nomisscumsum(plotdata$siosw) # cumulative sum excluding missing values
  p <- ggplot(data = plotdata) + theme_classic() +
    scale_fill_discrete(name = "") +
    scale_x_continuous(name = "Cumulative sum of weights") +
    scale_y_continuous(name = "Stabilized weights")
  p +
    geom_line(aes(x = place, y = siosw), size = .2)
}


### helper functions (prefixed with ".")

## weight calculation (copied from transport analyses program) ##
# This function calculates inverse odds weights based on a logistic model
# that estimates the "pseudo-probability" of target membership given a
# linear parameterization. Actual weight calculation is passed off to the helper
# function .calcwts, while data concatenation is handled by the .joindata function
.getweights <- function(trialdata,          # (name) The name of the trial data set
                       targetdata,         # (name) The name of the target data set
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

# cumulative sum, excluding missing values
.nomisscumsum <- function(variable) {
  tvar <- variable
  tvar[is.na(variable)] <- 0
  cumsum(tvar)
}

# This function joins two datasets in a clean fashion that avoids some of the
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
.calcwts <- function(concat,             # (name) The name of the trial + target concatenated data set
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
         siosw = intrial * numodds / denomodds
  )
}

#### 3. Analyses ####
# First, let's import the trial and target data
source(paste0(prgdir, "01_trial_and_target_creation.R"))
# create local copies
trialset <- analytictrialcohort
targset  <- analytictargcohort

# the skyscraper plot randomly sorts the dataset, so set a seed for reproducibility
set.seed(124)
# create plot using the skyscraper_plot function
p <- skyscraper_plot(trialset, targset,
                     modelparams = c("age", "female", "colon", "livermets", "wildkras")
)

# additional cleanup of the plot
p <- p + coord_cartesian(ylim = c(0, 20)) + geom_hline(aes(yintercept = 10), size = .1) # note this is extensible with ggplot2 functions
p
ggsave(filename = paste0(figdir, "Skyscraper.png"), plot = p, width = unit(4, "in"), height = unit(3.5, "in"))
