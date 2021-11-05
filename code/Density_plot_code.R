###############################################################################
# Author: Alex Keil
# Program: Density_plot_code.R
# Language: R
# Date: January 8, 2020
# Project: An introduction to transporting treatment effects from randomized
#  clinical trials to clinical practice
# Description: Visualization #3: A density plots of the predicted "sampling" probabilities and their overlap
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

# Main plotting function: most of the calculations are done in the helper function .calcprobs
dens_plot <- function(trialdata,       # (name) Name of trial / study population file
                      targetdata,      # (name) Name of target population file
                      catparams = NULL, # (character vector) Set of variables to be plotted in the figure - should be binary
                      modelparams = NULL, # (character vector) Set of variables to be plotted in the figure - should be continuous
                      figurename = NULL      # (not used in R) The name of the output data set that includes weights
) {
  concat <- .joindata(trialdata, targetdata)
  res <- .calcprobs(concat, modelparams)
  p <- ggplot() + theme_classic() +
    scale_fill_discrete(name = "") +
    scale_x_continuous(name = "Predicted probability of presence in the trial data") +
    scale_y_continuous(name = "Percent")
  p <- p +
    geom_histogram(aes(x = trialprobs, fill = "Trial", y = 100 * stat(count) / sum(count)), alpha = 0.5,
                   data = filter(res, intrial == 1 & !is.na(trialprobs)), binwidth = 0.005, color = "black", size = 0.1) +
    geom_histogram(aes(x = targprobs, fill = "Target", y = 100 * stat(count) / sum(count)), alpha = 0.5,
                   data = filter(res, intrial == 0), binwidth = 0.005, color = "black", size = 0.1)
  p[["data"]] <- res
  p
}



### helper functions (prefixed with ".")

# This function calculates inverse odds weights based on a logistic model
# that estimates the "pseudo-probability" of target membership given a
# linear parameterization with possible indicator functions for categorical
# variables
.calcprobs <- function(concat,           # (name) The name of the trial + target concatenated data set
                       modelparams = NULL # (character vector) All the variables to be used in the model
) {
  # build model formula in R syntax
  outcome <- "intrial"
  f_samp_str <- paste0(outcome, " ~ ", paste(c("1", modelparams), collapse = " + "))
  f_samp <- as.formula(f_samp_str)
  # fit model
  m_samp <- glm(formula = f_samp, data = concat, family = binomial())
  # get predicted probabilities
  concat$sampprob <- predict(m_samp, type = "response", newdata = concat)
  message(paste0("Sampling model: ", f_samp_str))
  mutate(.data = concat,
         trialprobs = case_when(
           intrial == 1 ~ sampprob,
           TRUE ~ as.numeric(NA)
         ),
         targprobs = case_when(
           intrial == 0 ~ sampprob,
           TRUE ~ as.numeric(NA)
         )
  )
}

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



#### 3. Analyses ####
# First, let's import the trial and target data
source(paste0(prgdir, "01_trial_and_target_creation.R"))
# create local copies
trialset <- analytictrialcohort
targset  <- analytictargcohort

# create plot using the dens_plot function
p <- dens_plot(trialdata = trialset,
               targetdata = targset,
               modelparams = c("age", "female", "colon", "livermets", "wildkras"),
               figurename = NULL
)
# additional cleanup of the plot
p <- p + coord_cartesian(xlim = c(0, .20)) # note this is extensible with ggplot2 functions
p
ggsave(filename = paste0(figdir, "Density_plot.png"), plot = p, width = unit(5, "in"), height = unit(3.5, "in"))
