###############################################################################
# Author: Alex Keil
# Program: VITT_plot_code.R
# Language: R
# Date: January 8, 2020
# Project: An introduction to transporting treatment effects from randomized
#  clinical trials to clinical practice
# Description: Visualization #2: VITT plot (variable importance)
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

# Main plotting function: most of the calculations are done in the helper function .lnor_vs_lnor_onesample
vitt_plot <- function(trialdata,       # (name) Name of trial / study population file
                      targetdata,      # (name) Name of target population file
                      modelparams = NULL, # (character vector) Set of variables to be plotted in the figure - should be continuous
                      treatmentinoutcome = TRUE, # (logical scalar) TRUE Specify whether treatment will be included in the outcome model or not - may want to examine both
                      nicenames = NULL      # (character vector or NULL) Names for figure legend (equal to number of binary + continuous covariates + #levels of categorical variables - 1)
) {
  # vitt_plot
  # This function works a little differently from the SAS version
  
  ## First, replicate the trial and target data 500 times to get a sense of uncertainty in estimates  ##
  
  ## Next, add back in the original trial and target data as replicate 0 ##
  ## Then, use the trial data to estimate multivariable odds ratios between variables and the outcome ##
  ## Next, we concatenate the study / trial and target data for estimating the sampling OR ##
  lnor_original <- .lnor_vs_lnor_onesample(trialdata, targetdata, treatmentinoutcome, modelparams)
  or_original <- exp(lnor_original)
  ## And we run another logistic regression, this time on the intrial variable ##
  B <- 500 # number of bootstrap samples
  set.seed(1231)
  seeds <- sample(seq_len(.Machine$integer.max), B + 1) # create random seed values, rather than sequential
  trial_n <- nrow(trialdata)
  target_n <- nrow(targetdata)
  gen_ors <- function(i, .trial_n = trial_n, .target_n = target_n) {
    set.seed(seeds[i])
    # bootstrap sample
    trial_boot <- trialdata[sample(seq_len(.trial_n), .trial_n, replace = TRUE), ]
    target_boot <- targetdata[sample(seq_len(.target_n), .target_n, replace = TRUE), ]
    # get log - odds ratios
    lnors <- .lnor_vs_lnor_onesample(trial_boot, target_boot, treatmentinoutcome, modelparams)
    ors <- exp(lnors)
    ors
  }
  # repeat for total number of bootstrap samples, using parallel processing
  # # note that this will run if run interactively (e.g. in Rstudio) but would
  reslist <- future.apply::future_lapply(seq_len(B), gen_ors, .trial_n = trial_n, .target_n = target_n, future.seed = seeds[B + 1])
  res <- data.frame(do.call("rbind", reslist))
  res2 <- data.frame(rbind(or_original, res)) # include original results
  vars <- gsub("sel_", "", grep("sel_", names(res), value = TRUE))
  vals <- gsub("\\.", " = ", gsub("factor\\.", "", vars))
  if (!is.null(nicenames)) {
    if (length(nicenames) !=  length(vals)) stop(paste0("nicenames is the wrong length (should have ", length(vals), " values)"))
    vals <- nicenames
  }
  print(vals)
  print(vars)
  p <- ggplot() + theme_classic() +
    scale_color_manual(name = "", breaks = waiver(), labels = vals, limits = vars, values = 1:10) +
    scale_shape_manual(name = "", breaks = c("Bootstrap draw", "Estimate"), values = c(20, 18)) +
    scale_fill_manual(name = "", breaks = waiver(), labels = vals, limits = vars, values = 1:10) +
    scale_y_log10(name = "Outcome odds ratios") +
    scale_x_log10(name = "'Sampling' odds ratios")
  for (v in vars) {
    p <- p + geom_point(aes_(x = as.name(paste0("sel_", v)), y = as.name(paste0("out_", v)), fill = v, color = v, shape = "Bootstrap draw"), stroke = 0, alpha = 0.25, data = res)
  }
  for (v in vars) {
    p <- p + geom_point(aes_(x = as.name(paste0("sel_", v)), y = as.name(paste0("out_", v)), shape = "Estimate"),
                        data = res2[1, , drop = FALSE], color = "black")
  }
  p <- p + geom_hline(aes(yintercept = 1), size = .3) + geom_vline(aes(xintercept = 1), size = .3)
  p[["data"]] <- res2
  p
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


.lnor_vs_lnor_onesample <- function(trialdata,       # (name) Name of trial / study population file
                                    targetdata,      # (name) Name of target population file
                                    treatmentinoutcome = TRUE, # (logical scalar) TRUE Specify whether treatment will be included in the outcome model or not - may want to examine both
                                    modelparams = NULL    # (character vector) Set of variables to be plotted in the figure - should be continuous or binary
) {
  # selection regression
  outcome <- "intrial"
  concat <- .joindata(trialdata, targetdata)
  # build model formulas in R syntax
  f_sel_str <- paste0(outcome, " ~ ", paste(c("1", modelparams), collapse = " + "))
  f_sel <- as.formula(f_sel_str)
  m_sel <- glm(formula = f_sel, data = concat, family = binomial())
  #
  # outcome regression
  outcome <- "outcome"
  if (treatmentinoutcome) modelparams <- c(modelparams, "treatment")
  # build model formulas in R syntax
  f_out_str <- paste0(outcome, " ~ ", paste(c("1", modelparams), collapse = " + "))
  f_out <- as.formula(f_out_str)
  m_out <- glm(formula = f_out, data = trialdata, family = binomial())
  #
  # output
  ff <- c(coef(m_out)[-1], coef(m_sel)[-1])
  names(ff) <- c(paste0("out_", names(coef(m_out)[-1])), paste0("sel_", names(coef(m_sel)[-1])))
  ff
}


#### 3. Analyses ####
# First, let's import the trial and target data
source(paste0(prgdir, "01_trial_and_target_creation.R"))
# create local copies
trialset <- analytictrialcohort %>%
  mutate(
    outcome = case_when(
      pfsdycr >=  365 ~ 0,
      pfscr == 1 ~ 1,
      TRUE ~ as.numeric(NA)
    ),
    age = age / 20
  )
targset  <- analytictargcohort %>%
  mutate(
    age = age / 20
  )


# set up parallel processing using all available cores (takes a few seconds)
future::plan(multisession)
# create plot using the vitt_plot function
p <- vitt_plot(trialdata = trialset,
               targetdata = targset,
               modelparams = c("age", "female", "colon", "livermets", "wildkras"),
               treatmentinoutcome = TRUE,
               nicenames = c("Age / 20", "Female sex", "Colon cancer", "Liver Metastases", "Wild - type - KRAS")
)
# additional cleanup of the plot
p <- p + coord_fixed(xlim = c(0.1, 10), ylim = c(0.1, 10), ratio = 1) # note this is extensible with ggplot2 functions
p
ggsave(filename = paste0(figdir, "VITT_plot.png"), plot = p, width = unit(6, "in"), height = unit(3.5, "in"))