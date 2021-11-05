###############################################################################
# Author: Alex Keil
# Program: Love_plot_code.R
# Language: R
# Date: January 8, 2020
# Project: An introduction to transporting treatment effects from randomized
#  clinical trials to clinical practice
# Description: Visualization #1: Love plot
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

# Main plotting function: most of the calculations are done in the helper function .calculate_smds
love_plot <- function(trialdata, # Name of trial / study population file
                      targetdata, # Name of target population file
                      variables, # names of variables used in loveplot (binary or continuous assumed)
                      nicenames = NULL # Names of variables to be printed on plot
) {
  if (is.null(nicenames)) nicenames <- variables
  smds <- .calculate_smds(trialdata, targetdata, variables)
  plotdata <- cbind(smds, data.frame(variable = variables, nicenames = nicenames))
  plotdata
  p <- ggplot(data = plotdata) + theme_classic() +
    geom_vline(aes(xintercept = c(0.0))) +
    geom_vline(aes(xintercept = c(-.1)), linetype = "dashed", color = "gray50") +
    geom_vline(aes(xintercept = c(.1)), linetype = "dashed", color = "gray50") +
    geom_point(aes(x = smd, y = variable)) +
    scale_y_discrete(name = "Covariate", labels = plotdata$nicenames, limits = plotdata$variable) +
    scale_x_continuous(name = "Standardized Mean Difference")
  p[["data"]] <- plotdata
  p
}




### helper functions (prefixed with ".")

# calculate standardized mean difference for all variables
.calculate_smds <- function(trialdata, # Name of trial / study population file
                            targetdata, # Name of target population file
                            variables # names of variables used in loveplot (binary or continuous assumed)
) {
  trialdata$intrial <- 1
  targetdata$intrial <- 0
  # combine datasets
  concat <- .joindata(trialdata, targetdata)
  # get standardized mean differenc for each variable
  smdsl <- lapply(variables, function(x) .get_smd(concat[[x]], concat$intrial, .is_this_binary(concat[[x]])))
  smds <- data.frame(do.call(rbind, smdsl))
  anymissing <- do.call(c, lapply(smdsl, function(x) attr(x, "anymissing")))
  if (any(anymissing)) {
    warning(paste("Missing / NA values (ignored in computations) detected in:", paste(variables[which(anymissing)], collapse = ", ")))
  }
  smds
}

# calculate standardized mean difference for a single variable
.get_smd <- function(v, intrial, isbinary = FALSE) {
  # warn if there are missing data (ideally we use analytic, complete - case data here)
  set_warning <- any(is.na(v))
  # 1. calculate mean of each variable
  meandiff <- mean(v[which(intrial == 1)], na.rm = TRUE) - mean(v[which(intrial == 0)], na.rm = TRUE)
  # 2. calculate the mean difference
  pooledsd <- .sd_pooled(v, intrial, isbinary) # equal weight to all data sets
  # 3. divide the difference by pooled standard deviation
  r <- c(smd = meandiff / pooledsd, meandiff = meandiff, poolesd = pooledsd)
  attr(r, which = "anymissing") <- set_warning
  r
}

# figure out whether a variable is binary or not
.is_this_binary <- function(v) {
  length(unique(v)) == 2
}

# calculate the pooled standard deviation
.sd_pooled <- function(v, intrial, isbinary = FALSE) {
  # following https: /  / www.ncbi.nlm.nih.gov / pmc / articles / PMC3472075 /
  # this applies equal weight to target and trial data (but not equal weight to all observations)
  # (not the only way to calculate a pooled SD)
  if (!isbinary) {
    var1 <- var(v[which(intrial == 1)], na.rm = TRUE)
    var0 <- var(v[which(intrial == 0)], na.rm = TRUE)
  }
  if (isbinary) {
    p1 <- mean(v[which(intrial == 1)], na.rm = TRUE)
    p0 <- mean(v[which(intrial == 0)], na.rm = TRUE)
    var1 <- p1 * (1 - p1)
    var0 <- p0 * (1 - p0)
    # FLAG
  }
  sqrt((var1 + var0) / 2)
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

# create plot using the love_plot function
p <- love_plot(trialdata = trialset,
               targetdata = targset,
               variables = c("age", "livermets", "wildkras", "female", "colon"),
               nicenames = c("Age (in years)", "Liver metastases", "Wild - type KRAS", "Female sex", "Colon cancer")
)
# additional cleanup of the plot
p <- p + coord_cartesian(xlim = c(-.75, .75)) # note this is extensible with ggplot2 functions
p
ggsave(filename = paste0(figdir, "Love_plot.png"), plot = p, width = unit(5, "in"), height = unit(3.5, "in"))
