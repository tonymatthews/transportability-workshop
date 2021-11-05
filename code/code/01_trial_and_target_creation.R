###############################################################################
# Author: Alex Keil
# Program: 01_trial_and_target_creation.R
# Language: R
# Date: January 8, 2020
# Project: An introduction to transporting treatment effects from randomized
#  clinical trials to clinical practice
# Description: Create trial and target datasets
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


# import and clean trial data
# limit to complete cases
analytictrialcohort <- read_sas(paste0(cohort, "adsl_pds2019.sas7bdat")) %>%
  merge(read_sas(paste0(cohort, "biomark_pds2019.sas7bdat")), by = "SUBJID") %>%
  rename_all(.funs = list(tolower)) %>%
  mutate(
    intrial = 1,
    wildkras = case_when(
      tolower(bmmtr1) %in% c("", "failure") ~ as.numeric(NA),
      tolower(bmmtr1) == "wild-type" ~ 1,
      TRUE ~ 0
    ),
    treatment = case_when(
      tolower(trt) == "folfox alone" ~ 0,
      TRUE ~ 1
    ),
    colon = case_when(
      diagtype == "Colon" ~ 1,
      TRUE ~ 0
    ),
    female = case_when(
      sex == "Female" ~ 1,
      TRUE ~ 0
    ),
    livermets = case_when(
      tolower(livermet) == "y" ~ 1,
      TRUE ~ 0
    )
  )  %>%
  select("intrial", "treatment", "wildkras", "female", "colon", "livermets", "age",
         "pfsdycr", "pfscr") %>%
  filter_all(complete.cases)


# import and clean target data, including simulation of KRAS for example
# limit to complete cases
set.seed(122355)
analytictargcohort <- read_sas(paste0(targ, "basepop.sas7bdat")) %>%
  rename_all(.funs = list(tolower)) %>%
  mutate(
    intrial = 0
  ) %>%
  select("intrial", "wildkras", "female", "colon", "livermets", "age") %>%
  filter_all(complete.cases)

message("Finished creating  'analytictrialcohort' and 'analytictargcohort'")