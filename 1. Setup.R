# =================== Package and Function Setup

# for data access
library(ETLUtils)
library(ffbase)
library(ff)
options("fftempdir"="D:/Data/Temp")
options(scipen=999)

# for data wrangling
library(doBy)
library(tidyr)
library(dplyr)
library(data.table)
library(xlsx)
library(stringr)

# for analyses
library(car)
library(tables)
library(Hmisc)
library(gmodels)
library(caret)
library(PerformanceAnalytics)
library(truncdist)
library(foreach)
library(doParallel)
library(bigstatsr)

# for modelling
library(survival)
library(markovchain)
library(fitdistrplus)
library(actuar)
library(goftest)

#for plots
library(ggplot2)
library(scales)
library(ggthemes)
library(extrafont)
library(RColorBrewer)
library(survminer)
library(ggfortify)
library(corrplot)
library(VennDiagram)

#font_import()

# -- Ensure critical functions are defined
source("DelinqM.R")
source("3.1a Function definitions for loss assesment.R")
source("3.1b Function definitions for forecasting.R")

# custom mode function
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

# recursive function that tests for non-zero interest rate within history of a given account
# if zero, then recursively test previous value.
getIntRate <- function(rateHist, period) {
  if (rateHist[period] == 0) {
    return (getIntRate(rateHist, period-1))
  } else {
    return(rateHist[period])
  }
}