#### Preparation ####

#install packages
renv::restore() #Press Y when prompted

#load packages
library(devtools)
library(BiocManager)
library(ggplot2)
library(impute)
library(methyLImp)
library(methyLImp2)
library(Metrics)
library(missMethods)
library(parallel)
library(pracma)
library(reshape2)
library(tictoc)
library(minfi)
library(scales)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(methylclockData)
library(methylclock)

# Set a longer timeout value for the download (in seconds)
options(timeout = 3600)  # 1 hour timeout

#Download data
#DNA methylation profiles 25 blood cell types 
download.file(url = "https://download.cncb.ac.cn/ewas/datahub/download/blood_methylation_v1.zip", 
              destfile = "blood_methylation_v1.zip", mode = "wb")

#Unpack data
unzip("blood_methylation_v1.zip", exdir = "blood_methylation_data")

#remove zip folder
file.remove("blood_methylation_v1.zip")

# load custom functions
scripts <- paste0("R/", dir("R"))
invisible({lapply(scripts, source)})
rm(scripts)

#load data into R
load("blood_methylation_data/blood_methylation_v1.RData")
