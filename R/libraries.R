
ncores<-4
options(BiocManager::repositories())
Check_packages<-function(packlist){
  library(utils)
  pack_not_available <- packlist[!(packlist %in% installed.packages()[,"Package"])]
  if(length(pack_not_available)>0) {

    install.packages("data/methyLImp.zip")
    if (!require("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("methylclock",force=TRUE)
    BiocManager::install("methylclockData",force=TRUE)
    BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19",force=TRUE)
    BiocManager::install("minfi",force=TRUE)
    BiocManager::install("impute",force=TRUE)
    BiocManager::install("methyLImp2",force=TRUE)
  }

}

Check_packages(c(
  "BiocManager",
  "IlluminaHumanMethylation450kanno.ilmn12.hg19",
  "methylclock",
  "methylclockData",
  "methyLImp",
  "methyLImp2",
  "impute",
  "minfi"
))

library(ggplot2)
library(Metrics)
library(missMethods)
library(parallel)
library(pracma)
library(reshape2)
library(scales)
library(tictoc)
library(BiocManager)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(methylclock)
library(methylclockData)
library(methyLImp)
library(methyLImp2)
library(impute)
library(minfi)