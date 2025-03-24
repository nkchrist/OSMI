
ncores<-4
#options(BiocManager::repositories())
# Check_packages<-function(packlist){
#   library(utils)
#   pack_not_available <- packlist[!(packlist %in% installed.packages()[,"Package"])]
#   if(length(pack_not_available)>0) {
#
#     install.packages("data/methyLImp.zip")
#     if (!require("BiocManager", quietly = TRUE))
#       install.packages("BiocManager")
#     BiocManager::install("methylclock",force=TRUE)
#     BiocManager::install("methylclockData",force=TRUE)
#     BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19",force=TRUE)
#     BiocManager::install("minfi",force=TRUE)
#     BiocManager::install("impute",force=TRUE)
#     BiocManager::install("methyLImp2",force=TRUE)
#   }
#
# }
#
# Check_packages(c(
#     "IlluminaHumanMethylation450kanno.ilmn12.hg19",
#     "impute",
#     "methylclock",
#     "methylclockData",
#     "methyLImp",
#     "methyLImp2",
#     "minfi"
# ))
#
# library(methyLImp)
# library(methyLImp2)
# library(minfi)
# library(methylclock)
# library(methylclockData)
# library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
# library(impute)
#
#
# # use_package("IlluminaHumanMethylation450kanno.ilmn12.hg19")
# # use_package("impute")
# # use_package("methylclock")
# # use_package("methylclockData")
# # use_package("methyLImp")
# # use_package("methyLImp2")
# # use_package("minfi")
