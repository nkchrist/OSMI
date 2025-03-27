
#Generate synthetic data and apply the 3 imputation methods
data <- impute_missing_CpGs(blood_download,ncores = 1, nRows = c(400, 1000, 3000, 5000), 
                            nCols = c(1,2,5,10, 50, 100, 500, 1000, 1500), 
                            nRepeat = 10) 
