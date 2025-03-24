# OSMI

The OSMI project has been implemented for missing values imputation on expression data set. The advantage of the method is that it can be applied on a single subject DNA-methylation data. The method is based on the localisation of CpGs on chromosomes. Basically it replaces the missing values with its nearest available CpGs value.

## Installation

The project is available on github <https://github.com/Kemda-Imsid/OSMI.git> and can be installed these ways :

``` r
# install.packages("devtools")
devtools::install_github("https://github.com/Kemda-Imsid/OSMI.git") ##recommanded

library(devtools)
install_github("https://github.com/Kemda-Imsid/OSMI.git")

###You can download the .zip data from github  

#Local installation with base 
 install.packages("/Yourpath/OSMI.zip",  repos = NULL, type = 'source')
# Local installation with devtools
 library(devtools)
devtools::install("/Yourpath/OSMI.zip")
```

## Main functions

``` r
Beside the main fucntion OSMI, others has been writen for an specific example of the purpose. 

OSMI: A function to impute the missing values on expression data at single and multiples samples data set.

impute_missing_CpGs: A function to impute missing values on expression data using impute.knn, methyLImp and OSMI. 

eliminateMissingsMulti: a function to eliminates stepwise rows and columns comtaining missing values

sample_generator: a function to subset a large data set into a given smaller size

OSMI_on_aging_clock
```

## 
