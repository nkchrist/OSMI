# OSMI

The OSMI project has been implemented for missing values imputation on expression data set. The advantage of the method is that it can be applied on a single subject DNA-methylation data. The method is based on the localisation of CpGs on chromosomes. Basically it replaces the missing values with its nearest available CpGs value.

The code in this repository can be used to reproduce the results presented in the manuscript *One-sample missing DNA-methylation value imputation* submitted to BMC Bioinformatics.

## How to run the code
Clone the repository to your local machine and run the code by sourcing the three following scripts consecutively:

### Step 1: Preparation
The script *01_preparation.R* installs the packages needed by creating a project specific library of R packages using the `renv` package. Then it loads the required packages and downloads and unpacks the DNA methylation profiles data from the web. The data will be saved to the project directory in the folder *blood_methylation_data*. Finally the script loads all custom defined functions from the *R/* directory.

**Note** If you rerun the script because you started a new R session and you want to load the packages, functions and data again, consider commenting out lines 28-37 to avoid unnecessary duplicate downloads of the dataset from the web.

### Step 2: Generation of synthetic data and imputation
After you've run the first script, run the script *02_imputation.R*. This script takes the dataset that has been loaded into the working space by the first script (`blood_download`) and uses it to generate synthetic data and to perform the 3 imputation methods discussed in the paper. Depending on the resources of your machine, this script may run for several days. If you are on a Unix machine, you can improve runtime by setting `ncores = 4` or higher for parallelisation. Note that this will not work on a windows machine! 

The function `impute_missing_CpGs()` will create a folder called *results* in you project directory, which will contain 10 subfolders, one for each repetition. Do not change the name of this directory or of any of the files in there!

### Step 3: Analysis
In the final script *03_analysis.R*, the results created in the step 2 are loaded from the directory and aggregated. Then, the plots for Figures 5-7 in the main manuscript and Figure 1 from the supplement are created consecutively. 




