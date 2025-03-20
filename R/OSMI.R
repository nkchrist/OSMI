#'Missing values imputation on single or multiple samples from expression data
#'
#' @description OSMI imputes missing DNA methylation values of CpGs on a chromosome within the strand.
#' It replaces missing values by methylation values of nearest CpGs (in base pairs). It can be applied in single-sample applications.
#' @param extendedBetaData An expression data frame containing CpGs in rows and samples and  CpG position in columns.
#' @param positionColumnName A string names corresponding of the column with CpG positions
#' @param subData A data frame (for experimental purposes only); a sub data frame of extendedBetaData. The imputation is only
#'                  carried out for this data frame in order to reduce the computational effort. The nearest neighbour CpG search
#'                  is always based on extendedBetaData
#'
#' @return A list with two components named 'imputed' and 'info'.
#'          'imputed'  a data frame, where missing values have been imputed
#'           'info'    a data frame with some additional information about each imputation
#'           (row, column, imputed value, distance in base pairs to  the nearest neighbour
#'           if detectable, kind of imputation (mean or nearest neighbour))
#' @export
#'
#' @examples
#' OSMI(extendedBetaData, "pos", subData)
#' @author author Lutz Leistritz
#   Institute of Medical Statistics, Computer and Data Sciences
#   Jena University Hospital, Jena, Germany
#
# Version
#   1.0 (initial release)
#
# Date
#   April 17, 2024
OSMI <-function(extendedBetaData, positionColumnName, subData=NULL) {
  extendedBetaData<- extendedBetaData[order(extendedBetaData[,positionColumnName]),]
  # to check the presence of subdata
  if (is.null(subData)) {
    subData <- extendedBetaData
  }

  # to init
  CpG_IDs <- rownames(subData)
  sample_IDs <- colnames(subData)
  replaced_Data <- subData

  # to check data and parameters
  if (!positionColumnName %in% colnames(extendedBetaData)) {stop("CpG positions not found in 'extendedBetaData'.")}
  if (!positionColumnName %in% sample_IDs) {stop("CpG positions not found in 'subData'.")}
  if (!all(sample_IDs %in% colnames(extendedBetaData))) {stop("'subData' is not part of 'extendedBetaData'.")}
  if (!all(CpG_IDs %in% rownames(extendedBetaData))) {stop("'subData' is not part of 'extendedBetaData'.")}

  # to check missing CpG positions in the small data set
  missing_indices <- which(is.na(subData), arr.ind=TRUE)
  position_column <- which(sample_IDs == positionColumnName)
  sample_columns <- which(sample_IDs != positionColumnName)

  # just to allocate memory of info dataframe
  info_Data <- data.frame(matrix(vector(), nrow=nrow(missing_indices), ncol=5))
  colnames(info_Data) <- c("row", "column", "imputed", "distance", "type")
  info_counter <- 1

  # to process rows with missing CpG positions
  missing_rows_mean <- missing_indices[missing_indices[,2] == position_column,1]

  # mean value imputation
  for (row in CpG_IDs[missing_rows_mean]) {
    affected_columns <- sample_columns[is.na(subData[row, sample_columns])]
    subData[row, positionColumnName] <- -9999 # to make CpG position column 'complete'
    for (col in sample_IDs[affected_columns]) {
      replaced_Data[row, col] <- mean(as.numeric(extendedBetaData[complete.cases(extendedBetaData[,col]),col]))
      info_Data[info_counter,] <- c(row, col, replaced_Data[row, col], NA, "mean")
      info_counter <- info_counter + 1
    }
  }

  # array indices of missing values
  missing_indices_NN <- which(is.na(subData) & is.na(replaced_Data), arr.ind=TRUE)

  # loop through samples with missings (column by column (sample by sample))
  for (sample in unique(missing_indices_NN[,2])) {

    # setup lookup table based on all CpG positions
    lookup <- extendedBetaData[,c(sample_IDs[sample], positionColumnName)]
    lookup <- lookup[complete.cases(lookup), ]
    lookup_length <- nrow(lookup)
    if (lookup_length > 0) {
      # replace methylation levels of sample: sample_IDs[sample]
      row_idx <- which(missing_indices_NN[,2] == sample) # row numbers in missing_indices_NN
      row_nmb <- missing_indices_NN[row_idx,1] # row numbers in subData
      CpGs_to_replace <- CpG_IDs[row_nmb]
      query <- subData[CpGs_to_replace, positionColumnName]

      # to find nearest neighbors
      left_index <- findInterval(query, lookup[,2])
      left_index[left_index == 0] <- 1 # just to avoid index out of bounds
      right_index <- left_index + 1
      right_index[right_index > lookup_length] <- lookup_length # just to avoid index out of bounds
      left_position <- lookup[left_index,2]
      right_position <- lookup[right_index,2]
      to_left <- (query - left_position) <= (right_position - query)
      to_right <- !to_left
      nearest_index <- left_index * to_left + right_index * to_right
      nearest_position <- lookup[nearest_index,2]
      replaced_Data[CpGs_to_replace,sample_IDs[sample]] <- lookup[nearest_index, 1]

      # to save imputation details
      n_imputed <- length(CpGs_to_replace)
      index <- info_counter : (info_counter + n_imputed - 1)
      
      info_Data[index,1] <- CpGs_to_replace
      info_Data[index,2] <- rep(sample_IDs[sample], n_imputed)
      info_Data[index,3] <- replaced_Data[CpGs_to_replace,sample_IDs[sample]]
      info_Data[index,4] <- abs(nearest_position - query)
      info_Data[index,5] <- rep("NN", n_imputed)
      info_counter <- info_counter + n_imputed
    }
    else {
      warning(paste("No complete data available for sample ", sample_IDs[sample]))
    }
  }

  # output
  sample_idx <- !sample_IDs %in% positionColumnName
  out <- list(imputed=replaced_Data[,sample_idx], info=info_Data[1:(info_counter-1),])
  return(out)

}

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("minfi")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("IlluminaHumanMethylationEPICv2anno.20a1.hg38")
