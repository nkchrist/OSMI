


#' Missing values imputation in expression data under consideration of islands structures.
#'
#' @param extendedBetaData An expression data frame containing CpGs in rows and samples and  CpG position in columns.
#' @param positionColumnName A string names corresponding of the column with CpG positions
#' @param islandColumnName A string names corresponding of the column with island names
#' @param subData A data frame (for experimental purposes only); a sub data frame of extendedBetaData. The imputation is only
#'                  carried out for this data frame in order to reduce the computational effort. The nearest neighbour CpG search
#'                  is always based on extendedBetaData
#'
#' @return A data frame with informations about each imputation
#'           (row, column, imputed value, distance in base pairs to  the nearest neighbour
#'           if detectable, kind of imputation (mean or nearest neighbour))
#' @export
#'
#' @examples
#' OSMI_island(blood_download, "pos", "Islands_Name",blood_download_complete)


OSMI_island<- function(extendedBetaData, positionColumnName,islandName, subData=NULL) {
  #extendedBetaData<-extendedBetaData[,!colnames(subData)%in%c("chr","strand")]
  extendedBetaData <- extendedBetaData[order(extendedBetaData[,positionColumnName]),]
 # subData<-subData[,colnames(subData)%in%c("chr","strand" ,"Islands_Name")]

  # to the presence of subdata
  if (is.null(subData)) {
    subData <- extendedBetaData
  }

  # to init
  CpG_IDs <- rownames(subData)
  sample_IDs <- colnames(subData)
  replaced_Data <- subData
  #replaced_Data[] <-lapply(replaced_Data, as.numeric)
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
  info_Data <- data.frame(matrix(vector(), nrow=nrow(missing_indices), ncol=7))
  colnames(info_Data) <- c("row", "column", "imputed", "distance", "type","true","Island")
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
  ###replaced_Data <- subData
  # array indices of missing values
  missing_indices_NN <- data.frame(which(is.na(subData) & is.na(replaced_Data), arr.ind=TRUE))
  # loop through samples with missings (column by column (sample by sample))
  for (sample in unique(missing_indices_NN[,1])) {

    # setup lookup table based on all CpG positions
    x<-which(is.na(subData[CpG_IDs[sample],]))
    missing_cpgs<-colnames(subData[,c(1,x)])
    extendedBetaData[CpG_IDs[sample],missing_cpgs[-1]]<-NA
    lookup <- extendedBetaData[extendedBetaData[,islandName]==extendedBetaData[CpG_IDs[sample],][c(1,islandName)],]
    x<-which(rownames(lookup)==CpG_IDs[sample])
    lookup <-lookup[-x,]
    lookup_length <- ncol(lookup)
    # replace methylation levels of sample: sample_IDs[sample]
    row_idx <- which(missing_indices_NN[,1] == sample) # row numbers in missing_indices_NN
    col_nmb <- missing_indices_NN[row_idx,islandName] # col numbers in subData
    samples_to_replace <- sample_IDs[col_nmb]
    query <- as.numeric(extendedBetaData[CpG_IDs[sample],][1,1])

    # to find nearest neighbors
    while (sum(is.na(replaced_Data[CpG_IDs[sample],]))>=2 && nrow(lookup)>1) {

      left_index <- findInterval(query, lookup[,1])
      left_index[left_index == 0] <- 1 # just to avoid index out of bounds
      right_index <- left_index + 1
      right_index[right_index > lookup_length] <- lookup_length # just to avoid index out of bounds
      left_position <- as.numeric(lookup[left_index,1])
      right_position <- as.numeric(lookup[right_index,1])
      to_left <- (query - left_position) <= (right_position - query)
      to_right <- !to_left
      nearest_index <- left_index * to_left + right_index * to_right
      nearest_position <- lookup[nearest_index,1]
      replaced_Data[CpG_IDs[sample],is.na(replaced_Data[CpG_IDs[sample],])]<-
        lookup[nearest_index, colnames(replaced_Data[CpG_IDs[sample],is.na(replaced_Data[CpG_IDs[sample],])])]
      lookup<-lookup[-nearest_index,]
      #samples_to_replace<-which(is.na(replaced_Data[CpG_IDs[sample],samples_to_replace]))
      #samples_to_replace<-colnames(replaced_Data[CpG_IDs[sample],samples_to_replace])
      print(CpG_IDs[sample])
    }

    while(sum(is.na(replaced_Data[CpG_IDs[sample],]))==1 && nrow(lookup)>1) {

      x<-which(is.na(replaced_Data[CpG_IDs[sample],]))
      missing_cpgs<-colnames(replaced_Data[,c(1,x)])


      left_index <- findInterval(query, lookup[,1])
      left_index[left_index == 0] <- 1 # just to avoid index out of bounds
      right_index <- left_index + 1
      right_index[right_index > lookup_length] <- lookup_length # just to avoid index out of bounds
      left_position <- as.numeric(lookup[left_index,1])
      right_position <- as.numeric(lookup[right_index,1])
      to_left <- (query - left_position) <= (right_position - query)
      to_right <- !to_left
      nearest_index <- left_index * to_left + right_index * to_right
      nearest_position <- lookup[nearest_index,1]
      replaced_Data[CpG_IDs[sample], missing_cpgs[-1]]<-
        lookup[nearest_index, missing_cpgs[-1]]
      lookup<-lookup[-nearest_index,]
      #samples_to_replace<-which(is.na(replaced_Data[CpG_IDs[sample],samples_to_replace]))
      #samples_to_replace<-colnames(replaced_Data[CpG_IDs[sample],samples_to_replace])
      print(CpG_IDs[sample])
    }



    if(sum(is.na(replaced_Data[CpG_IDs[sample],]))!=0){

      replaced_Data[CpG_IDs[sample],is.na(replaced_Data[CpG_IDs[sample],])] <-
        mean(as.numeric(replaced_Data[CpG_IDs[sample],-1]),na.rm=TRUE)

    }

    # to save imputation details
    n_imputed <- length(samples_to_replace)
    index <- info_counter : (info_counter + n_imputed - 1)
    info_Data[index,1] <- rep(CpG_IDs[sample], n_imputed)
    info_Data[index,2] <- samples_to_replace
    info_Data[index,3] <- as.numeric(replaced_Data[CpG_IDs[sample],samples_to_replace])
    info_Data[index,4] <- abs(as.numeric(nearest_position) - query)
    info_Data[index,5] <- rep("NN", n_imputed)

    #info_Data[index,6] <- as.numeric(blood_download_complete[CpG_IDs[sample],samples_to_replace])
    info_Data[index,7] <- lookup[nearest_index, islandName]
    info_counter <- info_counter + n_imputed
    print(nearest_index)

  }

  # output
  sample_idx <- !sample_IDs %in% positionColumnName
  out <- list(imputed=replaced_Data[,sample_idx], info=info_Data[1:(info_counter-1),])

  return(out)

}
