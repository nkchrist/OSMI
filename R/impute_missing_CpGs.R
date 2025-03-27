#'sample_generator subsets a large data set into a smaller and introduces missing values
#'
#'@description
#'sample_generator first resamples a large data set into a smaller according to the requested number of rows and columns,
#'then adds 10 % missing values to 90 % of rows.
#'
#' @param extendedBetaData An complete cases data frame.
#' @param rnum Number of rows to be selected
#' @param cnum Nmber of columns to be selected
#'
#' @return A list of data frames containing in 1. the subset without missing values, 2. The  missing values, 3. The missng values possitions.
#' @export
#'
#' @examples
#' x<-sample_generator(extendedBetaData, rnum=1, cnum=50)




sample_generator<-function(extendedBetaData=extendedBetaData,rnum=rnum,cnum=cnum){

  data.rnames<-base::sample(rownames(extendedBetaData),rnum)

  data.cnames<-base::sample(colnames(extendedBetaData),cnum)

  if(missingArg(rnum)){
    extendedBetaData<-extendedBetaData[,data.cnames]
  }else if(missingArg(cnum)){
    extendedBetaData<-extendedBetaData[data.rnames,]
  }else{
    extendedBetaData<-extendedBetaData[data.rnames,data.cnames]
  }
  extendedBetaData<-data.frame(extendedBetaData)
  colnames(extendedBetaData) <-data.cnames
  rownames(extendedBetaData) <-data.rnames

  if(ncol(extendedBetaData)==1){
    data.miss.pos<-sample(1:length(data.rnames),(rnum-rnum/10)/10)
    data.miss<-extendedBetaData
    data.miss[data.miss.pos,1]<-NA
    data.miss<-data.frame(data.miss)
    missing_cells <- data.frame(which(is.na(data.miss), arr.ind = TRUE))
    missing_cells$colnames<-colnames(data.miss)
    missing_cells$rownames<-rownames(missing_cells)
  }else{
    extendedBetaData[]<- lapply(extendedBetaData, as.numeric)
    data.miss<-extendedBetaData
    data.miss[1:(rnum-rnum/10),]<-missMethods::delete_MCAR(data.miss[1:(rnum-rnum/10),],0.1)
    missing_cells <- data.frame(which(is.na(data.miss), arr.ind = TRUE))
    missing_cells$colnames<-colnames(data.miss[,missing_cells[,2]])
    missing_cells$colnames<-gsub("\\..*","",missing_cells$colnames)
    missing_cells$rownames<-rownames(missing_cells)
    missing_cells$rownames<-gsub("\\..*","",missing_cells$rownames)
  }
  return(list(extendedBetaData,data.miss,missing_cells))
}

##
#' Use the function impute.knn to impute missing values
#'
#' @param datalist A list of data frames containing in 1. the subset without missing values, 2. The  missing values, 3. The missng values possitions.
#' @importFrom impue impute.knn
#' @returns a list containing a data frame with records of the impution, and the execution time
#'
#' @examples
#' x.knn<-use_imputeknn(x) ### x is the list resulted from sample_generator#'


use_imputeknn<-function(datalist=datalist){
  missing_cells<-datalist[[3]]
  data.inc<-as.matrix(datalist[[2]])
  data.cc<-data.frame(datalist[[1]])
  tic()
  data.knn<-impute::impute.knn(data.inc)$data
  exectime <- toc()
  data.knn<-data.frame(data.knn)
  true_val<-c()
  knn<-c()
  for(i in 1:nrow(missing_cells)){
    true_val<-c(true_val,data.cc[missing_cells[i,]$rownames,missing_cells[i,]$colnames])
    knn<-c(knn,data.knn[missing_cells[i,]$rownames,missing_cells[i,]$colnames])
  }
  data.knn.imp<-data.frame(cbind("rownames"=missing_cells$rownames,"colnames"=missing_cells$colnames,
                                 "true_val"=true_val,
                                 "knn"=knn))


  return(list(data.knn.imp,exectime))
}

#' Use the function methyLImp to impute missing values
#'
#' @param datalist A list of data frames containing in 1. the subset without missing values, 2. The  missing values, 3. The missng values possitions.
#' @importFrom methyLImp methyLImp
#' @returns A list containing a data frame with records of the imputation, and the execution time
#' @examples
#' use_methylimp(x) ### x is the list resulted from sample_generator

use_methylimp<-function(datalist=datalist){
  missing_cells<-datalist[[3]]
  data.cc<-as.matrix(datalist[[1]])
  data.inc<-t(datalist[[2]]) ### methyLImp expect CpGs in columns and samples in rows
  tic()
  data.meth<-methyLImp::methyLImp(data.inc, min = 0, max = 1)
  exectime <- toc()
  data.meth<-data.frame(t(data.meth))
  true_val<-c()
  meth<-c()
  for(i in 1:nrow(missing_cells)){
    true_val<-c(true_val,data.cc[missing_cells[i,]$rownames,missing_cells[i,]$colnames])
    meth<-c(meth,data.meth[missing_cells[i,]$rownames,missing_cells[i,]$colnames])
  }
  data.meth.imp<-data.frame(cbind("rownames"=missing_cells$rownames,"colnames"=missing_cells$colnames,
                                  "true_val"=true_val,
                                  "meth"=meth))


  return(list(data.meth.imp,exectime))
}

#' extendedBetaData_prep splits an expression  data frame over chromosomes information into a list of data frames
#'
#' @param extendedBetaData An expression data frame containing CpGs in rows and samples and  CpG position in columns.
#' @importFrom minfi getAnnotation
#' @export
#' @returns A list of data frames in which CpGs were grouped depending on their position by chromosomes and strand
#' @examples
#' extendedBetaData_prep(extendedBetaData)


extendedBetaData_prep<-function(extendedBetaData){
  df_annotation<-data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19))[,c("chr","strand","Name","pos","Islands_Name")]
  extendedBetaData<- cbind("chr" =  df_annotation$chr[match(rownames(extendedBetaData), df_annotation$Name)],
                     "pos" = df_annotation$pos[match(rownames(extendedBetaData), df_annotation$Name)],
                     "strand" =  df_annotation$strand[match(rownames(extendedBetaData), df_annotation$Name)],extendedBetaData)

  extendedBetaData <- split(extendedBetaData, f=list(extendedBetaData[,"chr"], extendedBetaData[,"strand"]))

  double_pos<-c()
  for (k in seq(1,length(extendedBetaData))) {
    extendedBetaData[[k]] <- extendedBetaData[[k]][order(extendedBetaData[[k]][,"pos"]),]
    diffs <- diff(extendedBetaData[[k]][,"pos"]) # to delete double CpG positions
    idx <- which(diffs == 0)
    if (length(idx) > 0) {
      double_pos<-rownames(extendedBetaData[[k]][idx,])
      extendedBetaData[[k]] <- extendedBetaData[[k]][-idx,]
    }
  }
  return(list(extendedBetaData,double_pos))
}


#' Use OSMI to impute missing values at a single and multiple subject levels.
#'
#' @param extendedBetaData An expression data frame containing CpGs in rows and samples and  CpG position in columns.
#
#' @param subData A data frame (for experimental purposes only); a sub data frame of extendedBetaData. The imputation is only
#                      carried out for this data frame in order to reduce the computational effort. The nearest neighbour CpG search
#                      is always based on extendedBetaData.
#' @param ncores Number of system cores to be used
#'
#' @return A list containing a data frame with records of the impution, and the execution time
#' @importFrom minfi getAnnotation
#' @export
#'
#' @examples
#' use_osmi(extendedBetaData, subData,df_annotation,ncores)



use_osmi<- function(extendedBetaData, subData, df_annotation, ncores){
  subData.cc<-subData[[1]]
  missing_cells<-subData[[3]]
  subData<-data.frame(subData[[2]])

  subData<- cbind("chr" =  df_annotation$chr[match(rownames(subData), df_annotation$Name)],
                  "pos" = df_annotation$pos[match(rownames(subData), df_annotation$Name)],
                  "strand" =  df_annotation$strand[match(rownames(subData), df_annotation$Name)],subData)


  if(sum(is.na(subData[,"chr"]))!=0 || sum(is.na(subData[,"strand"]))!=0){
    subData[is.na(subData[,"chr"]),][,"chr"]<-"chrNA"
    subData[is.na(subData[,"strand"]),][,"strand"]<-"+"
  }

  subData <- split(subData, f=list(subData[,"chr"], subData[,"strand"]))

  if("chrNA.+"%in%names(subData)){
    subData_NA<-subData[["chrNA.+"]]
    subData[c("chrNA.-","chrNA.+")]<-NULL
  }

  v_0<-c()
  for(i in 1:length(subData)){
    if(nrow(subData[[i]])==0 | sum(is.na(subData[[i]]))==0){
      v_0<-c(v_0,i)
    }
  }

  subData[v_0]<-NULL
  v<-names(extendedBetaData)[!(names(extendedBetaData) %in% names(subData))]
  extendedBetaData[v]<-NULL
  subData<-subData[order(names(subData))]
  extendedBetaData<-extendedBetaData[order(names(extendedBetaData))]
  for (i in names(subData)) {
    missing_indices_NN <- data.frame(which(is.na(subData[[i]]), arr.ind=TRUE))
    missing_indices_NN$colnames<-colnames(subData[[i]][,c(1,missing_indices_NN[,2])])[-1]
    missing_indices_NN$colnames<-gsub("\\..*","",missing_indices_NN$colnames)
    missing_indices_NN$rownames<-rownames(missing_indices_NN)
    missing_indices_NN$rownames<-gsub("\\..*","",missing_indices_NN$rownames)

    for(j in 1:nrow(missing_indices_NN)){
      extendedBetaData[[i]][missing_indices_NN[j,"rownames"],missing_indices_NN[j,"colnames"]]<-NA
    }

  }
  tic()
  subData.imp<-mcmapply(OSMI,extendedBetaData,"pos",subData,mc.cores=ncores)

  subData.impute<-subData.imp[1,][[1]]
  if(exists("subData_NA")){
    subData_NA<-subData_NA[,!colnames(subData_NA)==c("pos")]
    for(i in  1:ncol(subData.imp)){
      subData_NA<-rbind(subData.imp[1,][[i]],subData_NA)
    }
    for (i in colnames(subData_NA)[-1:-2]) {
      subData_NA[,i]<-as.numeric(subData_NA[,i])
      subData_NA[,i][is.na(subData_NA[,i])] <- round(mean(subData_NA[,i], na.rm = TRUE),2)
    }
    subData.impute<-subData_NA
  }else{
    for(i in 2:ncol(subData.imp)){
      subData.impute<-rbind(subData.impute,subData.imp[1,][[i]])
    }

  }
  exectime <- toc()

  true_val<-c()
  os<-c()
  for(i in 1:nrow(missing_cells)){
    true_val<-c(true_val,subData.cc[missing_cells[i,]$rownames,missing_cells[i,]$colnames])
    os<-c(os,subData.impute[missing_cells[i,]$rownames,missing_cells[i,]$colnames])
  }
  data.osmi.imp<-data.frame(cbind("rownames"=missing_cells$rownames,"colnames"=missing_cells$colnames,
                                 "true_val"=true_val,
                                 "osmi"=os))

  return(list(data.osmi.imp,exectime))
}



# subset<-eliminateMissingsMulti(extendedBetaData)

#' Use OSMI_island to impute missing values at a single and multiple subject levels.
#'
#' @param extendedBetaData An expression data frame containing CpGs in rows and samples and  CpG position in columns.
#' @param df_annotation A reference data frame containing CpG sites informations
#' @param ncores Number of system cores to be used
#'
#' @return A data frame with imputation scores
#' @importFrom minfi getAnnotation
#' @export

#' @examples
#' use_osmi_island(blood_download,4)

use_osmi_island<- function(extendedBetaData,ncores){
  
  df_annotation<-data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19))[,c("chr","strand","Name","pos","Islands_Name")]
  extendedBetaData<-extendedBetaData[!rownames(extendedBetaData)%in%c("sample_id","tissue"),]
  subData <-cbind("chr" =  df_annotation$chr[match(rownames(extendedBetaData), df_annotation$Name)],
                  "pos" = df_annotation$pos[match(rownames(extendedBetaData), df_annotation$Name)],
                  "strand" =  df_annotation$strand[match(rownames(extendedBetaData), df_annotation$Name)],
                  "Islands_Name" =  df_annotation$Islands_Name[match(rownames(extendedBetaData), df_annotation$Name)],
                  extendedBetaData)
  extendedBetaData <-extendedBetaData_prep(extendedBetaData)
  if(!is.null(extendedBetaData[[2]])){
    subset<-subset[!rownames(subset)%in%extendedBetaData[[2]],]
  }
  extendedBetaData<-extendedBetaData[[1]]
  print("extended data ready")
  
  subData <- subData[subData[,"Islands_Name"]!="",]
  subData.table<-data.frame(table(subData[,"Islands_Name"]))
  x<-subData.table[subData.table$Freq<=1,]$Var1
  subData.miss<-subData[!subData[,"Islands_Name"]%in%x,]
  missing_cells <- data.frame(which(is.na(subData.miss[,!colnames(subData.miss)%in%c("chr","pos","strand","Islands_Name")])
                                    , arr.ind = TRUE))
  missing_cells$colnames<-colnames(subData.miss[,missing_cells[,2]])
  missing_cells$colnames<-gsub("\\..*","",missing_cells$colnames)
  subData.miss.island <- split(subData.miss, f=list(subData.miss[,"Islands_Name"]))
  print("subdata island ready")
  subData.miss <-extendedBetaData_prep(subData.miss)
  if(!is.null(subData.miss[[2]])){
    subset<-subset[!rownames(subset)%in%subData.miss[[2]],]
  }
  subData.miss<-subData.miss[[1]]
  print("subdata ready")
  
  
  
  v_0<-c()
  for(i in 1:length(subData.miss)){
    if(nrow(subData.miss[[i]])==0 | sum(is.na(subData.miss[[i]]))==0){
      v_0<-c(v_0,i)
    }
  }
  
  subData.miss[v_0]<-NULL
  v<-names(extendedBetaData)[!(names(extendedBetaData) %in% names(subData.miss))]
  extendedBetaData[v]<-NULL
  subData.miss<-subData.miss[order(names(subData.miss))]
  extendedBetaData<-extendedBetaData[order(names(extendedBetaData))]
  
  subData.miss.imp_island<-mcmapply(OSMI,subData.miss.island,"pos",subData.miss.island,mc.cores=ncores)
  subData.miss.imp<-mcmapply(OSMI,extendedBetaData,"pos",subData.miss,mc.cores=ncores)
  
  return(list(subData.miss.imp,subData.miss.imp_island))
}

#subData.imp<-use_osmi_island(extendedBetaData, 4)



#'Impute_missing_CpGs includes OSMI, impute.knn, methyLImp to impute missing values
#'
#' @param extendedBetaData An expression data frame containing CpGs in rows and samples and  CpG position in columns
#'
#' @param nRows A list or Number of rows to be selected
#' @param nCols A list or Number of columns to be selected
#' @param nRepeat A Number of runs
#' @param ncores Number of system cores to be used
#'
#' @return A list of data  frames containing, 1. the subset without missing values,
#' 2. imputation records with impute.knn, 3. imputation records with methyLImp, 4. imputation records with OSMI
#' @importFrom minfi getAnnotation
#' @export
#'
#' @example
#' impute_missing_CpGs(extendedBetaData,nRows=c(400, 1000,3000,5000 ),nCols=nCols<-c(1,2,5,10,50, 100, 500, 1000, 1500),5,4)




impute_missing_CpGs<-function(extendedBetaData,nRows=nRows,nCols=nCols,
                    nRepeat,ncores){

   df_annotation<-data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19))[,c("chr","strand","Name","pos","Islands_Name")]
   extendedBetaData<-extendedBetaData[!rownames(extendedBetaData)%in%c("sample_id","tissue"),]
   subset<-eliminateMissingsMulti(extendedBetaData)
   subset<-extendedBetaData[subset[["rows"]],subset[["cols"]]]

   extendedBetaData <-extendedBetaData_prep(extendedBetaData)
   if(!is.null(extendedBetaData[[2]])){
     subset<-subset[!rownames(subset)%in%extendedBetaData[[2]],]
   }
   
  extendedBetaData<-extendedBetaData[[1]]

  #create results dir
  if(!dir.exists(file.path(paste(getwd(),"results",sep = "/")))){
    dir.create(file.path(paste(getwd(),"results",sep = "/")))
  }
  
  for(n in 1:nRepeat){

    ifelse(!dir.exists(file.path(paste(getwd(),"results",sep = "/"), paste("results",n,sep = "_"))),
           dir.create(file.path(paste(getwd(),"results",sep = "/"), paste("results",n,sep = "_"))), "directory exists")
    
    print(paste("Subfolder result",n,"created", sep=" "))
    
    for(i in nRows){
      for(j in nCols){
        subData<-sample_generator(extendedBetaData=subset,rnum=i,cnum=j)  ##change input data according
        names(subData)<-c("subset","subset_NA","missing_cells")
        print(dim(subData[[2]]))
        if(j==1){
          print("impute.knn and methyLImp not applicable")
          subData.osmi<-use_osmi(extendedBetaData, subData,df_annotation,ncores)

          data.meth.knn.imp<-subData.osmi

        }
        else if(j==2){

          print("methyLImp not applicable, imputation done with OSMI and impute.knn")
          subData.knn<-use_imputeknn(subData)
          subData.osmi<-use_osmi(extendedBetaData, subData,df_annotation,ncores)
          data.meth.knn.imp<-list(subData,subData.knn,subData.osmi)

        }
         else{
         subData.knn<-use_imputeknn(subData)
         subData.meth<-use_methylimp(subData)
         subData.osmi<-use_osmi(extendedBetaData,subData,df_annotation, ncores)
         data.meth.knn.imp<-list(subData,subData.knn,subData.meth,subData.osmi)
        }
        saveRDS(data.meth.knn.imp,file =paste(file.path(paste(getwd(),"results",sep = "/"),
                                                        paste("results",n,sep = "_")),paste(i,j,".rds",sep = "_"),sep = "/"))
        print("results available in results folder")
      }

    }

  }
  return(data.meth.knn.imp)
}




dataFrame.imputed<-function(subData.impute,missing_cells){
  subData.impute<-subData.miss.imp[[1]]
  for(i in 2:ncol(subData.miss.imp)){
    subData.impute<-rbind(subData.impute,subData.miss.imp[[i]])
  }

  imputed<-c()
  for(i in 1:nrow(missing_cells)){
    imputed<-c(imputed,subData.impute[missing_cells[i,]$rownames,missing_cells[i,]$colnames])
  }
  data.imputed<-data.frame(cbind("rownames"=missing_cells$rownames,"colnames"=missing_cells$colnames,
                                  "imputed"=imputed))
  return(data.imputed)
}
###Variables

positionColumnName <- "pos"
islandName<- "Islands_Name"
nRows<-c(400, 1000,3000,5000 )
nCols<-c(1,2,5,10,50, 100, 500, 1000, 1500)
nRepeat<-8

