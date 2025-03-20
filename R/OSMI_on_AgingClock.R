

#' A function to show the impact of missing values imputation methods on aging clock
#'
#' @param extendedBetaData An expression data frame containing CpGs in rows and samples and  CpG position in columns
#' @param metaData Meta data frame samples informations
#' @param coefHorvath Horvath reference file
#'
#' @return
#' @export
#'
#' @examples
#'Impact_of_imputation_on_horvath_clock(extendedBetaData, metaData, coefHorvath)
#' 
Impact_of_imputation_on_horvath_clock<-function(extendedBetaData, metaData, coefHorvath){
  
  
  df_annotation<-data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19))[,c("chr","strand","Name","pos","Islands_Name")]
  
  extendedBetaData <-extendedBetaData_prep(extendedBetaData,df_annotation)
  if(!is.null(extendedBetaData[[2]])){
    subset<-subset[!rownames(subset)%in%extendedBetaData[[2]],]
  }
  metaData<-metaData[!is.na(metaData$age),]
  subData<-extendedBetaData[coefHorvath$CpGmarker[-1],rownames(metaData)]
  subData[]<-lapply(subData, as.numeric)
  subData.h.knn<-impute.knn(as.matrix(subData))$data
  subData.h.knn<-data.frame(cbind("CpGs_Names"=rownames(subData.h.knn),subData.h.knn))
  subData.h.knn[,-1]<-lapply(subData.h.knn[,-1], as.numeric)
  age.subData.knn<-DNAmAge(subData.h.knn, clocks ="Horvath")
  
  subData.meth[]<-lapply(subData.meth,as.numeric)
  subData.meth<-t(subData.meth)
  subData.meth<-methyLImp(subData.meth, min = 0, max = 1, max.sv = NULL, col.list = NULL)
  subData.meth<-data.frame(t(subData.meth))
  subData.meth<-cbind("CpGs_Names"=rownames(subData.meth),subData.meth)
  subData.meth[,-1]<-lapply(subData.meth[,-1], as.numeric)
  age_subData.meth<-DNAmAge(subData.meth, clocks ="Horvath")
  
  extendedBetaData.osmi <-extendedBetaData_prep(extendedBetaData,df_annotation)
  subData.osmi <-subData_prep(extendedBetaData,df_annotation)
  v_0<-c()
  for(i in 1:length(subData.osmi)){
    if(nrow(subData.osmi[[i]])==0 | sum(is.na(subData.osmi[[i]]))==0){
      v_0<-c(v_0,i)
    }
  }
  subData.osmi[v_0]<-NULL
  v<-names(extendedBetaData.osmi.osmi)[!(names(extendedBetaData.osmi.osmi) %in% names(subData.osmi))]
  extendedBetaData.osmi.osmi[v]<-NULL
  subData.osmi<-subData.osmi[order(names(subData.osmi))]
  extendedBetaData.osmi.osmi<-extendedBetaData.osmi.osmi[order(names(extendedBetaData.osmi.osmi))]
  subData.osmi.imp<-mcmapply(OSMI,extendedBetaData.osmi,"pos",subData.osmi,mc.cores=ncores)
  subData.osmi.impute<-subData.osmi.imp[1,][[1]]
  for(i in 2:ncol(subData.osmi.imp)){
    subData.osmi.impute<-rbind(subData.osmi.impute,subData.osmi.imp[1,][[i]])
  }
  
  subData.osmi<-cbind("CpGs_Names"=rownames(subData.osmi),subData.osmi)
  subData.osmi[,-1]<-lapply(subData.osmi[,-1], as.numeric)
  age_subData.osmi<-DNAmAge(subData.osmi, clocks ="Horvath")
  
  
  return(list(age_subData.osmi,age.subData.knn,age_subData.meth))
}


