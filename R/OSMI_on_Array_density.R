#' Impact of the density of the expression data on OSMI
#'
#' @param blood_download An expression data frame containing CpGs in rows and samples in columns
#
#' @param ncores Number of system cores to be used
#' @description
#' OSMI_on_array_density figures out the impact of the density of an methylation array on OSMI imputation.
#' The denser the array, the should be the performance of OSMI.
#'
#' @return A data frame with 2 columns, the rmse of the imputation per samples
#' and the length of the array used for the imputation
#'
#' @export

#' @examples
#' OSMI_on_array_density(blood_download,4)

function (extendedBetaData, subset, df_annotation, ncores = ncores){
  
  df_annotation<-data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19))[,c("chr","strand","Name","pos","Islands_Name")]
  extendedBetaData<-extendedBetaData[!rownames(extendedBetaData)%in%c("sample_id","tissue"),]
  subset<-eliminateMissingsMulti(extendedBetaData)
  subset<-extendedBetaData[subset[["rows"]],subset[["cols"]]]
  subset.colnames <- sample(colnames(subset), 1)                                            
  extendedBetaData <- cbind(chr = df_annotation$chr[match(rownames(extendedBetaData),df_annotation$Name)], 
                            pos = df_annotation$pos[match(rownames(extendedBetaData),df_annotation$Name)], 
                            strand = df_annotation$strand[match(rownames(extendedBetaData),df_annotation$Name)], extendedBetaData)                                               
   subset <- cbind(chr = df_annotation$chr[match(rownames(subset),df_annotation$Name)], 
                   pos = df_annotation$pos[match(rownames(subset),df_annotation$Name)], 
                   strand = df_annotation$strand[match(rownames(subset),df_annotation$Name)], subset) 
    missing_cells<-sample(rownames(subset),5000)
    subset<-subset[missing_cells,c("chr","pos","strand",subset.colnames)]
    subset.miss<-subset
    subset.miss[missing_cells,subset.colnames]<-NA
    y_ext<-extendedBetaData[,c("chr","pos","strand",subset.colnames)]
    y_ext[missing_cells,subset.colnames]<-NA
    subset.rownames<-rownames(subset)
    subset.miss <- subset.miss[order(subset.miss$chr,subset.miss$pos),]
    subset.miss<- split(subset.miss, f=list(subset.miss[,"chr"], subset.miss[,"strand"]))
    subset.miss<-subset.miss[order(names(subset.miss))]
    set.seed(14184703)
    t<-colnames(extendedBetaData)[!colnames(extendedBetaData)%in%c("chr","pos","strand")]
    ncpgs_rmse<-data.frame(matrix(nrow = length(t),ncol = 2))
    rownames(ncpgs_rmse)<-t
    colnames(ncpgs_rmse)<-c("rmse","nCpGs")
    for(i in t){
      x<-extendedBetaData[,c("chr","pos","strand",i),]
      x.rownames<-rownames(x[complete.cases(x),])
      x.names<-c(x.rownames,missing_cells[!missing_cells%in%x.rownames])
      y<-y_ext[rownames(y_ext)%in%x.names,]
      y <- y[order(y[,"chr"],y[,"strand"],y[,"pos"]),]
      y<- split(y, f=list(y[,"chr"], y[,"strand"]))
      v<-names(y)[!(names(y) %in% names(subset.miss))]
      y[v]<-NULL
      y<-y[order(names(y))]
      subset.miss.imp<-mcmapply(OSMI,y,"pos",subset.miss,mc.cores=ncores)
      subset.miss.imp<-subset.miss.imp[1,]
      subset.imp<-subset.miss.imp[[1]]
      for(j in 2:length(subset.miss.imp)){
        subset.imp<-rbind(subset.imp,subset.miss.imp[[j]])
      }
      subset.imp<-subset.imp[!is.na(subset.imp[,subset.colnames]),]
      x<-rownames(subset.imp)[rownames(subset.imp)%in%rownames(subset)]
      subset<-subset[x,]
      ncpgs_rmse[i,]<-c(rmse(as.numeric(subset.imp[,subset.colnames]),as.numeric(subset[,subset.colnames])),length(x.names))
      print(ncpgs_rmse[i,])
    }
    p<- ggplot(res, aes(nCpGs, rmse))+
      geom_point(size=0.1) +
      ylab("RMSE (ß-values)")+
      xlab("Number of available CpGs (in thousands)")+
      theme(text = element_text(size = 8), axis.text.x = element_text(
        size=8),
        axis.text.y = element_text(
          size=8))
    ggsave("results/blood_download_ncpgs_vs_rmse170.png",
           plot = p,width =170,
           height = 109, units = "mm", dpi = 300)
    
    return(list(ncpgs_rmse,p))
    }


#res<-OSMI_on_array_density(blood_download,ncores=4)

