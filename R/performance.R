
### An example of test to estimate performance of the OSMI using the rmse and mae

#' Evaluate the performance of OSMI
#'
#' @param nRows A list or Number of rows to be selected
#' @param nCols A list or Number of rows to be selected
#' @param nRepeat A Number of runs
#'
#' @return A table or a list of tables containing the means values from the repetitions
#' @export perfomance_check
#'
#' @examples
#'perfomance_check(nRows=c(400, 1000,3000,5000 ),nCols=c(1,2,5,10,50, 100, 500, 1000, 1500),nRepeat=10)

perfomance_check<-function(nRows=nRows,nCols=nCols,nRepeat=nRepeat){

  tab<-array(,c(length(nRows),length(nCols),nRepeat))
  rownames(tab)<-nRows
 colnames(tab)<-nCols
  dlist<- vector(mode='list', length=9)
  names(dlist)<-c("rmse_table.osmi" ,"rmse_table.knn" , "rmse_table.meth" ,"mae_table.osmi" , "mae_table.knn" ,  "mae_table.meth",
                  "table.osmi.t" ,"table.knn.t", "table.meth.t")
  for(i in names(dlist)){
    dlist[[i]]<-tab
  }


  for(n in 1:nRepeat){

    paths<-ifelse(dir.exists(file.path(paste(getwd(),"results",sep = "/"), paste("results",n,sep = "_"))),
                  file.path(paste(getwd(),"results",sep = "/"), paste("results",n,sep = "_")), "directory does not exists")

    x<-list.files(path=paths, pattern=NULL, all.files=FALSE,
                  full.names=FALSE)
    print(paste("run",n,"started", sep=" "))
    for(i in x){
      sr<-gsub("_[^_]+$", "", i)
      sc<-strsplit(sr,split = "_")[[1]][2]
      sr<-strsplit(sr,split = "_")[[1]][1]
      if(sc=="1"){
        y<-readRDS(paste(paths,i,sep = "/"))[[1]]
        y_t<-readRDS(paste(paths,i,sep = "/"))[[2]][[1]]
        dlist[["rmse_table.osmi"]][sr,sc,n]<-round(rmse(as.numeric(y[,"true_val"]),as.numeric(y[,"osmi"])),3)
        dlist[["mae_table.osmi"]][sr,sc,n]<-round(mae(as.numeric(y[,"true_val"]),as.numeric(y[,"osmi"])),3)
        dlist[["table.osmi.t"]][sr,sc,n]<-as.numeric(y_t)
      }else if(sc=="2"){
        y<-readRDS(paste(paths,i,sep = "/"))[[3]][[1]]
        z<-readRDS(paste(paths,i,sep = "/"))[[2]][[1]]
        #y_t<-strsplit(readRDS(paste(paths,i,sep = "/"))[[3]][[2]]$callback_msg,split = " ")[[1]][1]
        y_t<-readRDS(paste(paths,i,sep = "/"))[[3]][[2]][[1]]
        z_t<-readRDS(paste(paths,i,sep = "/"))[[2]][[2]][[1]]
        dlist[["rmse_table.osmi"]][sr,sc,n]<-round(rmse(as.numeric(y[,"true_val"]),as.numeric(y[,"osmi"])),3)
        dlist[["mae_table.osmi"]][sr,sc,n]<-round(mae(as.numeric(y[,"true_val"]),as.numeric(y[,"osmi"])),3)
        dlist[["table.osmi.t"]][sr,sc,n]<-as.numeric(y_t)
        dlist[["rmse_table.knn"]][sr,sc,n]<-round(rmse(as.numeric(z[,"true_val"]),as.numeric(z[,"knn"])),3)
        dlist[["mae_table.knn"]][sr,sc,n]<-round(mae(as.numeric(z[,"true_val"]),as.numeric(z[,"knn"])),3)
        dlist[["table.knn.t"]][sr,sc,n]<-as.numeric(z_t)
      }else {
        y<-readRDS(paste(paths,i,sep = "/"))[[4]][[1]]
        z<-readRDS(paste(paths,i,sep = "/"))[[2]][[1]]
        p<-readRDS(paste(paths,i,sep = "/"))[[3]][[1]]
        y_t<-readRDS(paste(paths,i,sep = "/"))[[4]][[2]][[1]]
        z_t<-readRDS(paste(paths,i,sep = "/"))[[2]][[2]][[1]]
        p_t<-readRDS(paste(paths,i,sep = "/"))[[3]][[2]][[1]]
        dlist[["rmse_table.osmi"]][sr,sc,n]<-round(rmse(as.numeric(y[,"true_val"]),as.numeric(y[,"osmi"])),3)
        dlist[["mae_table.osmi"]][sr,sc,n]<-round(mae(as.numeric(y[,"true_val"]),as.numeric(y[,"osmi"])),3)
        dlist[["table.osmi.t"]][sr,sc,n]<-as.numeric(y_t)
        dlist[["rmse_table.knn"]][sr,sc,n]<-round(rmse(as.numeric(z[,"true_val"]),as.numeric(z[,"knn"])),3)
        dlist[["mae_table.knn"]][sr,sc,n]<-round(mae(as.numeric(z[,"true_val"]),as.numeric(z[,"knn"])),3)
        dlist[["table.knn.t"]][sr,sc,n]<-as.numeric(z_t)
        dlist[["rmse_table.meth"]][sr,sc,n]<-round(rmse(as.numeric(p[,"true_val"]),as.numeric(p[,"meth"])),3)
        dlist[["mae_table.meth"]][sr,sc,n]<-round(mae(as.numeric(p[,"true_val"]),as.numeric(p[,"meth"])),3)
        dlist[["table.meth.t"]][sr,sc,n]<-as.numeric(p_t)

      }

    }
  }

  ###Means and confidence intervals

  time_tab_mean<-data.frame(matrix(nrow=length(nRows)*3*3,ncol=length(nCols)))
  colnames(time_tab_mean)<-nCols
  rmse_tab_mean<-data.frame(matrix(nrow=length(nRows)*3*3,ncol=length(nCols)))
  colnames(rmse_tab_mean)<-nCols
  mae_tab_mean<-data.frame(matrix(nrow=length(nRows)*3*3,ncol=length(nCols)))
  colnames(mae_tab_mean)<-nCols
  n<-c()
  Print("Statistic started")
  for(i in nRows){
      n<-c(n,paste(i,"osmi_upper",sep="_"),paste(i,"osmi_mean",sep="_"),paste(i,"osmi_lower",sep="_"),
           paste(i,"knn_upper",sep="_"),paste(i,"knn_mean",sep="_"),paste(i,"knn_lower",sep="_"),
           paste(i,"meth_upper",sep="_"),paste(i,"meth_mean",sep="_"),paste(i,"meth_lower",sep="_")

           )
  }
  rownames(time_tab_mean)<-n
  rownames(rmse_tab_mean)<-n
  rownames(mae_tab_mean)<-n

  for (i in nRows) {
    for (j in nCols) {
      i<-as.character(i)
      j<-as.character(j)
      rmse_tab_mean[paste(i,"osmi_upper",sep="_"),j]<-mean(dlist[["rmse_table.osmi"]][i,j,1:nRepeat],na.rm = TRUE )+
        std_err(dlist[["rmse_table.osmi"]][i,j,1:nRepeat])
      rmse_tab_mean[paste(i,"osmi_mean",sep="_"),j]<-mean(dlist[["rmse_table.osmi"]][i,j,1:nRepeat],na.rm = TRUE )
      rmse_tab_mean[paste(i,"osmi_lower",sep="_"),j]<-mean(dlist[["rmse_table.osmi"]][i,j,1:nRepeat],na.rm = TRUE )-
        std_err(dlist[["rmse_table.osmi"]][i,j,1:nRepeat])
    rmse_tab_mean[paste(i,"knn_upper",sep="_"),j]<-mean(dlist[["rmse_table.knn"]][i,j,1:nRepeat],na.rm = TRUE )+
        std_err(dlist[["rmse_table.knn"]][i,j,1:nRepeat])
      rmse_tab_mean[paste(i,"knn_mean",sep="_"),j]<-mean(dlist[["rmse_table.knn"]][i,j,1:nRepeat],na.rm = TRUE )
      rmse_tab_mean[paste(i,"knn_lower",sep="_"),j]<-mean(dlist[["rmse_table.knn"]][i,j,1:nRepeat],na.rm = TRUE )-
        std_err(dlist[["rmse_table.knn"]][i,j,1:nRepeat])
     rmse_tab_mean[paste(i,"meth_upper",sep="_"),j]<-mean(dlist[["rmse_table.meth"]][i,j,1:nRepeat],na.rm = TRUE )+
        std_err(dlist[["rmse_table.meth"]][i,j,1:nRepeat])
      rmse_tab_mean[paste(i,"meth_mean",sep="_"),j]<-mean(dlist[["rmse_table.meth"]][i,j,1:nRepeat],na.rm = TRUE )
      rmse_tab_mean[paste(i,"meth_lower",sep="_"),j]<-mean(dlist[["rmse_table.meth"]][i,j,1:nRepeat],na.rm = TRUE )-
        std_err(dlist[["rmse_table.meth"]][i,j,1:nRepeat])


      mae_tab_mean[paste(i,"osmi_upper",sep="_"),j]<-mean(dlist[["mae_table.osmi"]][i,j,1:nRepeat],na.rm = TRUE )+
        std_err(dlist[["mae_table.osmi"]][i,j,1:nRepeat])
      mae_tab_mean[paste(i,"osmi_mean",sep="_"),j]<-mean(dlist[["mae_table.osmi"]][i,j,1:nRepeat],na.rm = TRUE )
      mae_tab_mean[paste(i,"osmi_lower",sep="_"),j]<-mean(dlist[["mae_table.osmi"]][i,j,1:nRepeat],na.rm = TRUE )-
        std_err(dlist[["mae_table.osmi"]][i,j,1:nRepeat])
      mae_tab_mean[paste(i,"knn_upper",sep="_"),j]<-mean(dlist[["mae_table.knn"]][i,j,1:nRepeat],na.rm = TRUE )+
        std_err(dlist[["mae_table.knn"]][i,j,1:nRepeat])
      mae_tab_mean[paste(i,"knn_mean",sep="_"),j]<-mean(dlist[["mae_table.knn"]][i,j,1:nRepeat],na.rm = TRUE )
      mae_tab_mean[paste(i,"knn_lower",sep="_"),j]<-mean(dlist[["mae_table.knn"]][i,j,1:nRepeat],na.rm = TRUE )-
        std_err(dlist[["mae_table.knn"]][i,j,1:nRepeat])
      mae_tab_mean[paste(i,"meth_upper",sep="_"),j]<-mean(dlist[["mae_table.meth"]][i,j,1:nRepeat],na.rm = TRUE )+
        std_err(dlist[["mae_table.meth"]][i,j,1:nRepeat])
      mae_tab_mean[paste(i,"meth_mean",sep="_"),j]<-mean(dlist[["mae_table.meth"]][i,j,1:nRepeat],na.rm = TRUE )
      mae_tab_mean[paste(i,"meth_lower",sep="_"),j]<-mean(dlist[["mae_table.meth"]][i,j,1:nRepeat],na.rm = TRUE )-
        std_err(dlist[["mae_table.meth"]][i,j,1:nRepeat])

      time_tab_mean[paste(i,"osmi_upper",sep="_"),j]<-mean(dlist[["table.osmi.t"]][i,j,1:nRepeat],na.rm = TRUE )+
        std_err(dlist[["table.osmi.t"]][i,j,1:nRepeat])
      time_tab_mean[paste(i,"osmi_mean",sep="_"),j]<-mean(dlist[["table.osmi.t"]][i,j,1:nRepeat],na.rm = TRUE )
      time_tab_mean[paste(i,"osmi_lower",sep="_"),j]<-mean(dlist[["table.osmi.t"]][i,j,1:nRepeat],na.rm = TRUE )-
        std_err(dlist[["table.osmi.t"]][i,j,1:nRepeat])
      time_tab_mean[paste(i,"knn_upper",sep="_"),j]<-mean(dlist[["table.knn.t"]][i,j,1:nRepeat],na.rm = TRUE )+
        std_err(dlist[["table.knn.t"]][i,j,1:nRepeat])
      time_tab_mean[paste(i,"knn_mean",sep="_"),j]<-mean(dlist[["table.knn.t"]][i,j,1:nRepeat],na.rm = TRUE )
      time_tab_mean[paste(i,"knn_lower",sep="_"),j]<-mean(dlist[["table.knn.t"]][i,j,1:nRepeat],na.rm = TRUE )-
        std_err(dlist[["table.knn.t"]][i,j,1:nRepeat])
      time_tab_mean[paste(i,"meth_upper",sep="_"),j]<-mean(dlist[["table.meth.t"]][i,j,1:nRepeat],na.rm = TRUE )+
        std_err(dlist[["table.meth.t"]][i,j,1:nRepeat])
      time_tab_mean[paste(i,"meth_mean",sep="_"),j]<-mean(dlist[["table.meth.t"]][i,j,1:nRepeat],na.rm = TRUE )
      time_tab_mean[paste(i,"meth_lower",sep="_"),j]<-mean(dlist[["table.meth.t"]][i,j,1:nRepeat],na.rm = TRUE )-
        std_err(dlist[["table.meth.t"]][i,j,1:nRepeat])


    }
  }
  Print("Done!")
  return(list(rmse_tab_mean,mae_tab_mean,time_tab_mean))
}

#tab_mean<-perfomance_check(nRows=nRows,nCols=nCols,nRepeat=nRepeat)


#' Plot function to present the result in a graphic
#'
#' @param tab_mean A data table with samples id on the columns and fucntion name on rows
#'
#' @return A list of plots
#' @export
#'
#' @examples
#' perfomance_plot<-function(tab_mean=tab_mean)
#'

perfomance_plot<-function(tab_mean=tab_mean){

tab_mean$ID <- factor(rownames(tab_mean))

tab_mean.melt_400<-cbind(melt(tab_mean[c("400_osmi_mean","400_knn_mean" ,"400_meth_mean"),], id='ID'),
                         "CI_upper"=melt(tab_mean[c("400_osmi_upper","400_knn_upper" ,"400_meth_upper"),], id='ID')[,3],
                         "CI_lower"=melt(tab_mean[c("400_osmi_lower","400_knn_lower" ,"400_meth_lower"),], id='ID')[,3])

tab_mean.melt_1000<-cbind(melt(tab_mean[c("1000_osmi_mean","1000_knn_mean" ,"1000_meth_mean"),], id='ID'),
                         "CI_upper"=melt(tab_mean[c("1000_osmi_upper","1000_knn_upper" ,"1000_meth_upper"),], id='ID')[,3],
                         "CI_lower"=melt(tab_mean[c("1000_osmi_lower","1000_knn_lower" ,"1000_meth_lower"),], id='ID')[,3])

tab_mean.melt_3000<-cbind(melt(tab_mean[c("3000_osmi_mean","3000_knn_mean" ,"3000_meth_mean"),], id='ID'),
                          "CI_upper"=melt(tab_mean[c("3000_osmi_upper","3000_knn_upper" ,"3000_meth_upper"),], id='ID')[,3],
                          "CI_lower"=melt(tab_mean[c("3000_osmi_lower","3000_knn_lower" ,"3000_meth_lower"),], id='ID')[,3])
tab_mean.melt_5000<-cbind(melt(tab_mean[c("5000_osmi_mean","5000_knn_mean" ,"5000_meth_mean"),], id='ID'),
                          "CI_upper"=melt(tab_mean[c("5000_osmi_upper","5000_knn_upper" ,"5000_meth_upper"),], id='ID')[,3],
                          "CI_lower"=melt(tab_mean[c("5000_osmi_lower","5000_knn_lower" ,"5000_meth_lower"),], id='ID')[,3])



p1<- ggplot(tab_mean.melt_400, aes(variable, value,shape =as.factor(ID),color=as.factor(ID),group = as.factor(ID)))+
  geom_point() +
  scale_shape_manual(values=c(15, 16, 17), name="",
                     labels=c("Impute.knn", "MethylImp", "OSMI"))+
  scale_color_manual(values=c("#ff0000", "#ff7514", "blue"), name="",
                     labels=c("Impute.knn", "MethylImp", "OSMI"))+
  geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper,colour = as.factor(ID)),width=.1) +
  geom_line(aes(colour = as.factor(ID))) +
  xlab("")+
  ylab("RMSE (ß-value)")+
  ylim(0,0.4)+
  ggtitle("400 CpGs")+
  theme(legend.position = "bottom",legend.text=element_text(size=10),text=element_text(size=8))

p2<- ggplot(tab_mean.melt_1000, aes(variable, value,shape =as.factor(ID),color=as.factor(ID),group = as.factor(ID)))+
  geom_point() +
  scale_shape_manual(values=c(15, 16, 17), name="",
                     labels=c("Impute.knn", "MethylImp", "OSMI"))+
  scale_color_manual(values=c("#ff0000", "#ff7514", "blue"), name="",
                     labels=c("Impute.knn", "MethylImp", "OSMI"))+
  geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper,colour = as.factor(ID)),width=.1) +
  geom_line(aes(colour = as.factor(ID))) +
  xlab("")+
  ylab("")+
  ylim(0,0.4)+
  ggtitle("1000 CpGs")+
  theme(legend.position = "bottom",legend.text=element_text(size=10),text=element_text(size=8))



p3<- ggplot(tab_mean.melt_3000, aes(variable, value,shape =as.factor(ID),color=as.factor(ID),group = as.factor(ID)))+
  geom_point() +
  scale_shape_manual(values=c(15, 16, 17), name="",
                     labels=c("Impute.knn", "MethylImp", "OSMI"))+
  scale_color_manual(values=c("#ff0000", "#ff7514", "blue"), name="",
                     labels=c("Impute.knn", "MethylImp", "OSMI"))+
  geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper,colour = as.factor(ID)),width=.1) +
  geom_line(aes(colour = as.factor(ID))) +
  xlab("Number of samples")+
  ylab("RMSE (ß-value)")+
  ylim(0,0.4)+
  ggtitle("3000 CpGs")+
  theme(legend.position = "bottom",legend.text=element_text(size=10),text=element_text(size=8))
p4<- ggplot(tab_mean.melt_5000, aes(variable, value,shape =as.factor(ID),color=as.factor(ID),group = as.factor(ID)))+
  geom_point() +
  scale_shape_manual(values=c(15, 16, 17), name="",
                     labels=c("Impute.knn", "MethylImp", "OSMI"))+
  scale_color_manual(values=c("#ff0000", "#ff7514", "blue"), name="",
                     labels=c("Impute.knn", "MethylImp", "OSMI"))+
  geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper,colour = as.factor(ID)),width=.1) +
  geom_line(aes(colour = as.factor(ID))) +
  xlab("Number of samples")+
  ylab("")+
  ylim(0,0.4)+
  ggtitle("5000 CpGs")+
  theme(legend.position = "bottom",legend.text=element_text(size=10),text=element_text(size=8))

png(filename = "rmse_blood_download_complete_400_5000_2_1500.10timesrep.point.knn.meth.png", width = 170, height = 130,units = "mm", res=300)
cowplot::plot_grid(p1,p2, p3,p4,ncol =2,nrow=2,label_colour ="Blue3", label_x = '0', label_y = '1')
dev.off()
return(list(p1,p2,p3,p4))
}

#p_mae<-perfomance_plot(tab_mean=mae_tab_mean)
#p_rmse<-perfomance_plot(tab_mean=rmse_tab_mean)



#' A plot for the test on methyLImp 400 CpGs and 1000 Cpgs by 100 repeat
#'
#' @param tab_mean A table with the means results from the imputation
#'
#' @return A list of plots
#' @export
#'
#' @examples
perfomance_plot.t<-function(tab_mean=tab_mean){

  tab_mean$ID <- factor(rownames(tab_mean))

  tab_mean.melt_400.t<-cbind(melt(tab_mean[c("400_osmi_mean.t","400_knn_mean.t" ,"400_meth_mean.t"),], id='ID'),
                           "CI_upper.t"=melt(tab_mean[c("400_osmi_upper.t","400_knn_upper.t" ,"400_meth_upper.t"),], id='ID')[,3],
                           "CI_lower.t"=melt(tab_mean[c("400_osmi_lower.t","400_knn_lower.t" ,"400_meth_lower.t"),], id='ID')[,3])

  tab_mean.melt_1000.t<-cbind(melt(tab_mean[c("1000_osmi_mean.t","1000_knn_mean.t" ,"1000_meth_mean.t"),], id='ID'),
                            "CI_upper.t"=melt(tab_mean[c("1000_osmi_upper.t","1000_knn_upper.t" ,"1000_meth_upper.t"),], id='ID')[,3],
                            "CI_lower.t"=melt(tab_mean[c("1000_osmi_lower.t","1000_knn_lower.t" ,"1000_meth_lower.t"),], id='ID')[,3])

  tab_mean.melt_3000.t<-cbind(melt(tab_mean[c("3000_osmi_mean.t","3000_knn_mean.t" ,"3000_meth_mean.t"),], id='ID'),
                            "CI_upper.t"=melt(tab_mean[c("3000_osmi_upper.t","3000_knn_upper.t" ,"3000_meth_upper.t"),], id='ID')[,3],
                            "CI_lower.t"=melt(tab_mean[c("3000_osmi_lower.t","3000_knn_lower.t" ,"3000_meth_lower.t"),], id='ID')[,3])
  tab_mean.melt_5000.t<-cbind(melt(tab_mean[c("5000_osmi_mean.t","5000_knn_mean.t" ,"5000_meth_mean.t"),], id='ID'),
                            "CI_upper.t"=melt(tab_mean[c("5000_osmi_upper.t","5000_knn_upper.t" ,"5000_meth_upper.t"),], id='ID')[,3],
                            "CI_lower.t"=melt(tab_mean[c("5000_osmi_lower.t","5000_knn_lower.t" ,"5000_meth_lower.t"),], id='ID')[,3])


  p1<- ggplot(tab_mean.melt_400.t, aes(variable, value,shape =as.factor(ID),color=as.factor(ID),group = as.factor(ID)))+
    geom_point() +
    scale_shape_manual(values=c(15, 16, 17), name="",
                       labels=c("Impute.knn", "MethylImp", "OSMI"))+
    scale_color_manual(values=c("#ff0000", "#ff7514", "blue"), name="",
                       labels=c("Impute.knn", "MethylImp", "OSMI"))+
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)),
                  limits=c(10^-4, 10^4)) +
    geom_errorbar(aes(ymin=CI_lower.t, ymax=CI_upper.t,colour = as.factor(ID)),width=.1) +
    geom_line(aes(colour = as.factor(ID))) +
    xlab("Number of samples")+
    ylab("time (s)")+
    ggtitle("400 CpGs")+
    theme(legend.position = "bottom",legend.text=element_text(size=10),text=element_text(size=8))

  p2<- ggplot(tab_mean.melt_1000.t, aes(variable, value,shape =as.factor(ID),color=as.factor(ID),group = as.factor(ID)))+
    geom_point() +
    scale_shape_manual(values=c(15, 16, 17), name="",
                       labels=c("Impute.knn", "MethylImp", "OSMI"))+
    scale_color_manual(values=c("#ff0000", "#ff7514", "blue"), name="",
                       labels=c("Impute.knn", "MethylImp", "OSMI"))+
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)),
                  limits=c(10^-4, 10^4)) +
    geom_errorbar(aes(ymin=CI_lower.t, ymax=CI_upper.t,colour = as.factor(ID)),width=.1) +
    geom_line(aes(colour = as.factor(ID))) +
    xlab("Number of samples")+
    ylab("time (s)")+
    ggtitle("1000 CpGs")+
    theme(legend.position = "bottom",legend.text=element_text(size=10),text=element_text(size=8))



  p3<- ggplot(tab_mean.melt_3000.t, aes(variable, value,shape =as.factor(ID),color=as.factor(ID),group = as.factor(ID)))+
    geom_point() +
    scale_shape_manual(values=c(15, 16, 17), name="",
                       labels=c("Impute.knn", "MethylImp", "OSMI"))+
    scale_color_manual(values=c("#ff0000", "#ff7514", "blue"), name="",
                       labels=c("Impute.knn", "MethylImp", "OSMI"))+
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)),
                  limits=c(10^-4, 10^4)) +
    geom_errorbar(aes(ymin=CI_lower.t, ymax=CI_upper.t,colour = as.factor(ID)),width=.1) +
    geom_line(aes(colour = as.factor(ID))) +
    xlab("Number of samples")+
    ylab("time (s)")+
    ggtitle("3000 CpGs")+
    theme(legend.position = "bottom",legend.text=element_text(size=10),text=element_text(size=8))

  p4<- ggplot(tab_mean.melt_5000.t, aes(variable, value,shape =as.factor(ID),color=as.factor(ID),group = as.factor(ID)))+
    geom_point() +
    scale_shape_manual(values=c(15, 16, 17), name="",
                       labels=c("Impute.knn", "MethylImp", "OSMI"))+
    scale_color_manual(values=c("#ff0000", "#ff7514", "blue"), name="",
                       labels=c("Impute.knn", "MethylImp", "OSMI"))+
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)),
                  limits=c(10^-4, 10^4)) +
    geom_errorbar(aes(ymin=CI_lower.t, ymax=CI_upper.t,colour = as.factor(ID)),width=.1) +
    geom_line(aes(colour = as.factor(ID)) ) + ##,linewidth=0.25
    xlab("Number of samples")+
    ylab("time (s)")+
    ggtitle("5000 CpGs")+
    theme(legend.position = "bottom",legend.text=element_text(size=10),text=element_text(size=8))


  png(filename = "blood_download_complete_400_5000_2_1500.5timesrep.point.knn.meth_timelog10.png",
      width = 170, height = 130,units = "mm", res=300)
  cowplot::plot_grid(p1,p2, p3,p4,ncol =2,nrow=2,label_colour ="Blue3", label_x = '0', label_y = '1')
  dev.off()

  return(list(p1,p2,p3,p4))
}


#p_rmse<-perfomance_plot.t(tab_mean=time_tab_mean)


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19",force = TRUE)

