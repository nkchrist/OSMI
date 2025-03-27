##### Analysis ######

#aggregate performance measures
performance_result <- performance_check(nRows = c(400, 1000, 3000, 5000), 
                           nCols = c(1,2,5,10, 50, 100, 500, 1000, 1500), 
                           nRepeat = 10)


#Create Plots for manuscript
#Plots for Figure 4: Accuracy of missing DNA-methylation value imputation (RMSE)
plots_figure_4 <- performance_plot(performance_result[[1]])

#Plots for Figure 1 Supplement: Accuracy of missing DNA-methylation value imputation (MAE)
plots_figure_1_supplement <- performance_plot(performance_result[[2]])

#Plots for Figure 5: Computational time in seconds of missing DNA-methylation value imputations 
plots_figure_5 <- performance_plot(performance_result[[3]])

#Data for Figure 6
result_osmi_island <- use_osmi_island(blood_download, ncores = 1)

#Plot for Figure 7 
OSMI_on_array_density(blood_download, ncores = 1)
