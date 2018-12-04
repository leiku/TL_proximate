rm(list=ls())

require(miscTools)
require(mvtnorm)
require(lmodel2)
library(moments)
library(nlme)
source('Fn_distributions.R')
source('Fn_tls.R')
source('Fn_simulation.R')
source('Fn_averaged_igood.R')
source('Fn_plot_check_hold.R')
source('Fn_plot_barplot.R')


set.seed(10)
n.runs<-1000
nN=c(20,40,60,80,100)   # location
rhos<-seq(0, 0.9, 0.1)    

nT=60               # time

#get simulation results for all blocks, and save as "Res_simulation.RDS"
filename="Res_simulation"
Fn_simulation(nT, nN, rhos, n.runs, filename)


#get average values for each block, and save as "res_simulation_average.RDS"
#get the index of good blocks (enough replicates which have TL hold), and save as "res_simulation_indexgood"
Res<-readRDS(paste0(filename,".RDS"))
Fn_averagd_igood(Res, file.name1="res_simulation_average", 
                 file.name2="res_simulation_indexgood")


#output the figure of testing TL hold (Fig S4, save as "Fig_test_simulation.tiff") 
#and the barplot of TL slopes (Fig 1, save as "Fig_barplot_simulation.tiff")
res<-readRDS("res_simulation_average.RDS")
i.good.comb<-readRDS("res_simulation_indexgood.RDS")
Fn_plot_check_hold(res, file.name="Fig_test_simulation")
Fn_plot_barplot(res, i.good.comb, nN, file.name="Fig_barplot_simulation")


#plot showing average temporal TL slopes vs. spatial TL slopes
#output: "Fig_simulation_temporal_vs._spatial_average.tif"
res<-readRDS("res_simulation_average.RDS")
i.good.comb<-readRDS("res_simulation_indexgood.RDS")
Fn_simulation_temp_vs_spat(res, i.good.comb)



#get R2 for different regression models, save as "ans_simulation_regression_R2.csv"
#meanwhile save the Results and average results of the good blocks as 
#"Res_simulation_good.RDS" and "res_simulation_average_good.RDS"
Res<-readRDS("Res_simulation.RDS")
res<-readRDS("res_simulation_average.RDS")
i.good.comb<-readRDS("res_simulation_indexgood.RDS")
Fn_simulation_regression_analysis(Res, res, i.good.comb, file.name="ans_simulation_regression_R2")


#final regression model
res.good<-readRDS("res_simulation_average_good.RDS")
z<-lm(ave.slope.spat~ave.J+omega1, data=res.good)
summary(z)

z1<-lm(ave.slope.temp~ave.J+ave.syn.spat, data=res.good)
summary(z1)

