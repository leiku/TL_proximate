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
n.runs<-100
nN=c(10,50,100,150,200)   # location
rhos<-seq(0, 0.9, 0.1)    

nT=40               # time
filename=paste0("Res_test_extended_dim_",nT)
Fn_simulation(nT, nN, rhos, n.runs, filename)
Res<-readRDS(paste0(filename,".RDS"))
Fn_averagd_igood(Res, file.name1=paste0("res_test_extended_dim_averaged_",nT), 
                 file.name2=paste0("res_test_extended_dim_igood_",nT))


nT=60               # time
filename=paste0("Res_test_extended_dim_",nT)
Fn_simulation(nT, nN, rhos, n.runs, filename)
Res<-readRDS(paste0(filename,".RDS"))
Fn_averagd_igood(Res, file.name1=paste0("res_test_extended_dim_averaged_",nT), 
                 file.name2=paste0("res_test_extended_dim_igood_",nT))

nT=80               # time
filename=paste0("Res_test_extended_dim_",nT)
Fn_simulation(nT, nN, rhos, n.runs, filename)
Res<-readRDS(paste0(filename,".RDS"))
Fn_averagd_igood(Res, file.name1=paste0("res_test_extended_dim_averaged_",nT), 
                 file.name2=paste0("res_test_extended_dim_igood_",nT))




#plot check TL
nT<-40
res<-readRDS(paste0("res_test_extended_dim_averaged_",nT,".RDS"))
i.good.comb<-readRDS(paste0("res_test_extended_dim_igood_",nT,".RDS"))
Fn_plot_check_hold(res, file.name=paste0("Fig_test_extended_dim_",nT))
Fn_plot_barplot(res, i.good.comb, file.name=paste0("Fig_barplot_extended_dim_",nT))

nT<-60
res<-readRDS(paste0("res_test_extended_dim_averaged_",nT,".RDS"))
i.good.comb<-readRDS(paste0("res_test_extended_dim_igood_",nT,".RDS"))
Fn_plot_check_hold(res, file.name=paste0("Fig_test_extended_dim_",nT))
Fn_plot_barplot(res, i.good.comb, file.name=paste0("Fig_barplot_extended_dim_",nT))

nT<-80
res<-readRDS(paste0("res_test_extended_dim_averaged_",nT,".RDS"))
i.good.comb<-readRDS(paste0("res_test_extended_dim_igood_",nT,".RDS"))
Fn_plot_check_hold(res, file.name=paste0("Fig_test_extended_dim_",nT))
Fn_plot_barplot(res, i.good.comb, file.name=paste0("Fig_barplot_extended_dim_",nT))




