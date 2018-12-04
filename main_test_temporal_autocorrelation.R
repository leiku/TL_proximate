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
source('Fn_simulation_temp_vs_spat.R')


set.seed(10)
n.runs<-1000
nN=c(20,40,60,80,100)   # location
rhos<-seq(0, 0.9, 0.1)    
rhos.r=0.8
nT=60               # time

#get simulation results for all blocks, and save as "Res_simulation.RDS"
n.all<-n.runs*5*length(nN)*length(rhos) # all entries

dim.ratio<-nT/nN    
name.model<-c("Poisson","Negative binomial","Exponential","Gamma","Lognormal")

models<-list(model_NonIDPois,model_NonIDNegBinom,model_NonIDGamma,model_NonIDGamma,model_NonIDLognorm)
#notice that the third one is actually Exponential distribution
params<-list(c(1,1),c(5,0.4,1),c(1,1,1),c(4,1,1),c(1,1,1))


getE<-function(T, N, rho.c, rho.r){
  mu<-numeric(N)
  Sigma<-matrix(rho.c*(1-rho.r^2),N,N)+diag((1-rho.c)*(1-rho.r^2),N,N)
  E0<-rmvnorm(T-1, mu, Sigma, method='chol')
  
  E<-matrix(NA, T, N)
  Sigma0<-matrix(rho.c, N, N) + diag(1-rho.c, N, N)
  E[1,]<-rmvnorm(1, mu, Sigma0, method='chol')
  for(t in 1:(T-1)){
    E[t+1, ]<- rho.r*E[t,]+E0[t,]
  }
  return(E)
}

Res<-data.frame(rho=rep(NA,n.all), dim.ratio=rep(NA,n.all), distribution=rep(NA,n.all), replicate=rep(NA,n.all),
                slope.spat=rep(NA,n.all), intercept.spat=rep(NA,n.all), p.line.spat=rep(NA,n.all),
                p.hete.spat=rep(NA,n.all), RMSE.spat=rep(NA,n.all),  r.spat=rep(NA,n.all),
                slope.temp=rep(NA,n.all), intercept.temp=rep(NA,n.all), p.line.temp=rep(NA,n.all),
                p.hete.temp=rep(NA,n.all), RMSE.temp=rep(NA,n.all),  r.temp=rep(NA,n.all),
                syn.spat=rep(NA,n.all),syn.temp=rep(NA,n.all),J=rep(NA,n.all),J.all=rep(NA,n.all),
                p.fit.spat=rep(NA, n.all), p.fit.temp=rep(NA, n.all))
a<-1
for(ii in 1:length(name.model)){
  for(jj in 1:length(nN)){
    for(i in 1:length(rhos)){
      for(j in 1:n.runs){
        
        E<-getE(nT, nN[jj], rho.c = rhos[i], rho.r = rhos.r)
        
        temp.para<-t(matrix(rep(params[[ii]],times=nN[jj]),ncol=nN[jj]))
        y<-models[[ii]](E, temp.para)
        #plot.TL(y)
        
        tmp1<-cor(y,use = "pairwise.complete.obs")
        omega.spat<-mean(tmp1,na.rm=T)
        
        tmp1<-cor(t(y),use = "pairwise.complete.obs")
        omega.temp<-mean(tmp1,na.rm=T)
        
        Res$rho[a]<-rhos[i]
        Res$dim.ratio[a]<-nT/nN[jj]  # time / location
        Res$distribution[a]<-name.model[ii]
        Res$replicate[a]<-j
        
        z.spat<-tls(y)
        z.temp<-tls(t(y))
        Res[a,5:16]<-c(z.spat[1:6], z.temp[1:6])  # spatial + temporal (slope, intercept, p.quad, p.homo, RMSE, r)
        Res$p.fit.spat[a]<-z.spat[7]
        Res$p.fit.temp[a]<-z.temp[7]
        
        Res$syn.spat[a]<-omega.spat
        Res$syn.temp[a]<-omega.temp
        
        tmp<-skewness(y)*apply(y,2,mean)/apply(y,2,sd)
        Res$J[a]<-mean(tmp)
        y1<-as.vector(y);  Res$J.all[a]<-skewness(y1)/sd(y1)*mean(y1)
        
        a<-a+1
      }
    }
    show(paste(ii,jj))
  }
}
saveRDS(Res,"Res_test_temp_autocor_one_rhor.RDS")


#get average values for each block, and save as "res_simulation_average.RDS"
#get the index of good blocks (enough replicates which have TL hold), and save as "res_simulation_indexgood"
Res<-readRDS("Res_test_temp_autocor_one_rhor.RDS")
Fn_averagd_igood(Res, file.name1="res_test_temp_autocor_average_one_rhor", 
                 file.name2="res_test_temp_autocor_indexgood_one_rhor")


#output the figure of testing TL hold (Fig S4, save as "Fig_test_simulation.tiff") 
#and the barplot of TL slopes (Fig 1, save as "Fig_barplot_simulation.tiff")
res<-readRDS("res_test_temp_autocor_average_one_rhor.RDS")
i.good.comb<-readRDS("res_test_temp_autocor_indexgood_one_rhor.RDS")
Fn_plot_check_hold(res, file.name="Fig_test_temp_autocor_one_rhor")
Fn_plot_barplot(res, i.good.comb, nN, file.name="Fig_barplot_temp_autocorr_one_rhor")


#plot showing average temporal TL slopes vs. spatial TL slopes
#output: "Fig_simulation_temporal_vs._spatial_average.tif"
res<-readRDS("res_test_temp_autocor_average_one_rhor.RDS")
i.good.comb<-readRDS("res_test_temp_autocor_indexgood_one_rhor.RDS")
Fn_simulation_temp_vs_spat(res, i.good.comb)



#get R2 for different regression models, save as "ans_simulation_regression_R2.csv"
#meanwhile save the Results and average results of the good blocks as 
#"Res_simulation_good.RDS" and "res_simulation_average_good.RDS"
Res<-readRDS("Res_test_temp_autocor_one_rhor.RDS")
res<-readRDS("res_test_temp_autocor_average_one_rhor.RDS")
i.good.comb<-readRDS("res_test_temp_autocor_indexgood_one_rhor.RDS")
name.model<-dimnames(i.good.comb)[[1]]
dim.ratios<-dimnames(i.good.comb)[[2]]
rhos<-dimnames(i.good.comb)[[3]]
for(ii in 1:dim(i.good.comb)[1]){
  for(jj in 1:dim(i.good.comb)[2]){
    for(i in 1:dim(i.good.comb)[3]){
      if(is.na(i.good.comb[ii,jj,i])){
        tmp<-which(Res$distribution==name.model[ii] & Res$dim.ratio==dim.ratios[jj] & Res$rho==rhos[i])
        if(length(tmp)>0){Res<-Res[-tmp,]}
        
        tmp<-which(res$distribution==name.model[ii] & res$dim.ratio==dim.ratios[jj] & res$rho==rhos[i])
        if(length(tmp)>0){res<-res[-tmp,]}
      }
    }
  }
}
Res$omega1<-1/Res$syn.spat
res$omega1<-1/res$ave.syn.spat

saveRDS(Res,"Res_test_temp_autocor_one_rhor_good.RDS")
saveRDS(res,"res_test_temp_autocor_average_one_rhor_good.RDS")


Res1<-data.frame(slope.temp=Res$slope.temp, slope.spat=Res$slope.spat, distribution=Res$distribution,
                 omega=Res$syn.spat, J=Res$J, DR=Res$dim.ratio, omega.reciprocal=Res$omega1,rho=Res$rho)
res1<-data.frame(slope.temp=res$ave.slope.temp, slope.spat=res$ave.slope.spat, distribution=res$distribution,
                 omega=res$ave.syn.spat, J=res$ave.J, DR=res$dim.ratio, omega.reciprocal=res$omega1,rho=res$rho)

eqs.spat<-list(slope.spat~J+omega.reciprocal+DR, slope.spat~J*omega.reciprocal,
               slope.spat~J+omega.reciprocal, slope.spat~J+DR, slope.spat~omega.reciprocal+DR,
               slope.spat~J, slope.spat~omega.reciprocal, slope.spat~DR)
eqs.temp<-list(slope.temp~J+omega+DR, slope.temp~J*omega,
               slope.temp~J+omega, slope.temp~J+DR, slope.temp~omega+DR,
               slope.temp~J, slope.temp~omega, slope.temp~DR)

eqs<-list(eqs.spat, eqs.temp)
datas<-list(Res1,res1)
response<-c("spatial","temporal")
types<-c("all","average")
name.eqs<-c("J+omega+DR","J*omega","J+omega","J+DR","omega+DR","J","omega","DR")

res.regression<-data.frame(type=rep(NA,32), response=rep(NA,32), model=rep(NA,32), r2=rep(NA,32))
a<-1
for(i in 1:2){   # spatial, temporal, or difference
  for(j in 1:2){ # all data, averaged data
    for(jj in 1:8){
      z<-lm(eqs[[i]][[jj]],data=datas[[j]])
      res.regression$type[a]<-types[j]
      res.regression$response[a]<-response[i]
      res.regression$model[a]<-name.eqs[jj]
      res.regression$r2[a]<-round(summary(z)$r.squared,3)
      a<-a+1
    }
  }
}
write.csv(res.regression, "ans_test_temp_autocor_one_rhor.csv",row.names = F)




#final regression model
res.good<-readRDS("res_test_temp_autocor_average_one_rhor_good.RDS")
z<-lm(ave.slope.spat~ave.J+omega1, data=res.good)
summary(z)

z1<-z<-lm(ave.slope.temp~ave.J+ave.syn.spat, data=res.good)
summary(z1)

