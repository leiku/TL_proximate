#input: nT: number of time (rows); it has to be a single number
#       nN: number of locations (columns); it can be a vector 
#       rhos: a vector which contains the gradients of rho (controling spatial synchrony)
#       n.runs: number of replicates for each block
#       filename: the file name for the output RDS file
#output: a RDS file containing the calculated results

Fn_simulation<-function(nT, nN, rhos, n.runs, filename){
  n.all<-n.runs*5*length(nN)*length(rhos) # all entries
  
  dim.ratio<-nT/nN    
  name.model<-c("Poisson","Negative binomial","Exponential","Gamma","Lognormal")
  
  models<-list(model_NonIDPois,model_NonIDNegBinom,model_NonIDGamma,model_NonIDGamma,model_NonIDLognorm)
  #notice that the third one is actually Exponential distribution
  params<-list(c(1,1),c(5,0.4,1),c(1,1,1),c(4,1,1),c(1,1,1))
  
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
          mu<-numeric(nN[jj])
          Sigma<-matrix(rhos[i],nN[jj],nN[jj])+diag(1-rhos[i],nN[jj],nN[jj])
          E<-rmvnorm(nT,mu,Sigma,method='chol')
          
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
  saveRDS(Res,paste0(filename,".RDS"))
} 


