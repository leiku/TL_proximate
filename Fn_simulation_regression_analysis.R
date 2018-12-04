
#input: Res: the RDS file got from Fn_simulation
#       res: the RDS file got from Fn_averaged_igood
#       i.good.comb: the RDS file got from Fn_averaged_igood
#output: a csv file showing R2 of different regression models

Fn_simulation_regression_analysis<-function(Res, res, i.good.comb, file.name){
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

  saveRDS(Res,"Res_simulation_good.RDS")
  saveRDS(res,"res_simulation_average_good.RDS")
  
  
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
  write.csv(res.regression, paste0(file.name, ".csv"),row.names = F)
}

