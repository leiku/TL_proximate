#input: Res: the RDS file got from Fn_simulation
#output: a RDS file containing averaged result for each block; 
#        a RDS file indicating the good blocks

Fn_averagd_igood<-function(Res, file.name1, file.name2){
  name.models<-unique(Res$distribution)
  dim.ratios<-unique(Res$dim.ratio)
  rhos<-unique(Res$rho)
  n.runs<-length(unique(Res$replicate))
  
  n.all<-length(name.models)*length(dim.ratios)*length(rhos)
  
  res<-data.frame(rho=rep(NA,n.all), dim.ratio=rep(NA,n.all), distribution=rep(NA,n.all),
                  frac.line.spat=rep(NA,n.all), frac.line.temp=rep(NA,n.all), frac.hete.spat=rep(NA,n.all),
                  frac.hete.temp=rep(NA,n.all), ave.RMSE.spat=rep(NA,n.all), sd.RMSE.spat=rep(NA,n.all),
                  ave.RMSE.temp=rep(NA,n.all), sd.RMSE.temp=rep(NA,n.all),ave.slope.spat=rep(NA,n.all),
                  sd.slope.spat=rep(NA,n.all),ave.slope.temp=rep(NA,n.all),sd.slope.temp=rep(NA,n.all),
                  ave.syn.spat=rep(NA,n.all),sd.syn.spat=rep(NA,n.all),ave.J=rep(NA,n.all),sd.J=rep(NA,n.all),
                  frac.fit.spat=rep(NA, n.all), frac.fit.temp=rep(NA, n.all))
  i.good.comb<-array(NA, c(length(name.models),length(dim.ratios),length(rhos)),dimnames=list(name.models,dim.ratios,rhos))
  i.good.comb1<-array(NA, c(length(name.models),length(dim.ratios),length(rhos)),dimnames=list(name.models,dim.ratios,rhos))
  a<-1
  for(ii in 1:length(name.models)){ # five models
    d1<-subset(Res, distribution==name.models[ii])
    for(jj in 1:length(dim.ratios)){
      d2<-subset(d1,dim.ratio==dim.ratios[jj])
      for(i in 1:length(rhos)){
        d3<-subset(d2, rho==rhos[i])
        
        res$rho[a]<-rhos[i]
        res$dim.ratio[a]<-dim.ratios[jj]
        res$distribution[a]<-name.models[ii]
        
        res$frac.line.spat[a]<-length(which(d3$p.line.spat<=0.01))/n.runs
        res$frac.line.temp[a]<-length(which(d3$p.line.temp<=0.01))/n.runs
        res$frac.hete.spat[a]<-length(which(d3$p.hete.spat<=0.01))/n.runs
        res$frac.hete.temp[a]<-length(which(d3$p.hete.temp<=0.01))/n.runs
        res$frac.fit.spat[a]<- length(which(d3$p.fit.spat>0.05))/n.runs #percent of bad fitting
        res$frac.fit.temp[a]<- length(which(d3$p.fit.temp>0.05))/n.runs 
        if(res$frac.line.spat[a]<0.1 & res$frac.line.temp[a]<0.1 & res$frac.hete.spat[a]<0.1 & res$frac.hete.temp[a]<0.1){
          i.good.comb[ii,jj,i]<-1
        }else{i.good.comb[ii,jj,i]<-NA}
        
        if(res$frac.line.spat[a]<0.1 & res$frac.line.temp[a]<0.1 & res$frac.hete.spat[a]<0.1 & res$frac.hete.temp[a]<0.1
           & res$frac.fit.spat[a]<0.1 & res$frac.fit.temp[a]<0.1){
          i.good.comb1[ii,jj,i]<-1
        }else{i.good.comb1[ii,jj,i]<-NA}
        
        
        res$ave.RMSE.spat[a]<-mean(d3$RMSE.spat,na.rm=T)
        res$sd.RMSE.spat[a]<-sd(d3$RMSE.spat,na.rm=T)
        res$ave.RMSE.temp[a]<-mean(d3$RMSE.temp,na.rm=T)
        res$sd.RMSE.temp[a]<-sd(d3$RMSE.temp,na.rm=T)
        res$ave.slope.spat[a]<-mean(d3$slope.spat,na.rm=T)
        res$sd.slope.spat[a]<-sd(d3$slope.spat,na.rm=T)
        res$ave.slope.temp[a]<-mean(d3$slope.temp,na.rm=T)
        res$sd.slope.temp[a]<-sd(d3$slope.temp,na.rm=T)
        res$ave.syn.spat[a]<-mean(d3$syn.spat,na.rm=T)
        res$sd.syn.spat[a]<-sd(d3$syn.spat,na.rm=T)
        res$ave.slope.diff[a]<-mean(d3$slope.temp-d3$slope.spat,na.rm=T)
        res$sd.slope.diff[a]<-sd(d3$syn.temp-d3$syn.spat,na.rm=T)
        res$ave.J[a]<-mean(d3$J,na.rm=T)
        res$sd.J[a]<-sd(d3$J,na.rm=T)
        
        a<-a+1
      }
    }
  }
  saveRDS(res,paste0(file.name1,".RDS"))
  saveRDS(i.good.comb,paste0(file.name2,".RDS"))
  saveRDS(i.good.comb1,paste0("test_",file.name2,".RDS"))
}



