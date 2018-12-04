
#variation of TL slopes explained by different models

res<-readRDS("Res_empirical_TL_slopes_good.RDS")

slope.spat<-lapply(res, lapply, "[[", 1)
slope.spat<-lapply(slope.spat, unlist)

slope.temp<-lapply(res, lapply, "[[", 10)
slope.temp<-lapply(slope.temp, unlist)

omega<-lapply(res, lapply, "[[", 7)
omega<-lapply(omega, unlist)
receip<-function(x){1/x}
omega.rec<-lapply(omega, receip)

J<-lapply(res, lapply, "[[", 8)
J<-lapply(J, unlist)


N.matrix<-unlist(lapply(J, length))
dr<-c(35/11, 56/26, 28/19, 28/12, 28/12, 28/12)
DR<-J
for(i in 1:6){DR[[i]]<-rep(dr[i], N.matrix[i])}


###### all six datasets together   ########
slope.spat1<-unlist(slope.spat)
slope.temp1<-unlist(slope.temp)
J1<-unlist(J)
omega.rec1<-unlist(omega.rec)
omega1<-unlist(omega)
DR1<-unlist(DR)


models<-c("slope.spat1~J1+omega.rec1+DR1", "slope.spat1~J1+omega.rec1",
          "slope.spat1~J1+DR1", "slope.spat1~omega.rec1+DR1",
          "slope.spat1~J1", "slope.spat1~omega.rec1","slope.spat1~DR1")
models1<-c("slope.temp1~J1+omega1+DR1", "slope.temp1~J1+omega1",
          "slope.temp1~J1+DR1", "slope.temp1~omega1+DR1",
          "slope.temp1~J1", "slope.temp1~omega1","slope.temp1~DR1")
res.spat<-data.frame(type=rep("spatial", length(models)), models=models, 
                     r2=rep(0, length(models)), AIC=rep(0, length(models)))
res.temp<-data.frame(type=rep("temporal", length(models1)), models=models1, 
                     r2=rep(0, length(models1)), AIC=rep(0, length(models1)))

for(i in 1:length(models)){
  z<-lm(as.formula(models[i]))
  z1<-summary(z)
  res.spat$r2[i]<-z1$r.squared
  res.spat$AIC[i]<-AIC(z)
  
  z<-lm(as.formula(models1[i]))
  z1<-summary(z)
  res.temp$r2[i]<-z1$r.squared
  res.temp$AIC[i]<-AIC(z)
}

tmp1<-min(res.spat$AIC)
res.spat$det.AIC<-res.spat$AIC-tmp1
res.spat$weight<-exp(-0.5*(res.spat$AIC-tmp1))
res.spat$weight<-res.spat$weight/sum(res.spat$weight)

tmp1<-min(res.temp$AIC)
res.temp$det.AIC<-res.temp$AIC-tmp1
res.temp$weight<-exp(-0.5*(res.temp$AIC-tmp1))
res.temp$weight<-res.temp$weight/sum(res.temp$weight)


write.csv(rbind(res.spat, res.temp),"res_empirical_regression_all.csv",row.names = F)



######### each dataset  ##########

models<-c("slope.spat2~J2*omega.rec2", "slope.spat2~J2+omega.rec2",
          "slope.spat2~J2", "slope.spat2~omega.rec2")
models1<-c("slope.temp2~J2*omega2", "slope.temp2~J2+omega2",
           "slope.temp2~J2", "slope.temp2~omega2")
Res.spat<-data.frame(type="spatial", models="delete", r2=0, AIC=0, 
                     det.AIC=0, weight=0)
Res.temp<-Res.spat

for(j in 1:6){
  slope.spat2<-slope.spat[[j]]
  slope.temp2<-slope.temp[[j]]
  J2<-J[[j]]
  omega.rec2<-omega.rec[[j]]
  omega2<-omega[[j]]
  
  res.spat<-data.frame(type=rep("spatial", length(models)), models=models, 
                       r2=rep(0, length(models)), AIC=rep(0, length(models)))
  res.temp<-data.frame(type=rep("temporal", length(models1)), models=models1, 
                       r2=rep(0, length(models1)), AIC=rep(0, length(models1)))
  for(i in 1:length(models)){
    z<-lm(as.formula(models[i]))
    z1<-summary(z)
    res.spat$r2[i]<-z1$r.squared
    res.spat$AIC[i]<-AIC(z)
    
    z<-lm(as.formula(models1[i]))
    z1<-summary(z)
    res.temp$r2[i]<-z1$r.squared
    res.temp$AIC[i]<-AIC(z)
  }
  
  tmp1<-min(res.spat$AIC)
  res.spat$det.AIC<-res.spat$AIC-tmp1
  res.spat$weight<-exp(-0.5*(res.spat$AIC-tmp1))
  res.spat$weight<-res.spat$weight/sum(res.spat$weight)
  
  tmp1<-min(res.temp$AIC)
  res.temp$det.AIC<-res.temp$AIC-tmp1
  res.temp$weight<-exp(-0.5*(res.temp$AIC-tmp1))
  res.temp$weight<-res.temp$weight/sum(res.temp$weight)
  
  Res.spat<-rbind(Res.spat, res.spat)
  Res.temp<-rbind(Res.temp, res.temp)
}
Res.spat<-Res.spat[-1,]
Res.temp<-Res.temp[-1,]
Res.spat$dataset<-rep(names(slope.spat), each=length(models))
Res.temp$dataset<-rep(names(slope.temp), each=length(models))

write.csv(rbind(Res.spat, Res.temp),"res_empirical_regression_each.csv",row.names = F)



