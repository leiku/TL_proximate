rm(list=ls())

library(lmodel2)
source("Fn_tls.R")

#get clean Data, save as "data_original.RDS",
#as well as the name of species, save as "species_names.csv"

source('empirical_data.R')


#plot count distributions for each dataset
#input: "data_original.RDS"
#output: Figs of the distribution of counts: FigS_distribution_(name of datasets).tif
source('empirical_hist_count.R')


#get TL slopes, as well as J and omega
#input: "data_original.RDS"
#output: "Res_empirical_TL_slopes.RDS"
source("empirical_TL_slopes.R")


#TL scatter plot for each dataset
#input: "Res_empirical_TL_slopes.RDS"
#output: "Figs_empirical_scatter_spatial_(dataset name).tif"
#        "Figs_empirical_scatter_temporal_(dataset name).tif"
source("empirical_scatter_plot.R")



#omit the matrix which rejected the TL tests
res<-readRDS("Res_empirical_TL_slopes.RDS")

for(i in 1:6){
  res1<-res[[i]]
  ii<-rep(NA, length(res1))
  for(j in 1:length(res1)){
    if(all(res1[[j]][3:4,]>0.01)){ii[j]<-1}
  }
  res[[i]]<-res1[!is.na(ii)]
}
saveRDS(res,"Res_empirical_TL_slopes_good.RDS")


#plot temporal slope vs. spatial slope 
#input: "Res_empirical_TL_slopes_good.RDS"
#output: "Fig2_empirical_temp_vs_spat.tif"
source("empirical_plot_temp_vs_spat.R")


#get regression results for different models explaining variation of TL slopes
#input: "Res_empirical_TL_slopes_good.RDS"
#output: "res_empirical_regression_all.csv" for all six datasets
#        "res_empirical_regression_each.csv" for each dataset separately
source("empirical_regression.R")


#predict empirical TL slopes using theoretical coefficients
#input: "Res_empirical_TL_slopes.RDS", coefficients
#output: "Fig_empirical_predict.tif"
#        "res_empirical_predict.csv"
source("empirical_predict.R")



 
 
 
 
 
 
 
 
 
 
 #################
### spat1: destroy spatial synchrony
### temp1: destroy temporal synchrony
### spat2 and temp2: destroy both

res<-vector("list", 6)   # dataset -> spat1/temp1/spat2/temp2 -> slope/p.lm/p.q/p.het -> spec*n.perm
names(res)<-c("aphid","north.sea","calcofi.gr1","calcofi.gr2","calcofi.gr3","calcofi.gr4")
n.spec<-c(20,22,rep(10,4))
n.surr<-1
for(i in 1:6){
  res.spat<-matrix(NA, 6, n.spec[i])
  row.names(res.spat)<-c("slope", "intercept", "p_curve", "p_het", "r","syn")
  colnames(res.spat)<-names(D[[i]])
  res.temp<-res.spat
  res.sur.temp<-array(NA, c(6, n.spec[i], n.surr), 
                      dimnames = list(c("slope", "intercept", "p_curve", "p_het", "r","syn"),NULL,NULL))
  res.sur.spat<-res.sur.temp
  res.sur.all.spat<-res.sur.temp
  res.sur.all.temp<-res.sur.temp
  for(sp in 1:n.spec[i]){
    ## output of tls function
    x<-D[[i]][[sp]]

    flags<-"lm"
    res.spat[,sp]<-tls(t(x), fitting=flags)
    res.temp[,sp]<-tls(x, fitting=flags)
    res.sur.spat[,sp,]<-tls.surrogate(t(x),n.surr, fitting=flags)
    res.sur.temp[,sp,]<-tls.surrogate(x,n.surr, fitting=flags)
    tmp<-tls.surrogate.all(t(x),n.surr, fitting=flags)

    res.sur.all.spat[,sp,]<-tmp[[1]]
    res.sur.all.temp[,sp,]<-tmp[[2]]
    
    show(sp)
  }
  res1<-list(res.spat, res.temp, res.sur.spat, res.sur.temp, res.sur.all.spat, res.sur.all.temp)
  names(res1)<-c("spat", "temp", "sur.spat1", "sur.temp1", "sur.spat2", "sur.temp2")
  res[[i]]<-res1
}

saveRDS(res,"Res_all_spat_vs_temp.RDS")



########## analysis (different factors) ##########
library(asbio)
res<-readRDS("Res_all_spat_vs_temp.RDS")

#find population matrices in which TL holds
II<-matrix(NA,6,3)
colnames(II)<-c("spat","temp","both")
for(i.data in 1:6){
  x<-res[[i.data]]
  y1<-x[[2]] #spatial
  y2<-x[[1]] #temporal
  
  #spatial
  II[i.data,1]<-length(which(y1[3,]>0.01 & y1[4,]>0.01))
  II[i.data,2]<-length(which(y2[3,]>0.01 & y2[4,]>0.01))
  II[i.data,3]<-length(which(y1[3,]>0.01 & y1[4,]>0.01 & y2[3,]>0.01 & y2[4,]>0.01))
}

n.spec<-c(20,22,rep(10,4))
name.data<-names(res)
dat.analysis<-data.frame(spat.slope=rep(NA,82), temp.slope=rep(NA,82), sk.cv.all=rep(NA,82), sk.cv.each=rep(NA,82), 
                         syn=rep(NA,82), dim.ratio=rep(NA,82), dataset=rep(NA,82))
a<-1
for(i.data in 1:6){
  x<-res[[i.data]]
  for(j in 1:n.spec[i.data]){
    dat.analysis$spat.slope[a]<-x[[1]][1,j]
    dat.analysis$temp.slope[a]<-x[[2]][1,j]
    
    y<-D[[i.data]][[j]]
    y<-t(y)
    y1<-as.vector(y)
    y1<-y1[!is.na(y1)]
    dat.analysis$sk.cv.all[a]<-skew(y1)/sd(y1)*mean(y1)
    
    sk.cv2<-rep(NA,ncol(y))
    for(jj in 1:ncol(y)){
      y1<-y[,jj]
      y1<-y1[!is.na(y1)]
      sk.cv2[jj]<-skew(y1)/sd(y1)*mean(y1)
    }
    dat.analysis$sk.cv.each[a]<-mean(sk.cv2,na.rm=T)
    
    dat.analysis$syn[a]<-mean(cor(y,use="pairwise.complete.obs"),na.rm=T)
    dat.analysis$dim.ratio[a]<-ncol(y)/nrow(y) # site / year
    dat.analysis$dataset[a]<-name.data[i.data]
    a<-a+1
  }
}
saveRDS(dat.analysis,"Res_all_skewness_syn.RDS")

dat.analysis<-readRDS("Res_all_skewness_syn.RDS")
dat.analysis$omega1<-1/dat.analysis$syn
nn<-6*2*3
dat.analysis$diff.slope<-dat.analysis$temp.slope-dat.analysis$spat.slope
res.aic<-data.frame(dataset=rep(NA,nn), response=rep(NA,nn), factor=rep(NA,nn), 
                    p.all=rep(NA,nn), r2=rep(NA,nn),
                    AIC=rep(NA,nn), detAIC=rep(NA,nn), weight=rep(NA,nn))
res.cov<-data.frame(dataset=rep(NA,nn), response=rep(NA,nn), factor=rep(NA,nn), 
                    cov1=rep(NA,nn), cov2=rep(NA,nn), cov3=rep(NA,nn), cov4=rep(NA,nn))
responses<-c("spatial","temporal")
i.resp<-c(1,2)
for(j in 1:2){
  for(i in 1:6){
    x<-subset(dat.analysis,dataset==name.data[i])
    if(j==1){jj=8}else{jj=5}
    z1<-lm(x[,i.resp[j]]~x[,4]+x[,jj])
    z2<-lm(x[,i.resp[j]]~x[,4])
    z3<-lm(x[,i.resp[j]]~x[,jj])
    z10<-summary(z1)

    a<-((3*i-2):(3*i))+3*6*(j-1)
    res.aic$dataset[a]<-name.data[i]
    res.aic$response[a]<-responses[j]
    res.aic$factor[a]<-c("J+omega","J","omega")
    res.aic$p.all[a]<-round(c(pf(z10$fstatistic[1],z10$fstatistic[2],z10$fstatistic[3],lower.tail = F),
                              summary(z2)$coef[2,4],summary(z3)$coef[2,4]),3)
    res.aic$r2[a]<-round(c(summary(z1)$r.squared, summary(z2)$r.squared, summary(z3)$r.squared),3)
    aics<-AIC(z1,z2,z3)[,2]
    res.aic$AIC[a]<-round(aics,3)
    tmp<-min(aics)
    res.aic$detAIC[a]<-round(aics-tmp,3)
    weights<-exp(-0.5*(aics-tmp))
    res.aic$weight[a]<-round(weights/sum(weights),3)
  }
}
write.csv(res.aic,"Res_model_selection_skew_each.csv",row.names = F)


# summary
res.aic<-read.csv("Res_model_selection_skew_each.csv")
res.sum<-data.frame(response=rep(NA,12), dataset=rep(NA,12), model=rep(NA,12), r2=rep(NA,12), p=rep(NA,12),
                    importance.J=rep(NA,12), importance.omega=rep(NA,12))
for(i in 1:12){
  tmp<-res.aic[(3*i-2):(3*i),]
  res.sum$response[i]<-as.character(tmp$response[1])
  res.sum$dataset[i]<-as.character(tmp$dataset[1])
  ii<-which.min(tmp$detAIC)
  res.sum$model[i]<-as.character(tmp$factor[ii])
  res.sum[i,4:5]<-tmp[ii,c(5,4)]
  res.sum$importance.J[i]<-sum(tmp$weight[1:2])
  res.sum$importance.omega[i]<-sum(tmp$weight[c(1,3)])
}
write.csv(res.sum,"Res_model_selection_summary_all.csv",row.names = F)


# cov
dat.analysis<-readRDS("Res_all_skewness_syn.RDS")
dat.analysis$omega1<-1/dat.analysis$syn
nn<-6
res.cov<-data.frame(dataset=rep(NA,nn),cov.J=rep(NA,nn), cov.omega=rep(NA,nn), cov.plus=rep(NA,nn), cov.prod=rep(NA,nn))
for(i in 1:6){
  x<-subset(dat.analysis,dataset==name.data[i])
  D<-sd(x[,1])*sd(x[,2])
  
  z1<-lm(x[,1]~x[,4])
  z2<-lm(x[,2]~x[,4])
  z11<-predict(z1)
  z21<-predict(z2)
  
  res.cov$dataset[i]<-name.data[i]
  res.cov$cov.J[i]<-cov(z11,z21)/D
  
  z1<-lm(x[,1]~x[,8])
  z2<-lm(x[,2]~x[,5])
  z11<-predict(z1)
  z21<-predict(z2)
  res.cov$cov.omega[i]<-cov(z11,z21)/D
  
  z1<-lm(x[,1]~x[,8]+x[,4])
  z2<-lm(x[,2]~x[,5]+x[,4])
  z11<-predict(z1)
  z21<-predict(z2)
  res.cov$cov.plus[i]<-cov(z11,z21)/D
  
  z1<-lm(x[,1]~x[,8]*x[,3])
  z2<-lm(x[,2]~x[,5]*x[,3])
  z11<-predict(z1)
  z21<-predict(z2)
  res.cov$cov.prod[i]<-cov(z11,z21)/D
}
write.csv(res.cov,"Res_model_selection_cov_each.csv",row.names = F)

############### hist ####################
dat.analysis<-readRDS("Res_all_skewness_syn.RDS")
tiff(paste0("FigS_hist_slopes.tif"),width=6, height=3, units="in",compression="zip",res=600)
op<-par(mfrow=c(1,2),oma=c(1,1,1,1), mar=c(2.5,2.5,0,0),mgp=c(1.5,0.6,0),pty="s")
hist(dat.analysis$spat.slope,15,xlab="spatial TL slopes",ylab="count", main=NA)
mtext(paste0("(a)"), side=3, line=-1.5, adj=0.05, cex=1)
hist(dat.analysis$temp.slope,15,xlab="temporal TL slopes",ylab="count", main=NA)
mtext(paste0("(b)"), side=3, line=-1.5, adj=0.05, cex=1)
par(op)
dev.off()
########### analysis all datasets ###############

dat.analysis<-readRDS("Res_all_skewness_syn.RDS")
res1<-data.frame(slope.spat=dat.analysis$spat.slope, slope.temp=dat.analysis$temp.slope, 
                 J=dat.analysis$sk.cv.each, omega=dat.analysis$syn, DR=dat.analysis$dim.ratio,
                 omega1=1/dat.analysis$syn)
eqs.spat<-list(slope.spat~J+omega1+DR, 
               slope.spat~J+omega1, slope.spat~J+DR, slope.spat~omega1+DR,
               slope.spat~J, slope.spat~omega1, slope.spat~DR)
eqs.temp<-list(slope.temp~J+omega+DR,
               slope.temp~J+omega, slope.temp~J+DR, slope.temp~omega+DR,
               slope.temp~J, slope.temp~omega, slope.temp~DR)

eqs<-list(eqs.spat, eqs.temp)
response<-c("spatial","temporal")
name.eqs<-c("J+omega+DR","J+omega","J+DR","omega+DR","J","omega","DR")

res.regression<-data.frame(response=rep(NA,14), model=rep(NA,14), r2=rep(NA,14), p=rep(NA,14), 
                           AIC=rep(NA,14), det.AIC=rep(NA,14), weight=rep(NA,14),weight.cor=rep(NA,14))
a<-1
for(i in 1:2){   # spatial, temporal, 
  AICs<-rep(0,7)
  for(j in 1:7){
    z<-lm(eqs[[i]][[j]],data=res1)
    z1<-summary(z)
    res.regression$response[a]<-response[i]
    res.regression$model[a]<-name.eqs[j]
    res.regression$r2[a]<-round(summary(z)$r.squared,3)
    res.regression$p[a]<-round(pf(z1$fstatistic[1],z1$fstatistic[2],z1$fstatistic[3],lower.tail = F),3)
    AICs[j]<-AIC(z)
    res.regression$AIC[a]<-round(AICs[j],3)
    
    a<-a+1
  }
  b<-(i-1)*7+(1:7)
  min.AIC<-min(AICs)
  res.regression$det.AIC[b]<-round(AICs-min.AIC,3)
  res.regression$weight[b]<-round(exp(-0.5*(AICs-min.AIC)),3)
  res.regression$weight.cor[b]<-round(res.regression$weight[b]/sum(res.regression$weight[b]),3)
}
write.csv(res.regression,"ans_regression_factors_all.csv",row.names = F)


######### Fig 1 temp vs spat (empirical) ##########

res<-readRDS("Res_all_spat_vs_temp.RDS")
n.spec<-c(20,22,rep(10,4))

axis.min<-c(1.8, 1.1, 1.5, 1.4, -0.5, -0.5)
axis.max<-c(3, 2.7, 2.9, 6.5, 3.2, 4.1)
tiff(paste0("Fig1_empirical_spat_vs_temp.tif"),width=3.22, height=4.5, units="in",compression="zip",res=300)
op<-par(mfrow=c(3,2),oma=c(3,3,0,1), mar=c(1.5,1.5,1,0),mgp=c(1,0.6,0),pty="s")

for(i.data in 1:6){
  x<-res[[i.data]]
  
  spat.slope<-x[[1]][1,1:n.spec[i.data]]
  temp.slope<-x[[2]][1,1:n.spec[i.data]]

  z<-lmodel2(temp.slope~spat.slope,data=x,nperm=99)
  z.res<-as.numeric(z$regression.results[2,])
  x.fitted<-c(min(spat.slope),max(spat.slope))
  y.fitted<-z.res[2]+z.res[3]*x.fitted
  
  y.range<-range(temp.slope)
  if(i.data<5){
    plot(temp.slope~spat.slope, pch=1, ann=F, ylim=c(y.range[1]-0.2*(y.range[2]-y.range[1]),y.range[2]))
  }else{plot(temp.slope~spat.slope, pch=1, ann=F)}
  #plot(temp.slope~spat.slope, pch=1, ann=F, xlim=c(axis.min[i.data],axis.max[i.data]),ylim=c(axis.min[i.data],axis.max[i.data]))
  
  if(z.res[5]>0.05){
    lines(x.fitted, y.fitted,lty="dashed")
  }else{lines(x.fitted, y.fitted,lty="solid")}
  
  mtext(paste0("(",letters[i.data],")"), side=3, line=-1.5, adj=0.05, cex=0.8)
  mtext(paste0("slope=",round(z.res[3],2)), side=1, line=-1.8, adj=0.95, cex=0.7)
  mtext(paste0("r=",round(cor(spat.slope,temp.slope),2)), side=1, line=-1, adj=0.95, cex=0.7)
}
# xlab + ylab
par(fig = c(0, 1, 0, 1), oma=c(0,0,0,0), mar = c(3, 3, 0, 0), mgp=c(2,1,0), new = TRUE)
plot(NA, xlim=c(0,1),ylim=c(0,1), ann=F,
     type = "n", bty = "n", xaxt = "n", yaxt = "n", cex.lab=1, font.lab=2)
mtext("spatial TL slope", side=1, line=6, cex=0.8)
mtext("temporal TL slope", side=2, line=0.6,cex=0.8)
par(op)

dev.off()



## all


tiff(paste0("Fig_ppt_empirical_spat_vs_temp.tif"),width=6, height=6, units="in",compression="zip",res=300)
name.data=c("aphid","north.sea","calcofi.gr1","calcofi.gr2","calcofi.gr3","calcofi.gr4")
plot(NA, xlim=range(dat.analysis$spat.slope),ylim=range(dat.analysis$temp.slope),xlab="spatial TL slope",ylab="temporal TL slope")
cols<-rainbow(6)
for(i in 1:6){
  points(temp.slope~spat.slope,data=subset(dat.analysis,dataset==name.data[i]),col=cols[i])
}
legend("topleft",name.data,col=cols,pch=1)
dev.off()

######### Fig 2 temp vs spat (permutation) ##########

res<-readRDS("Res_all_spat_vs_temp.RDS")
n.spec<-c(20,22,rep(10,4))
#axis.min<-c(1.7, 1.3, 1.5, 1.5, -0.5, -0.5)
#axis.max<-c(6, 5, 4, 6.5, 4.5, 5.5)
tiff(paste0("Fig1_permutation_spat_vs_temp.tif"),width=6, height=9, units="in",compression="zip",res=300)
op<-par(mfrow=c(3,2),oma=c(5,3,0,1), mar=c(3,3,1,0),mgp=c(2,1,0),pty="s")

for(i.data in 1:6){
  x<-res[[i.data]]
  
  tmp<-x[[5]][1,,]
  slope.temp<-apply(tmp[1:n.spec[i.data],],1,mean,na.rm=T)
  SD.temp<-apply(tmp[1:n.spec[i.data],],1,sd,na.rm=T)
  
  tmp<-x[[6]][1,,]
  slope.spat<-apply(tmp[1:n.spec[i.data],],1,mean,na.rm=T)
  SD.spat<-apply(tmp[1:n.spec[i.data],],1,sd,na.rm=T)

  z<-lmodel2(slope.temp~slope.spat,data=x,nperm=99)
  z.res<-as.numeric(z$regression.results[2,])
  x.fitted<-c(min(slope.spat-SD.spat),max(slope.spat+SD.spat))
  y.fitted<-z.res[2]+z.res[3]*x.fitted
  
  plot(slope.temp~slope.spat, pch=15, ann=F, xlim=range(c(slope.spat-SD.spat,slope.spat+SD.spat)),
       ylim=range(c(slope.temp-SD.temp, slope.temp+SD.temp)))
  if(z.res[5]>0.05){
    lines(x.fitted, y.fitted,lty="dashed")
  }else{lines(x.fitted, y.fitted,lty="solid")}
  
  arrows(slope.spat, slope.temp-SD.temp, slope.spat, slope.temp+SD.temp,
         length=0.03, angle=90, code=3,col="grey")
  arrows(slope.spat-SD.spat, slope.temp, slope.spat+SD.spat, 
         slope.temp, length=0.03, angle=90, code=3, col="grey")
  points(slope.temp~slope.spat,pch=15)
  
  mtext(paste0("(",letters[i.data],")"), side=3, line=-1.5, adj=0.05, cex=1)
  mtext(paste0("slope = ",round(z.res[3],2)), side=1, line=-2.2, adj=0.95, cex=0.9)
  mtext(paste0("r = ",round(cor(slope.spat,slope.temp),2)), side=1, line=-1, adj=0.95, cex=0.9)
}
# xlab + ylab
par(fig = c(0, 1, 0, 1), oma=c(0,0,0,0), mar = c(3, 3, 0, 0), mgp=c(2,1,0), new = TRUE)
plot(NA, xlim=c(0,1),ylim=c(0,1), ann=F,
     type = "n", bty = "n", xaxt = "n", yaxt = "n", cex.lab=1.2, font.lab=2)
mtext("spatial TL slope", side=1, line=11, cex=1)
mtext("temporal TL slope", side=2, line=0.5,cex=1)
par(op)

dev.off()


####### CIs of slopes  ###########
res<-readRDS("Res_all_spat_vs_temp.RDS")
n.spec<-c(20,22,rep(10,4))
ci.slope<-matrix(NA,6,4)
for(i.data in 1:6){
  x<-res[[i.data]]
  n.surr<-dim(x[[3]])[3]
  slopes<-vector("list",6) # spat, temp, spat1, temp1, spat2, temp2
  CIs<-vector("list",6)
  for(i in 1:6){
    if(i<=2){
      slopes[[i]]<-x[[i]][1,1:n.spec[i.data]]
      CIs[[i]]<-NA
    }else{
      tmp<-x[[i]][1,,]
      tmp.p<-x[[i]][2,,]
      tmp[tmp.p>0.05]<-NA
      slopes[[i]]<-apply(tmp[1:n.spec[i.data],],1,mean,na.rm=T)
      CIs[[i]]<-apply(tmp[1:n.spec[i.data],],1,sd,na.rm=T)/sqrt(n.surr)
    }
  }
  
  index<-matrix(c(1,5,2,6), 2, 2)
  
  for(i in 1:2){
    spat.slope<-slopes[[index[i,1]]]
    temp.slope<-slopes[[index[i,2]]]
    spat.ci<-CIs[[index[i,1]]]
    temp.ci<-CIs[[index[i,2]]]
    
    if(i==1){
      ii<-which(res[[i.data]][[1]][2,]<0.05 & res[[i.data]][[2]][2,]<0.05)
    }else{ii<-1:length(spat.slope)}
    pchs<-rep(0,length(spat.slope))
    pchs[ii]<-15
    
    z<-lmodel2(temp.slope[ii]~spat.slope[ii],data=x,nperm=999)
    ci.slope[i.data, ((i-1)*2+1):(i*2)]<-as.numeric(z$confidence.intervals[2,c(4,5)])
  }
}
write.csv(ci.slope,"res_ci_slope.csv")


################# cor
dat.analysis<-readRDS("Res_all_skewness_syn.RDS")
name.data<-c("aphid","north.sea","calcofi.gr1","calcofi.gr2","calcofi.gr3","calcofi.gr4")
res1<-data.frame(slope.spat=dat.analysis$spat.slope, slope.temp=dat.analysis$temp.slope, 
                 slope.diff=dat.analysis$temp.slope-dat.analysis$spat.slope,
                 J=dat.analysis$sk.cv.each, omega=dat.analysis$syn, DR=dat.analysis$dim.ratio,
                 omega1=1/dat.analysis$syn,dataset=dat.analysis$dataset)
cor.test(res1$slope.spat,res1$slope.temp,method="spearman")
z1<-lm(slope.spat~J,data=res1)
z2<-lm(slope.temp~J,data=res1)
cor.test(residuals(z1),residuals(z2),method="spearman")
plot(residuals(z1),residuals(z2))
z3<-lm(residuals(z1)~res1$omega1)
z4<-lm(residuals(z2)~res1$omega)
cor.test(residuals(z3),residuals(z4),method="spearman")
plot(residuals(z3),residuals(z4))

#each
cor.rho<-matrix(NA,6,3)
colnames(cor.rho)<-c("rho.original","rho.removingJ","rho.removingOmega")
cor.p<-cor.rho
colnames(cor.p)<-c("p.original","p.removingJ","p.removingOmega")
for(i in 1:6){
  d<-subset(res1,dataset==name.data[i])
  z0<-cor.test(d$slope.spat,d$slope.temp,method="spearman")
  cor.rho[i,1]<-z0$estimate
  cor.p[i,1]<-z0$p.value
  
  z1<-lm(slope.spat~J,data=d)
  z2<-lm(slope.temp~J,data=d)
  z0<-cor.test(residuals(z1),residuals(z2),method="spearman")
  cor.rho[i,2]<-z0$estimate
  cor.p[i,2]<-z0$p.value
  
  z3<-lm(residuals(z1)~d$omega1)
  z4<-lm(residuals(z2)~d$omega)
  z0<-cor.test(residuals(z3),residuals(z4),method="spearman")
  cor.rho[i,3]<-z0$estimate
  cor.p[i,3]<-z0$p.value
}
ans<-cbind(cor.rho,cor.p)
rownames(ans)<-name.data
write.csv(ans,"cor_redisuals.csv",row.names = T)

#plot
op<-par(mfrow=c(2,2))
i=5
d.test<-subset(res1,dataset==name.data[i])
z1<-lm(slope.spat~J,d.test)
z2<-lm(slope.temp~J,d.test)
plot(slope.temp~slope.spat,d.test)
points(predict(z1),predict(z2),col="red")
plot(residuals(z2)~residuals(z1))
plot(slope.spat~J,d.test)
points(d.test$J,predict(z1),col="red")
plot(slope.temp~J,d.test)
points(d.test$J,predict(z2),col="red")


################ skewness ###############  useless
library(asbio)
res<-readRDS("Res_all_spat_vs_temp.RDS")
n.spec<-c(20,22,rep(10,4))
R1<-data.frame(skew=rep(NA,82), CV=rep(NA,82), data=rep(NA,82),letter=rep(NA,82), prop=rep(NA,82), 
               temp=rep(NA,82), spat=rep(NA,82),temp.sur=rep(NA,82), spat.sur=rep(NA,82),
               syn.temp=rep(NA,82),syn.spat=rep(NA,82))
a<-1
R2<-data.frame(data=rep(NA,8200), temp=rep(NA,8200), spat=rep(NA,8200))
b<-0
for(i.dataset in 1:6){
  for(s in 1:length(D[[i.dataset]])){
    x<-D[[i.dataset]][[s]]
    R1$syn.temp[a]<-mean(cor(x,use = "pairwise.complete.obs"), na.rm=T)
    R1$syn.spat[a]<-mean(cor(t(x),use = "pairwise.complete.obs"), na.rm=T)
    
    x<-as.vector(x)
    x<-x[!is.na(x)]
    R1$skew[a]<-as.numeric(skew(x))
    R1$CV[a]<-sd(x)/mean(x)
    R1$data[a]<-names(D)[i.dataset]
    R1$letter[a]<-letters[s]
    
    R1$spat[a]<-res[[i.dataset]][[1]][1,s]
    R1$temp[a]<-res[[i.dataset]][[2]][1,s]
    
    y1<-res[[i.dataset]][[5]][2,s,]
    y2<-res[[i.dataset]][[6]][2,s,]
    z<-which(y1<0.05&y2<0.05)
    R1$prop[a]<-length(z)/length(y1)
    
    z1<-res[[i.dataset]][[5]][1,s,]
    z2<-res[[i.dataset]][[6]][1,s,]
    R2$data[(b+1):(b+length(z1))]<-names(D)[i.dataset]
    R2$spat[(b+1):(b+length(z1))]<-z1
    R2$temp[(b+1):(b+length(z1))]<-z2
    b<-b+length(z1)
    
    R1$spat.sur[a]<-mean(z1)
    R1$temp.sur[a]<-mean(z2)

    a<-a+1
  }
}
R2<-R2[-(b:8200),]

R1$b<-R1$skew/R1$CV

# spatial
names.data<-c("aphid","north.sea","calcofi.gr1","calcofi.gr2","calcofi.gr3","calcofi.gr4")
y.names<-c("spatial TL slope",  "spatial TL slope after permutation")
x.names<-c("skewness/CV",expression(Omega))
for(i in 1:6){
  r1<-subset(R1,data==names.data[i])
  tiff(paste0("Fig_skewness_spat_",i,"_",names.data[i],".tif"),width=8, height=8, units="in",compression="zip",res=300)
  op<-par(mfrow=c(2,2),oma=c(3,3,0,1), mar=c(3.5,3.5,2,2),mgp=c(2.5,1,0))
  for(jj in 1:2){
    y<-r1[,5+2*jj]   # spat, spat.surr
    for(j in 1:2){
      x<-r1[,13-j]   # skew/cv, omega
      plot(y~x,xlab=x.names[j],ylab=y.names[jj])
      z<-lm(y~x)
      z.fitted<-fitted(z)
      if(summary(z)$coef[2,4]>0.05){ltys<-"dashed"}else{ltys<-"solid"}
      tmp<-c(which.min(x),which.max(x))
      lines(x[tmp],z.fitted[tmp], lty=ltys)
      mtext(paste0("(",letters[j+(jj-1)*2],")"), side=3, line=-1.5, adj=0.05, cex=1)
      mtext(paste0("slope=",round(coef(z)[2],3)), side=1, line=-3, adj=0.95, cex=1)
      mtext(paste0("r=",round(cor(x,y),3)), side=1, line=-1.5, adj=0.95, cex=1)
    }
  }
  par(op)
  dev.off()
}


# temporal
names.data<-c("aphid","north.sea","calcofi.gr1","calcofi.gr2","calcofi.gr3","calcofi.gr4")
y.names<-c("temporal TL slope",  "temporal TL slope after permutation")
x.names<-c("skewness/CV",expression(Omega))
for(i in 1:6){
  r1<-subset(R1,data==names.data[i])
  tiff(paste0("Fig_skewness_temp_",i,"_",names.data[i],".tif"),width=8, height=8, units="in",compression="zip",res=300)
  op<-par(mfrow=c(2,2),oma=c(3,3,0,1), mar=c(3.5,3.5,2,2),mgp=c(2.5,1,0))
  for(jj in 1:2){
    y<-r1[,4+2*jj]   # spat, spat.surr
    for(j in 1:2){
      x<-r1[,13-j]   # skew/cv, omega
      plot(y~x,xlab=x.names[j],ylab=y.names[jj])
      z<-lm(y~x)
      z.fitted<-fitted(z)
      if(summary(z)$coef[2,4]>0.05){ltys<-"dashed"}else{ltys<-"solid"}
      tmp<-c(which.min(x),which.max(x))
      lines(x[tmp],z.fitted[tmp], lty=ltys)
      mtext(paste0("(",letters[j+(jj-1)*2],")"), side=3, line=-1.5, adj=0.05, cex=1)
      mtext(paste0("slope=",round(coef(z)[2],3)), side=1, line=-3, adj=0.95, cex=1)
      mtext(paste0("r=",round(cor(x,y),3)), side=1, line=-1.5, adj=0.95, cex=1)
    }
  }
  par(op)
  dev.off()
}

############## ratio ###################
res<-readRDS("Res_all_spat_vs_temp.RDS")
n.spec<-c(20,22,rep(10,4))
R1<-data.frame(data=rep(NA, 82), ratio<-rep(NA,82))
a<-0
for(i.data in 1:6){
  x<-res[[i.data]]
  tmp1<-x[[5]][1,,]
  tmp1[x[[5]][2,,]>0.05]<-NA
  tmp1<-apply(tmp1,1,mean,na.rm=T)
  
  tmp2<-x[[6]][1,,]
  tmp2[x[[6]][2,,]>0.05]<-NA
  tmp2<-apply(tmp2,1,mean,na.rm=T)
  
  ratio<-tmp1/tmp2
  
  R1$ratio[(a+1):(a+n.spec[i.data])]<-ratio
  R1$data[(a+1):(a+n.spec[i.data])]<-names(res)[i.data]
  a<-a+n.spec[i.data]
}

#### test ##########
i.dataset<-5
rr<-matrix(NA,length(D[[i.dataset]]),6)
colnames(rr)<-c("syn","skew","temp","spat","temp.surr","spat.surr")
for(s in 1:length(D[[i.dataset]])){
  x<-D[[i.dataset]][[s]]
  x<-t(x)
  rr[s,1]<-mean(cor(x,use = "pairwise.complete.obs"), na.rm=T)
  
  tmp<-tls(t(x))
  rr[s,3]<-tmp[1]
  tmp<-tls(x)
  rr[s,4]<-tmp[1]
  
  tmp<-tls.surrogate.all(x,100)
  rr[s,5]<-mean(tmp$temporal[1,])
  rr[s,6]<-mean(tmp$spatial[1,])
  
  x<-as.vector(x)
  x<-x[!is.na(x)]
  rr[s,2]<-as.numeric(skew(x))/sd(x)*mean(x)
  
  show(s)
}
rr<-as.data.frame(rr)
depth<-as.numeric(names(D[[i.dataset]]))
rr$depth<-depth
rr5<-rr
