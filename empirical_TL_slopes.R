########## get TL slopes ##########
library(asbio)

D<-readRDS("data_original.RDS")

n.spec<-c(20,22,rep(10,4))
N<-sum(n.spec)
names.data<-c("aphid","north.sea","calcofi.gr1","calcofi.gr2","calcofi.gr3","calcofi.gr4")

res<-vector("list", 6) 
names(res)<-names.data
for(i in 1:6){
  d<-D[[i]]
  res1<-vector("list", length(d))
  names(res1)<-names(d)
  for(j in 1:length(d)){
    d1<-t(d[[j]])
    z.spat<-tls(d1)
    z.temp<-tls(t(d1))
    tmp<-cbind(z.spat, z.temp)
    
    #synchrony
    omega.spat<-mean(cor(d1,use = "pairwise.complete.obs"),na.rm=T)
    omega.temp<-mean(cor(t(d1),use = "pairwise.complete.obs"),na.rm=T)
    
    #J
    J.tmp<-rep(NA, ncol(d1))
    for(jj in 1:ncol(d1)){
      y1<-d1[,jj]
      y1<-y1[!is.na(y1)]
      J.tmp[jj]<-skew(y1)/sd(y1)*mean(y1)
    }
    J.each<-mean(J.tmp)
    d1<-as.vector(d1)
    d1<-d1[!is.na(d1)]
    J.all<-skew(d1)*mean(d1)/sd(d1)
    
    tmp1<-matrix(c(omega.spat, omega.temp,rep(J.each,2), rep(J.all,2)), 3, 2, byrow = T)
    rownames(tmp1)<-c("omega","J.each", "J.all")
    
    tmp<-rbind(tmp, tmp1 )
    res1[[j]]<-tmp
  }
  res[[i]]<-res1
}
saveRDS(res, "Res_empirical_TL_slopes.RDS")


