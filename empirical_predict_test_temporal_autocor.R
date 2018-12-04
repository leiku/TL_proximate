
res<-readRDS("Res_empirical_TL_slopes.RDS")

X.spat<-vector("list", 6)
names(X.spat)<-c("aphid","north.sea","calcofi.gr1","calcofi.gr2","calcofi.gr3","calcofi.gr4")
X.temp<-X.spat
Y.spat<-X.spat
Y.temp<-X.spat
for(i in 1:6){
  res1<-res[[i]]
  x.spat<-rep(NA, length(res1))
  names(x.spat)<-names(res1)
  x.temp<-x.spat
  y.spat<-x.spat
  y.temp<-x.spat
  for(j in 1:length(res1)){
    res2<-res1[[j]]
    x.spat[j]<-res2[1,1]
    x.temp[j]<-res2[1,2]
    
    y.spat[j]<-1.1244*res2[8,1]+0.0124/res2[7,1]+2.3419*res2[7,2]-0.6067
    y.temp[j]<-1.6049*res2[8,1]+0.4332*res2[7,1]+0.0064/res2[7,2]-0.9341
  }
  X.spat[[i]]<-x.spat
  X.temp[[i]]<-x.temp
  Y.spat[[i]]<-y.spat
  Y.temp[[i]]<-y.temp
}


#get correlation

ans<-matrix(NA,7,6)
row.names(ans)<-c("all",names(res))
colnames(ans)<-c("r.spat","p.spat","r.temp","p.temp","eoor.spat","error.temp")
for(i in 1:7){
  if(i==1){
    x<-unlist(X.spat)
    y<-unlist(Y.spat)
    z<-cor.test(x,y)
    ans[i,1]<-z$estimate
    ans[i,2]<-z$p.value
    ans[i,5]<-mean(abs(x-y))
    
    x<-unlist(X.temp)
    y<-unlist(Y.temp)
    z<-cor.test(x,y)
    ans[i,3]<-z$estimate
    ans[i,4]<-z$p.value
    ans[i,6]<-mean(abs(x-y))
  }else{
    x<-X.spat[[i-1]]
    y<-Y.spat[[i-1]]
    z<-cor.test(x,y)
    ans[i,1]<-z$estimate
    ans[i,2]<-z$p.value
    
    x<-X.temp[[i-1]]
    y<-Y.temp[[i-1]]
    z<-cor.test(x,y)
    ans[i,3]<-z$estimate
    ans[i,4]<-z$p.value
    ans[i,5]<-mean(abs(X.spat[[i-1]]-Y.spat[[i-1]]))
    ans[i,6]<-mean(abs(X.temp[[i-1]]-Y.temp[[i-1]]))
  }
}

write.csv(ans,"res_empirical_predict_test.csv",row.names = T)


# plot
cols<-c("red","gold3","green4","cyan3","blue","purple")
tiff(paste0("Fig_empirical_predict_test.tif"),width=7, height=3.5, units="in",compression="zip",res=300)
op<-par(mfrow=c(1,2),oma=c(0.5,0.5,0.5,0.5), mar=c(3,3,0,0),mgp=c(1.8,0.6,0),pty="s")
x<-unlist(Y.spat)
y<-unlist(X.spat)
tmp<-rbind(range(x),range(y))
tmp<-range(tmp)
plot(NA,xlim=tmp, ylim=tmp,
     xlab="predicted spatial TL slopes",ylab="empirical spatial TL slopes")
lines(tmp, tmp,lty='dashed')
for(i in 1:6){
  points(Y.spat[[i]],X.spat[[i]],col=cols[i])
}
mtext(paste0("(",letters[1],")"), side=3, line=-1.5, adj=0.05, cex=1)
legend("bottomright",c("Aphid","Plankton","Chla1","Chla2","Chla3","Chla4"),col=cols,pch=1,cex=0.9)


x<-unlist(Y.temp)
y<-unlist(X.temp)
tmp<-rbind(range(x),range(y))
tmp<-range(tmp)
plot(NA,xlim=tmp, ylim=tmp,
     xlab="predicted temporal TL slopes",ylab="empirical temporal TL slopes")
lines(tmp, tmp, lty='dashed')
for(i in 1:6){
  points(Y.temp[[i]], X.temp[[i]],col=cols[i])  
}
mtext(paste0("(",letters[2],")"), side=3, line=-1.5, adj=0.05, cex=1)
par(op)
dev.off()





