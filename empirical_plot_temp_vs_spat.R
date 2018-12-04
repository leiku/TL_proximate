#generate a figure 

library(lmodel2)

res<-readRDS("Res_empirical_TL_slopes_good.RDS")
n.spec<-c(19,22, 9, 10, 5, 9)

axis.min<-c(1.8, 1.1, 1.5, 1.4, -0.5, -0.5)
axis.max<-c(3, 2.7, 2.9, 6.5, 3.2, 4.1)
tiff(paste0("Fig2_empirical_temp_vs_spat.tif"),width=3.22, height=4.5, units="in",compression="zip",res=300)
op<-par(mfrow=c(3,2),oma=c(3,3,0,1), mar=c(1.5,1.5,1,0),mgp=c(1,0.6,0),pty="s")

for(i.data in 1:6){
  x<-res[[i.data]]
  spat.slope<-rep(NA, length(x))
  temp.slope<-spat.slope
  for(j in 1:length(x)){
    spat.slope[j]<-x[[j]][1,1]
    temp.slope[j]<-x[[j]][1,2]
  }
  
  y<-data.frame(spat.slope=spat.slope, temp.slope=temp.slope)
  
  z<-lmodel2(temp.slope~spat.slope, data=y, nperm=99)
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


