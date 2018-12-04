#input: res: the RDS file got from Fn_averaged_igood
#       i.good.comb: the RDS file got from Fn_averaged_igood
#output: a tiff figure showing barplot of TL slopes as a function of rhos

Fn_plot_barplot<-function(res, i.good.comb, nN, file.name){
  
  
  name.model<-c("Poisson","Negative binomial","Exponential","Gamma","Lognormal")
  rhos<-unique(res$rho)
  dim.ratios<-unique(res$dim.ratio)
  
  tiff(paste0(file.name, ".tif"),width=6.7, height=7.5, units="in",compression="zip",res=600)
  op<-par(mfcol=c(5,2),oma=c(3,3,1,1), mar=c(3,2,0.5,0.5),mgp=c(1,0.5,0))
  
  i.test<-c(12,14,13,15)
  cols<-c("red","purple","gold3","darkgreen","blue")
  slope.predicted<-c("Poisson (1.0)","Negative binomial (1.6)","Exponential (2.0)","Gamma (2.0)","Lognormal (4.7)")
  ylims<-c(3.6,5.2)
  for(i in 1:2){     # different test
    for(j in 1:5){   # different distribution
      d<-subset(res, distribution==name.model[j])
      
      i.good<-i.good.comb[j,,]
      y<-t(matrix(d[,i.test[i]],length(rhos),length(dim.ratios)))
      z<-t(matrix(d[,i.test[i+2]],length(rhos),length(dim.ratios)))
      #z<-z/sqrt(n.runs)
      
      y<-y*i.good
      z<-z*i.good
      if(j==1 & i==2){
        x<-barplot(y,beside=T,col=cols,ylim=c(0,ylims[i]),cex.axis = 1.1)
        arrows(x,y+z,x,y,length=0.02, angle=90, code=3,col="black")
        legend("topright",legend=nN, fill=cols, horiz=T)
      }else{x<-barplot(y,beside=T,col=cols,ylim=c(0,ylims[i]),cex.axis = 1.1)}
      arrows(x,y+z,x,y,length=0.02, angle=90, code=3,col="black")
      #box()
      
      mtext(paste0("(",letters[2*j-2+i],")"), side=3, line=0.1,adj=0.0,cex=0.8)
      if(i==1){mtext(slope.predicted[j], side=3, line=-1, adj=0.5, cex=0.8)}
    }
  }
  # xlab + ylab
  par(fig = c(0, 1, 0, 1), oma=c(0,0,0,0), mar = c(0, 0, 0, 0), mgp=c(1,0,0), new = TRUE)
  plot(NA, ann=F,xlim=c(0,1),ylim=c(0,1), bty = "n", xaxt = "n", yaxt = "n")
  mtext(expression(rho), side=1, line=-2, adj=0.5, cex=1)
  mtext("TL slopes", side=2, line=-2, adj=0.5, cex=1)
  par(op)
  dev.off()
  
}

