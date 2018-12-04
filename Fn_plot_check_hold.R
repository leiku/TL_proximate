#input: res: the RDS file got from Fn_averaged_igood
#       i.good.comb: the RDS file got from Fn_averaged_igood
#output: a tiff figure showing fraction of TL holdings

Fn_plot_check_hold<-function(res,  file.name){
  name.model<-c("Poisson","Negative binomial","Exponential","Gamma","Lognormal")
  dim.ratios<-unique(res$dim.ratio)
  
  tiff(paste0(file.name,".tif"),width=8, height=10, units="in",compression="zip",res=600)
  op<-par(mfrow=c(6,5),oma=c(4,4,3,1), mar=c(2,2,0.5,0.5),mgp=c(1.5,0.8,0),pty="s")
  
  name.test<-c("spatial linearity","temporal linearity","spatial heteroscedasticity","temporal heteroscedasticity",
               "spatial RMSE", "temporal RMSE")
  i.test<-c(4:7,8,10,9,11)
  cols<-c("red","purple","gold3","darkgreen","blue")
  a<-1
  for(i in 1:6){     # different test
    for(j in 1:5){     # different distribution
      d<-subset(res, distribution==name.model[j])
      
      if(i<5){
        plot(NA,xlim=c(0,0.9),ylim=c(0,1),ann=F)
        lines(c(-1,2),c(0.1,0.1), lty="dashed",col="black")
        for(jj in 1:5){
          d1<-subset(d, dim.ratio==dim.ratios[jj])    # different dim ratio
          lines(d1$rho,d1[,i.test[i]],col=cols[jj])
        }
      }else{
        plot(NA,xlim=c(0,0.9),ylim=c(0,1),ann=F)
        for(jj in 1:5){
          d1<-subset(d, dim.ratio==dim.ratios[jj])    # different dim ratio
          lines(d1$rho,d1[,i.test[i]],col=cols[jj])
          points(d1$rho,d1[,i.test[i]],col=cols[jj],pch=20)
          arrows(d1$rho, d1[,i.test[i]]-d1[,i.test[i+2]], d1$rho, d1[,i.test[i]]+d1[,i.test[i+2]],
                 length=0.05, angle=90, code=3,col=cols[jj])
        }
      }
      mtext(paste0("(",a,")"), side=3, line=-1.2,adj=0.05,cex=0.8)
      a<-a+1
    }
  }
  # xlab + ylab
  par(fig = c(0, 1, 0, 1), oma=c(0,0,0,0), mar = c(0, 0, 0, 0), mgp=c(1,0,0), new = TRUE)
  plot(NA, ann=F,xlim=c(0,1),ylim=c(0,1), bty = "n", xaxt = "n", yaxt = "n")
  mtext(expression(rho), side=1, line=5, adj=0.5, cex=1)
  mtext("fraction of p < 0.01", side=2, line=-3, adj=0.76, cex=1)
  mtext("root mean squared residuals", side=2, line=-3, adj=0.01, cex=1)
  
  par(op)
  dev.off()
}


