#plot a figure showing temporal vs. spatial TL slopes

Fn_simulation_temp_vs_spat<-function(res, i.good.comb){
  name.model<-unique(res$distribution)
  dim.ratios<-unique(res$dim.ratio)
  rhos<-unique(res$rho)
  
  cols<-c("red","purple","gold3","darkgreen","blue")
  
  tiff("Fig_simulation_temporal_vs._spatial_average.tif",width=5, height=5, units="in",compression="zip",res=300)
  x.range<-range(res$ave.slope.spat)
  y.range<-range(res$ave.slope.temp)
  plot(NA,xlim=x.range,ylim=y.range,xlab="spatial TL slope", ylab="temporal TL slope")
  xx<-min(c(x.range, y.range))
  yy<-max(c(x.range, y.range))
  lines(c(xx,yy), c(xx,yy), col="black", lwd=2)
  
  for(ii in 1:5){ #different distributions
    res1<-subset(res, distribution==name.model[ii])
    
    x<-0;
    y<-0;
    for(jj in 1:5){ # different DR
      res2<-subset(res1, dim.ratio==dim.ratios[jj])
      for(i in 1:10){ # different omega
        res3<-subset(res2, rho==rhos[i])
        if(!is.na(i.good.comb[ii,jj,i])){
          x<-c(x, res3$ave.slope.spat)
          y<-c(y, res3$ave.slope.temp)
        }
      }
    }
    x<-x[-1]
    y<-y[-1]
    
    points(y~x,col=cols[ii])
  }
  legend("topleft",  name.model, col=cols, pch=1)
  dev.off()
  
}

