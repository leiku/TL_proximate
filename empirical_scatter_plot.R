
#################### FigS  spatial TL scatter plots ######################
res<-readRDS("Res_empirical_TL_slopes.RDS")
D<-readRDS("data_original.RDS")

titles<-c("a_aphid","b_northsea","c_calcofi_gr1","d_calcofi_gr2",
          "e_calcofi_gr3","f_calcofi_gr4")
mfrows<-rbind(c(5,4),c(6,4),c(3,4),c(3,4),c(3,4),c(3,4))
sizes<-rbind(c(8,11),c(8,13),c(8,6),c(8,6),c(8,6),c(8,6))

for(i.dataset in 1:6){
  tiff(paste0("FigS_empirical_scatter_spatial_",titles[i.dataset],".tif"),
       width=sizes[i.dataset,1], height=sizes[i.dataset,2], units="in",compression="zip",res=300)
  op<-par(mfrow=mfrows[i.dataset,],oma=c(3,3,0,1), mar=c(3,3,1,0),mgp=c(2,1,0))
  for(s in 1:length(res[[i.dataset]])){
    d<-D[[i.dataset]][[s]]
    d<-t(d)
    
    ms<-apply(d, 1, mean, na.rm=T)
    vs<-apply(d, 1, var, na.rm=T)
    lms<-log10(ms)
    lvs<-log10(vs)
    
    tmp<-which(is.finite(lms) & is.finite(lvs))
    lms<-lms[tmp]
    lvs<-lvs[tmp]
    
    z<-lm(lvs~lms)
    z.fitted<-fitted(z)
    tmp<-c(which.min(lms),which.max(lms))
    
    plot(NA,xlim=range(lms),ylim=c(min(lvs), max(lvs)+0.2),ann=F)
    points(lms,lvs,cex=0.9)
    if(summary(z)$coefficient[2,4]>0.05){
      lines(lms[tmp],z.fitted[tmp],lty="dashed")
    }else{lines(lms[tmp],z.fitted[tmp],lty="solid")}
    
    res1<-res[[i.dataset]][[s]]
    mtext(paste0("(",letters[s],")"),side=3,line=-1.5,adj=0.05,cex=0.8,font=2)
    mtext(bquote(p[q]==.(round(res1[3,1],3))), side=3, line=-3, adj=0.05, cex=0.8)
    mtext(bquote(p[het]==.(round(res1[4,1],3))), side=3, line=-4.5, adj=0.05, cex=0.8)
    mtext(bquote(slope==.(round(res1[1,1],2))), side=1, line=-3, adj=0.95, cex=0.8)
    mtext(bquote(r^2==.(round(summary(z)$r.squared,2))), side=1, line=-1.5, adj=0.95, cex=0.8)
  }
  # xlab + ylab
  par(fig = c(0, 1, 0, 1), oma=c(0,0,0,0), mar = c(5, 5, 0, 0), new = TRUE)
  plot(NA, xlim=c(0,1),ylim=c(0,1), xlab=bquote(log[10]~"mean"),ylab=bquote(log[10]~"variance"),
       type = "n", bty = "n", xaxt = "n", yaxt = "n", cex.axis=1.5, cex.lab=1.5, font.lab=2)
  par(op)
  
  dev.off()
}


#################### FigS temporal TL scatter plots ######################

for(i.dataset in 1:6){
  tiff(paste0("FigS_empirical_scatter_temporal_",titles[i.dataset],".tif"),
       width=sizes[i.dataset,1], height=sizes[i.dataset,2], units="in",compression="zip",res=300)
  op<-par(mfrow=mfrows[i.dataset,],oma=c(3,3,0,1), mar=c(3,3,1,0),mgp=c(2,1,0))
  for(s in 1:length(res[[i.dataset]])){
    d<-D[[i.dataset]][[s]]
    d<-t(d)
    
    ms<-apply(d, 2, mean, na.rm=T)
    vs<-apply(d, 2, var, na.rm=T)
    lms<-log10(ms)
    lvs<-log10(vs)
    
    tmp<-which(is.finite(lms) & is.finite(lvs))
    lms<-lms[tmp]
    lvs<-lvs[tmp]
    
    z<-lm(lvs~lms)
    z.fitted<-fitted(z)
    tmp<-c(which.min(lms),which.max(lms))
    
    plot(NA,xlim=range(lms),ylim=c(min(lvs), max(lvs)+0.2),ann=F)
    points(lms,lvs,cex=0.9)
    if(summary(z)$coefficient[2,4]>0.05){
      lines(lms[tmp],z.fitted[tmp],lty="dashed")
    }else{lines(lms[tmp],z.fitted[tmp],lty="solid")}
    
    res1<-res[[i.dataset]][[s]]
    mtext(paste0("(",letters[s],")"),side=3,line=-1.5,adj=0.05,cex=0.8,font=2)
    mtext(bquote(p[q]==.(round(res1[3,2],3))), side=3, line=-3, adj=0.05, cex=0.8)
    mtext(bquote(p[het]==.(round(res1[4,2],3))), side=3, line=-4.5, adj=0.05, cex=0.8)
    mtext(bquote(slope==.(round(res1[1,2],2))), side=1, line=-3, adj=0.95, cex=0.8)
    mtext(bquote(r^2==.(round(summary(z)$r.squared,2))), side=1, line=-1.5, adj=0.95, cex=0.8)
  }
  # xlab + ylab
  par(fig = c(0, 1, 0, 1), oma=c(0,0,0,0), mar = c(5, 5, 0, 0), new = TRUE)
  plot(NA, xlim=c(0,1),ylim=c(0,1), xlab=bquote(log[10]~"mean"),ylab=bquote(log[10]~"variance"),
       type = "n", bty = "n", xaxt = "n", yaxt = "n", cex.axis=1.5, cex.lab=1.5, font.lab=2)
  par(op)
  
  dev.off()
}

