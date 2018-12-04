############### Figs: all distribution ################

library(asbio)
D<-readRDS("data_original.RDS")
n.spec<-c(20,22,rep(10,4))

titles<-c("a_aphid","b_northsea","c_calcofi_gr1","d_calcofi_gr2",
          "e_calcofi_gr3","f_calcofi_gr4")
mfrows<-rbind(c(5,4),c(6,4),c(3,4),c(3,4),c(3,4),c(3,4))
sizes<-rbind(c(8,11),c(8,13),c(8,6),c(8,6),c(8,6),c(8,6))
xlabs<-c("number of aphids", "counts of plankton taxa", rep(expression(paste("concentration (",mu,"g/L)")),4))

for(i.dataset in 1:6){
  
  tiff(paste0("FigS_distribution_",titles[i.dataset],".tif"),
       width=sizes[i.dataset,1], height=sizes[i.dataset,2], units="in",compression="zip",res=300)
  op<-par(mfrow=mfrows[i.dataset,],oma=c(3,3,0,1), mar=c(3,3,1,0),mgp=c(2,1,0))
  for(s in 1:length(D[[i.dataset]])){
    x<-D[[i.dataset]][[s]]
    x<-as.vector(x)
    x<-x[!is.na(x)]
    hist(x,40,ann=F)
    
    CV<-sd(x)/mean(x)
    mtext(paste0("(",letters[s],")"),side=3,line=-1.5,adj=0.95,cex=0.8,font=2)
    mtext(paste0("skewness=",round(skew(x),2)),side=3,line=-3,adj=0.95,cex=0.8,font=2)
    mtext(paste0("CV=",round(CV,2)),side=3,line=-5,adj=0.95,cex=0.8,font=2)
  }
  # xlab + ylab
  par(fig = c(0, 1, 0, 1), oma=c(0,0,0,0), mar = c(5, 5, 0, 0), new = TRUE)
  plot(NA, xlim=c(0,1),ylim=c(0,1), xlab=xlabs[i.dataset], ylab="count",
       type = "n", bty = "n", xaxt = "n", yaxt = "n", cex.axis=1.5, cex.lab=1.5, font.lab=2)
  par(op)
  
  dev.off()
}


