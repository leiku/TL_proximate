#Computes Taylor's law slope (spatial if each column is a timeseries) and related information.
#
#Args
#m          A matrix (get spatial TL if each column is a timeseries)
#
#Output
#A vector with these elements:

#The Taylor's law regression slope
#The Taylor's law regression intercept 
#The p-value for a quadratic comparison
#The p-value for a homoskedasticity test
#The root mean squared resdiuals from the Taylor regression 
#The r-value from the Taylor plot


#Means and variances are computed across the rows
#
tls<-function(m)
{
  #Prepare the points for the Taylor's law plot
  ms<-apply(X=m,MARGIN=1,FUN=mean,na.rm=T)
  vs<-apply(X=m,MARGIN=1,FUN=var,na.rm=T)
  
  tmp<-which(is.finite(ms) & is.finite(vs))
  ms<-ms[tmp]
  vs<-vs[tmp]
  
  inds<-(ms>0 & vs>0)
  lms<-log10(ms[inds])
  lvs<-log10(vs[inds])
  lms2<-lms^2
  
  #Do various regressions and tests
  mod1<-lm(lvs~lms)
  mod2<-lm(lvs~lms2+lms)
  hq<-anova(mod1,mod2) #make the quadratic comparison
  absresids<-abs(resid(mod1))
  predicts<-predict(mod1)
  mod3<-lm(absresids~predicts)
  hhet<-anova(mod3)
  corr_r<-cor(lvs,lms)
  
  ans<-c(unname(mod1$coefficients[2]),unname(mod1$coefficients[1]),
         unname(hq[6][2,1]),unname(hhet[5][1,1]), sqrt(mean(resid(mod1)^2)), corr_r, summary(mod1)$coefficients[2,4])
  names(ans)<-c("slope","intercept","p.quad","p.het","rms","r","p.fit")
  return(ans)
}
