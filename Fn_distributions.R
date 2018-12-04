#***
#utility functions
#***

#A more convenient negative binomial - the one implemented
#in R is one of the alternative formulations from what
#I used in my math docs, and this does the conversion.
#
#Args
#We use the pmf phi(k)=C(k+r-1,k)(1-p)^rp^k, where C()
#is the binomial coefficient. This is the number of successes
#in a sequence of iid Bernoilli trials before r failures occur,
#where p is the probability of a success.
#
rnbinom_my<-function(n,r,p)
{
  if (r==0)
  { #A special case, for convenience
    return(numeric(n))
  }
  return(rnbinom(n,r,1-p))
}

dnbinom_my<-function(x,r,p,log=F)
{
  return(dnbinom(x,size=r,prob=1-p,log=log))
}

qnbinom_my<-function(x,r,p,log=F)
{
  return(qnbinom(x,r,p,log.p=log))  
}

#***
#models
#***

#***The below allow non-identically distributed versions of cases
#we have considered previously for identically distributed versions

#A multivariate random variable with gamma marginals, but with
#marginals not necessarily identically distributed.
#
#Args
#E        An array of standard normal draws, n=dim(E)[2]
#params   An n by 3 matrix, with columns giving alpha and beta
#           (shape and rate) for the gamma marginals
#
model_NonIDGamma<-function(E,params)
{
  alphas<-params[,1]
  betas<-params[,2]
  mult<-params[,3]
  if (length(alphas)!=dim(E)[2]) { stop("Error in model_NonIDGamma") }
  Y<-qgamma(pnorm(E),
            shape=array(rep(alphas,each=dim(E)[1]),dim=dim(E)),
            rate=array(rep(betas,each=dim(E)[1]),dim=dim(E)))*
    array(rep(mult,each=dim(E)[1]),dim=dim(E))
  return(Y)
}

#A multivariate random variable with Poisson marginals, but with
#marginals not necessarily identically distributed.
#
#Args
#E        An array of standard normal draws, n=dim(E)[2]
#params   An n by 2 matrix, lambda values for the Poisson
#           marginals you want, and the mulitpliers
#
model_NonIDPois<-function(E,params)
{
  lambdas<-params[,1]
  mult<-params[,2]
  if (length(lambdas)!=dim(E)[2]) { stop("Error in model_NonIDPoiss") }
  Y<-qpois(pnorm(E),
           lambda=array(rep(lambdas,each=dim(E)[1]),dim=dim(E)))*
    array(rep(mult,each=dim(E)[1]),dim=dim(E))
  return(Y)
}

#A multivariate random variable with negative-binomial marginals, 
#but with marginals not necessarily identically distributed.
#
#Args
#E        An array of standard normal draws, n=dim(E)[2]
#params   An n by 3 matrix, with columns giving r and p
#
model_NonIDNegBinom<-function(E,params)
{
  rs<-params[,1]
  ps<-params[,2]
  mult<-params[,3]
  if (length(rs)!=dim(E)[2]) { stop("Error in model_NonIDNegBinom") }
  Y<-qnbinom_my(pnorm(E),
                r=array(rep(rs,each=dim(E)[1]),dim=dim(E)),
                p=array(rep(ps,each=dim(E)[1]),dim=dim(E)))*
    array(rep(mult,each=dim(E)[1]),dim=dim(E))
  return(Y)
}

#A multivariate random variable with normal marginals, but with
#marginals not necessarily identically distributed.
#
#Args
#E        An array of standard normal draws, n=dim(E)[2]
#params   An n by 3 matrix, with columns giving mu and sigma
#           for the marginals
#
model_NonIDNorm<-function(E,params)
{
  mus<-params[,1]
  sigmas<-params[,2]
  mult<-params[,3]
  if (length(mus)!=dim(E)[2]) { stop("Error in model_NonIDNorm") }
  Y<-(E*array(rep(sigmas,each=dim(E)[1]),dim=dim(E))+
        array(rep(mus,each=dim(E)[1]),dim=dim(E)))*
    array(rep(mult,each=dim(E)[1]),dim=dim(E))
  return(Y)
}

#A multivariate random variable with lognormal marginals, but with
#marginals not necessarily identically distributed.
#
#Args
#E        An array of standard normal draws, n=dim(E)[2]
#params   An n by 3 matrix, with columns giving mu and sigma
#           for the lognormal marginals
#
model_NonIDLognorm<-function(E,params)
{
  mus<-params[,1]
  sigmas<-params[,2]
  mult<-params[,3]
  if (length(mus)!=dim(E)[2]) { stop("Error in model_NonIDLognormal") }
  Y<-qlnorm(pnorm(E),
            meanlog=array(rep(mus,each=dim(E)[1]),dim=dim(E)),
            sdlog=array(rep(sigmas,each=dim(E)[1]),dim=dim(E)))*
    array(rep(mult,each=dim(E)[1]),dim=dim(E))
  return(Y)
}




