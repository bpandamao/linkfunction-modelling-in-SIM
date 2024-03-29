library(jsonlite)

####
#function
####

library(mgcv);library(mvtnorm);library(truncnorm);library(GIGrvg);library(invgamma);library(bsplinePsd);library(beyondWhittle)


gibbs_bspline_est_linkfunction= function(X,Y,
                                        S=40000,
                                        L=20,
                                        Kmax=500,
                                        M=1,
                                        thetak=0.01,
                                        K=20,
                                        d=3,
                                        nu0=0.01,
                                        s0=0.01,
                                        A=5,
                                        d_a=0.5,
                                        thin=10,
                                        sigma=0.5,
                                        alpha1=0.01,
                                        alpha2=0.01,
                                        burnin=20000,
                                        baseb=exp(1),
                                        logstate=TRUE,
                                        YM){

 
  
  n=nrow(X)
  p=ncol(X)
  ee=1e-20 
  
  ###
  # The initialization
  ###
  
  # V Z for the weights
  V= beyondWhittle:::vFromP(rep(1 / (L + 1), L))
  Z= seq(from=1 / (2 * K), to = 1 - 1 / (2 * K), length.out = L + 1)
  U= beyondWhittle:::vFromP(rep(1 / (L + 1), L))
  C= seq(from = 1 / (2 * K), to = 1 - 1 / (2 * K), length.out = L + 1)
  
  P = bsplinePsd:::pFromV(V)
  G = bsplinePsd:::mixtureWeight(P, Z, K)
  
  
  # U C for the knots
  Q = bsplinePsd:::pFromV(U)  # u here
  newK = K - d  # CAUTION: This is the denominator for the second DP. 
  knot.diffs = bsplinePsd:::mixtureWeight(Q, C, newK)
  knots = c(0, cumsum(knot.diffs))
  
  # B-spline density matrix
  #db.list <- dbspline(w, knots, degree) 
  
  # Initialise r
  r=rbinom(p,1,1/2)
  if(all(r==0)){r[sample(1:p,1)]=1}
  
  #Initialise beta
  mcmc_r_beta_ratio=rexp(1,3)
  mcmc_beta_r_rho=0.2
  beta=rep(0,p)
  beta[which(r==1)]=rnorm(length(which(r==1)),0,1) #unit norm
  beta=beta/c(sqrt(crossprod(beta)));
  if(beta[which(r==1)[1]]<0){beta=-beta} #beta1>0
  
  #Initialise kc
  kc=runif(n,0,3) 
  
  #if log
  if(logstate){
    mu_f=function(x,baseb){
      return(logb(x,baseb))
    } 
  }else{
    mu_f=function(x,baseb){
      return(x)
    }
  }
  
  ###
  # Parameter space
  ###
  KK=AA=Sigma=sr=sr.star=Raa=NULL
  R=Beta=matrix(nr=S,nc=p)
# just use GG and knots in out
#   VV=matrix(nr=S,nc=L)
#   ZZ=matrix(nr=S,nc=L+1)
#   UU=matrix(nr=S,nc=L)
#   CC=matrix(nr=S,nc=L+1)
#   KC=matrix(nr=S,nc=n)
  GG=matrix(nr=S,nc=Kmax)
  GG_knots=matrix(nr=S,nc=Kmax)
  Rho=R_beta_ratio=NULL
  
  
  ptime = proc.time()[1]
  
  ###
  # target parameters updating process
  ###
  for(s in 1:S){
    
    if (s %% 1000 == 0) {
      print(paste("Iteration", s, ",", "Time elapsed", 
                  round(as.numeric(proc.time()[1] - ptime) / 60, 2),
                  "minutes"))
    }
    
    ### step1: K by Metropolis-within-Gibbs
    
    # use random walk for K
    ss = runif(1)
    if (ss < 0.75) {  
      jump = sample(-1:1, 1, prob = rep(1 / 3, 3))
    }else 
    {
      jump = round(rt(1, 1))  #discrete Cauchy
    }
    K.star = K + jump
    while (K.star < (d+2) || K.star > Kmax) {  # A bit hacky to ensure k doesn't go out of bounds
      if (ss < 0.75) {
        jump = sample(-1:1, 1, prob = rep(1 / 3, 3))  # Ordinary proposal
      }
      else {
        jump = round(rt(1, 1))  # Bold proposal
      }
      K.star = K + jump
    }
    
    if(K.star!=K){# the weights will change
      G.star= bsplinePsd:::mixtureWeight(P, Z, K.star)
      
      ### knots.star
      newK.star = K.star - d  # CAUTION: This is the denominator for the second DP. 
      knot.diffs.star = bsplinePsd:::mixtureWeight(Q, C, newK.star)
      knots.star = c(0, cumsum(knot.diffs.star))
      
      ### calculate the log posterior 
      XB=X%*%beta
      s1=s2=0
      w=exp(XB)/(1+exp(XB))
      dblist = dbspline(w, knots, d) 
      dblist.star = dbspline(w, knots.star, d)
      
      sr <- bsplinePsd:::densityMixture(G, dblist)
      sr <- pmax(sr, ee)
      
      sr.star <- bsplinePsd:::densityMixture(G.star, dblist.star)
      sr.star <- pmax(sr.star, ee)
      
      s1=sum(dnorm(Y,mu_f(A*sr.star,baseb),sqrt(8*kc*sigma),log=TRUE))
      s2=sum(dnorm(Y,mu_f(A*sr,baseb),sqrt(8*kc*sigma),log=TRUE))
      
      ### the accept or reject for K.star
      lnr=s1-s2-thetak*(K.star^2-K^2)
      if(min(exp(lnr),1)>runif(1)){
        # store the outcome
        K=K.star
        G=G.star
        knots=knots.star
        sr=sr.star
        newK=newK.star
        dblist=dblist.star
      }
    } # K.star!=K
    ### outcome store
    KK[s]=K 
    
    
    
    
    ### step3: r by Metropolis-within-Gibbs (q is sampled directly)
    
    # update rj in random order
    for(j in sample(1:p,p)){
      
      # sample q directly
      q=rbeta(1,1+r[j],2-r[j]) 
      rj.star=rbinom(1,1,q) # the proposal
      r.star=NULL
      r.star=r
      r.star[j]=rj.star
      
      ### if Bspline mixture density is NULL
      if(length(sr)==0){ 
        XB=X%*%beta
        w=exp(XB)/(1+exp(XB))
        dblist = dbspline(w, knots, d) 
        sr <- bsplinePsd:::densityMixture(G, dblist)
        sr <- pmax(sr, ee)
      }
      
      if(rj.star!=r[j] & sum(r.star)>=1){
        # updating q
        q.star=rbeta(1,1+rj.star,2-rj.star)
        # updating beta accordingly
        beta.star=NULL
        if(sum(r.star)==1){# one dimension
          beta.star=rep(0,p)
          beta.star[which(r.star==1)]=1
        }else
        {
          beta.star=beta
          non_zeron=which(Beta[1:(s-1),j]!=0)
          if(length(non_zeron)<=1){
            beta.star[j]=rj.star*rnorm(1,0,1)
          }else
          {
            beta.star[j]=rj.star*rnorm(1,mean(Beta[non_zeron,j]),sd(Beta[non_zeron,j]))
            #beta.star[j]=rj.star*rnorm(1,mean(Beta[non_zeron,j]),1)
          }
          # the normalisation
          beta.star=beta.star/c(sqrt(crossprod(beta.star)))
          if(beta.star[which(r.star==1)[1]]<0){beta.star=-beta.star}
        }# sum(r.star)>1
        nr=sum(r)
        nr.star=sum(r.star)
        
        ### calculate the log posterior 
        s1=s2=0
        XB.star=X%*%beta.star
        w.star=exp(XB.star)/(1+exp(XB.star))
        dblist.star = dbspline(w.star, knots, d)
        
        sr.star = bsplinePsd:::densityMixture(G, dblist.star)
        sr.star = pmax(sr.star, ee)
        
        s1=sum(dnorm(Y,mu_f(A*sr.star,baseb),sqrt(8*kc*sigma),log=TRUE))
        s2=sum(dnorm(Y,mu_f(A*sr,baseb),sqrt(8*kc*sigma),log=TRUE))
        
        ### the accept or reject for r.star
        lnr=s1-s2+(r[j])*log(q.star/q)+(1-r[j])*log((1-q.star)/(1-q))+(nr-nr.star)*log(pi)/2+log(gamma(nr.star/2))-log(gamma(nr/2))
        if(min(exp(lnr),1)>runif(1)){
          r[j]=rj.star
          beta=beta.star
          sr=sr.star
          dblist=dblist.star
        }
      }# rj.star!=r[j] & sum(r.star)>=1
    } # j ends
    ### outcome store
    R[s,]=r
    
    
    ### step4: beta by Metropolis-within-Gibbs (the rho is updated; rho is the multivariate normal distribution of non zero beta)
    
    non_zerob=which(r==1)
    if(length(non_zerob)>1){
      
      # the beta updating
      beta.star=rep(0,p)
      # the rho updating
      mcmc_beta_r_rho=max(0,mcmc_beta_r_rho+(0.234-mcmc_r_beta_ratio)/sqrt(s))
      Rho[s]=mcmc_beta_r_rho
      
      # the multivariate normal distribution
      beta.star[non_zerob]=as.vector(mvtnorm:::rmvnorm(1,sqrt(2)*mcmc_beta_r_rho*beta[non_zerob],diag(length(non_zerob))))
      beta.star=beta.star/sqrt(c(crossprod(beta.star[non_zerob])))
      if(beta.star[which(beta.star!=0)[1]]<0){beta.star=-beta.star}
      
      
      ### calculate the log posterior 
      s1=s2=0
      XB.star=X%*%beta.star
      w.star=exp(XB.star)/(1+exp(XB.star))
      dblist.star = dbspline(w.star, knots, d)
      
      sr.star <- bsplinePsd:::densityMixture(G, dblist.star)
      sr.star <- pmax(sr.star, ee)
      
      s1=sum(dnorm(Y,mu_f(A*sr.star,baseb),sqrt(8*kc*sigma),log=TRUE))
      s2=sum(dnorm(Y,mu_f(A*sr,baseb),sqrt(8*kc*sigma),log=TRUE))
      
      ### the accept or reject for beta.star
      lnr=s1-s2+(t(beta.star-sqrt(2)*mcmc_beta_r_rho*beta)%*%diag(rep(p))%*%(beta.star-sqrt(2)*mcmc_beta_r_rho*beta)-t(beta-sqrt(2)*mcmc_beta_r_rho*beta.star)%*%diag(rep(p))%*%(beta-sqrt(2)*mcmc_beta_r_rho*beta.star))/2
      # updating ratio immediately
      mcmc_r_beta_ratio=min(exp(lnr),1)
      R_beta_ratio[s]=mcmc_r_beta_ratio
      
      if(mcmc_r_beta_ratio>runif(1)){
        beta=beta.star;
        sr=sr.star
        dblist=dblist.star
      }
    }# length(non_zerob)>1
    Beta[s,]=beta
    
    ### the A
    A.s=runif(1,max(0,A-d_a),A+d_a)
    #A.s=rtruncnorm(1,0,Inf,A,2.38/sqrt(p))
    s1=s2=0
    
    s1=sum(dnorm(Y,mu_f(A.s*sr,baseb),sqrt(8*kc*sigma),log=TRUE))
    s2=sum(dnorm(Y,mu_f(A*sr,baseb),sqrt(8*kc*sigma),log=TRUE))
    
    lnr=s1-s2+alpha2*(A-A.s)+(alpha1-1)*(log(A.s)-log(A))
    raa=min(exp(lnr),1);if(raa>runif(1)){A=A.s}
    if(raa<0.15){d_a=d_a/2}else if(raa>0.65){d_a=d_a*2}
    AA[s]=A;Raa[s]=raa
    
    
    
    ### step5: V by Metropolis-within-Gibbs (V[l] l is randomly updated)
    
    for(l in sample(1:L,L)){
      
      # random walk of V[]
      dv=l/(l+2*sqrt(n))
      Vl.star=runif(1,V[l]-dv,V[l]+dv)
      if (Vl.star>1){
        Vl.star=Vl.star-1
      }else if (Vl.star<0){
        Vl.star=Vl.star+1
      }
      V.star=V
      V.star[l]=Vl.star
      
      # G is updated
      P.star = bsplinePsd:::pFromV(V.star)
      G.star = bsplinePsd:::mixtureWeight(P.star, Z, K)
      
      sr.star <- bsplinePsd:::densityMixture(G.star, dblist)
      sr.star <- pmax(sr.star, ee)
      
      # the log posterior 
      s1=s2=0
      s1=sum(dnorm(Y,mu_f(A*sr.star,baseb),sqrt(8*kc*sigma),log=TRUE))
      s2=sum(dnorm(Y,mu_f(A*sr,baseb),sqrt(8*kc*sigma),log=TRUE))
      
      ### the accept or reject for V.star
      lnr=s1-s2+(M-1)*(log(1-Vl.star)-log(1-V[l]))
      if(min(exp(lnr),1)>runif(1)){
        V[l]=Vl.star
        G=G.star
        sr=sr.star
        P=P.star
      }
    }# l ends
    
    ### step6: Z by Metropolis-within-Gibbs (Z[l] l is randomly updated)
    
    for(l in sample(1:(L+1),L+1)){
      
      dz=l/(l+2*sqrt(n))
      Zl.star=runif(1,Z[l]-dz,Z[l]+dz)
      if (Zl.star>1){
        Zl.star=Zl.star-1
      }else if (Zl.star<0){
        Zl.star=Zl.star+1
      }
      Z.star=Z
      Z.star[l]=Zl.star
      
      G.star = bsplinePsd:::mixtureWeight(P, Z.star, K)
      
      sr.star = bsplinePsd:::densityMixture(G.star, dblist)
      sr.star = pmax(sr.star, ee)
      
      # the log posterior 
      s1=s2=0
      s1=sum(dnorm(Y,mu_f(A*sr.star,baseb),sqrt(8*kc*sigma),log=TRUE))
      s2=sum(dnorm(Y,mu_f(A*sr,baseb),sqrt(8*kc*sigma),log=TRUE))
      
      lnr=s1-s2
      if(min(exp(lnr),1)>runif(1)){
        Z=Z.star
        G=G.star
        sr=sr.star
      }
    }# l ends
    # with V,Z ready G.star is stored
    GG[s,1:K]=G
    
    
    ### step7: U by Metropolis-within-Gibbs (The knots updates))
    
    for(l in sample(1:L,L)){
      
      du=l/(l+2*sqrt(n))
      Ul.star=runif(1,U[l]-du,U[l]+du)
      if (Ul.star>1){
        Ul.star=Ul.star-1
      }else if (Ul.star<0){
        Ul.star=Ul.star+1
      }
      U.star=U
      U.star[l]=Ul.star
      
      ### knots.star
      Q.star = bsplinePsd:::pFromV(U.star)
      knot.diffs.star = bsplinePsd:::mixtureWeight(Q.star, C, newK)
      knots.star = c(0, cumsum(knot.diffs.star))
      
      ### calculate the log posterior 
      s1=s2=0
      XB=X%*%beta
      w=exp(XB)/(1+exp(XB))
      dblist.star = dbspline(w, knots.star, d)
      
      sr.star <- bsplinePsd:::densityMixture(G, dblist.star)
      sr.star <- pmax(sr.star, ee)
      
      s1=sum(dnorm(Y,mu_f(A*sr.star,baseb),sqrt(8*kc*sigma),log=TRUE))
      s2=sum(dnorm(Y,mu_f(A*sr,baseb),sqrt(8*kc*sigma),log=TRUE))
      
      lnr=s1-s2+(M-1)*(log(1-Ul.star)-log(1-U[l]))                                                                                                                                                  
      if(min(exp(lnr),1)>runif(1)){
        U=U.star
        Q=Q.star
        knots=knots.star
        dblist=dblist.star
        sr=sr.star
      }
    }## l ends
    
    
    ### step8: C by Metropolis-within-Gibbs (The knots updates))
    
    for(l in sample(1:(L+1),L+1)){
      
      dc=l/(l+2*sqrt(n))
      Cl.star=runif(1,C[l]-dc,C[l]+dc)
      if (Cl.star>1){
        Cl.star=Cl.star-1
      }else if (Cl.star<0){
        Cl.star=Cl.star+1
      }
      C.star=C
      C.star[l]=Cl.star
      
      ### knots.star
      knot.diffs.star = bsplinePsd:::mixtureWeight(Q, C.star, newK)
      knots.star = c(0, cumsum(knot.diffs.star))
      
      ### calculate the log posterior 
      s1=s2=0
      dblist.star = dbspline(w, knots.star, d)
      
      sr.star <- bsplinePsd:::densityMixture(G, dblist.star)
      sr.star <- pmax(sr.star, ee)
      
      s1=sum(dnorm(Y,mu_f(A*sr.star,baseb),sqrt(8*kc*sigma),log=TRUE))
      s2=sum(dnorm(Y,mu_f(A*sr,baseb),sqrt(8*kc*sigma),log=TRUE))
      
      lnr=s1-s2
      if(min(exp(lnr),1)>runif(1)){
        C=C.star
        knots=knots.star
        dblist=dblist.star
        sr=sr.star
      }
    } # l ends
    GG_knots[s,1:(K-d+1)]=knots
    
    
    ### step9: sigma& kc (sampled directly)
    
    mu=t(Y)-mu_f(A*sr,baseb)
    sye=mu%*%diag(1/kc)%*%t(mu)
    sigma=rinvgamma(1,nu0+3*n/2,s0+sum(kc)+sye/16)
    Sigma[s]=sigma
    
    
    eta2=2/sigma
    for (i in 1:n){
      eta1=((mu[i])^2)/(8*sigma)
      kc[i]=rgig(1,1/2,eta1,eta2)
    } 
    

    # END: MCMC loop
  }
  
  ###
  # Which iterations to keep
  ###
  
  keep=seq(burnin+1,S,by=thin)
  KK=KK[keep]
  AA=AA[keep]
  R=R[keep,]
  Beta=Beta[keep,]
  GG=GG[keep,]
  GG_knots=GG_knots[keep,]
  
  
  ###
  # the estimate
  ###
  
  R.pred=Beta.pred=NULL
  # variable selection
  R.pred=apply(R,2,median)
  Beta.pred=apply(Beta,2,median)
  
  # Y.est
  Y.est=matrix(nr=n,nc=length(keep))
  w.pred=exp(X%*%Beta.pred)/(1+exp(X%*%Beta.pred))
                             
  for (isample in 1:length(keep)){
    dblist.est = dbspline(w.pred, GG_knots[isample,1:(KK[isample]-d+1)], d)
    sr.est = bsplinePsd:::densityMixture(GG[isample,1:KK[isample]], dblist.est)
    sr.est = pmax(sr.est, ee)
    Y.est[,isample]= mu_f(AA[isample]*sr.est,baseb)
  }
  
  
  # IAE 
  Y.p05 = apply(Y.est, 1, quantile, probs=0.05)
  Y.p95 = apply(Y.est, 1, quantile, probs=0.95)
  Y.median = apply(Y.est, 1, median)

  # Transformed versions of these for uniform CI construction
  madn = apply(Y.est, 1, stats::mad)
  helps = apply(Y.est, 1, uniformmax)
  Cvalue = stats::quantile(helps, 0.9)
  
  # Compute Uniform CIs
  Y.u95 = Y.median + Cvalue * madn
  Y.u05 = Y.median - Cvalue * madn
  
  IAE=sum(abs(Y.median-YM)/n)
  
  ###
  # the output
  ###
  
  out=list(
    R.pred=R.pred,
    Beta.pred=Beta.pred,
    Y.p05 = Y.p05,
    Y.p95 = Y.p95,
    Y.u05 = Y.u05,
    Y.u95 = Y.u95,
    Y.median = Y.median,
    IAE=IAE,
    k=KK,
    G=GG[,1:max(KK)],
    In_knots.est=GG_knots[,1:max(KK)],
    AA=AA)
  
  out
  }

gibbs_bernstein_polynomial_est_linkfunction= function(X,Y,
                                                      S=200000,
                                                      L= 20,
                                                      Kmax=500,
                                                      M=1,
                                                      thetak=0.01,
                                                      K=20,
                                                      nu0=0.01,
                                                      s0=0.01,
                                                      A=5,
                                                      d_a=0.5,
                                                      thin=10,
                                                      sigma=0.5,
                                                      alpha1=0.01,
                                                      alpha2=0.01,
                                                      burnin=100000,
                                                      baseb=exp(1),
                                                      logstate=TRUE,
                                                      YM){
  
  n=nrow(X)
  p=ncol(X)
  
  ###
  # The initialization
  ###
  
  # V Z for the weights
  V= beyondWhittle:::vFromP(rep(1 / (L + 1), L))
  Z= seq(from=1 / (2 * K), to = 1 - 1 / (2 * K), length.out = L + 1)
  
  P = bsplinePsd:::pFromV(V)
  G = bsplinePsd:::mixtureWeight(P, Z, K)
  
  
  # B-spline density matrix
  #db.list <- dbspline(w, knots, degree) 
  
  # Initialise r
  r=rbinom(p,1,1/2)
  if(all(r==0)){r[sample(1:p,1)]=1}
  
  #Initialise beta
  mcmc_r_beta_ratio=rexp(1,3)
  mcmc_beta_r_rho=0.2
  beta=rep(0,p)
  beta[which(r==1)]=rnorm(length(which(r==1)),0,1) #unit norm
  beta=beta/c(sqrt(crossprod(beta)));
  if(beta[which(r==1)[1]]<0){beta=-beta} #beta1>0
  
  #Initialise kc
  kc=runif(n,0,3) 
  
  #if log
  if(logstate){
    mu_f=function(x,baseb){
      return(logb(x,baseb))
    } 
  }else{
    mu_f=function(x,baseb){
      return(x)
    }
  }
  
  ###
  # Parameter space
  ###
  KK=Sigma=br=br.star=denbeta=AA=Raa=NULL
  R=Beta=matrix(nr=S,nc=p)
  GG=matrix(nr=S,nc=Kmax)
  Rho=R_beta_ratio=NULL
  
  ptime = proc.time()[1]
  
  ###
  # target parameters updating process
  ###
  for(s in 1:S){
    
    if (s %% 1000 == 0) {
      print(paste("Iteration", s, ",", "Time elapsed", 
                  round(as.numeric(proc.time()[1] - ptime) / 60, 2),
                  "minutes"))
    }
    
    ### step1: K by Metropolis-within-Gibbs
    
    # use random walk for K
    ss = runif(1)
    if (ss < 0.75) {  
      jump = sample(-1:1, 1, prob = rep(1 / 3, 3))
    }else 
    {
      jump = round(rt(1, 1))  #discrete Cauchy
    }
    K.star = K + jump
    while (K.star < (3) || K.star > Kmax) {  # A bit hacky to ensure k doesn't go out of bounds
      if (ss < 0.75) {
        jump = sample(-1:1, 1, prob = rep(1 / 3, 3))  # Ordinary proposal
      }
      else {
        jump = round(rt(1, 1))  # Bold proposal
      }
      K.star = K + jump
    }
    
    if(K.star!=K){# the weights will change
      G.star= bsplinePsd:::mixtureWeight(P, Z, K.star)
      
      ### calculate the log posterior 
      XB=X%*%beta
      s1=s2=0
      w=exp(XB)/(1+exp(XB))
      denbeta = beyondWhittle:::betaBasis_k(w,K,coarsened=0)
      denbeta.star = beyondWhittle:::betaBasis_k(w,K.star,coarsened=0) 
      
      br <- bsplinePsd:::densityMixture(G, denbeta)
      
      br.star <- bsplinePsd:::densityMixture(G.star, denbeta.star)
      
      s1=sum(dnorm(Y,mu_f(A*br.star,baseb),sqrt(8*kc*sigma),log=TRUE))
      s2=sum(dnorm(Y,mu_f(A*br,baseb),sqrt(8*kc*sigma),log=TRUE))
      
      ### the accept or reject for K.star
      lnr=s1-s2-thetak*(K.star^2-K^2)
      if(min(exp(lnr),1)>runif(1)){
        # store the outcome
        K=K.star
        G=G.star
        br=br.star
        denbeta=denbeta.star
      }
    } # K.star!=K
    ### outcome store
    KK[s]=K 
    
    
    ### step2: lambda sample directly
    
    
    ### step3: r by Metropolis-within-Gibbs (q is sampled directly)
    
    # update rj in random order
    for(j in sample(1:p,p)){
      
      # sample q directly
      q=rbeta(1,1+r[j],2-r[j]) 
      rj.star=rbinom(1,1,q) # the proposal
      r.star=NULL
      r.star=r
      r.star[j]=rj.star
      
      ### if Bspline mixture density is NULL
      if(length(denbeta)==0){ 
        XB=X%*%beta
        w=exp(XB)/(1+exp(XB))
        denbeta = beyondWhittle:::betaBasis_k(w,K,coarsened=0)
        br <- bsplinePsd:::densityMixture(G, denbeta)
      }
      
      if(rj.star!=r[j] & sum(r.star)>=1){
        # updating q
        q.star=rbeta(1,1+rj.star,2-rj.star)
        # updating beta accordingly
        beta.star=NULL
        if(sum(r.star)==1){# one dimension
          beta.star=rep(0,p)
          beta.star[which(r.star==1)]=1
        }else
        {
          beta.star=beta
          non_zeron=which(Beta[1:(s-1),j]!=0)
          if(length(non_zeron)<=1){
            beta.star[j]=rj.star*rnorm(1,0,1)
          }else
          {
            beta.star[j]=rj.star*rnorm(1,mean(Beta[non_zeron,j]),sd(Beta[non_zeron,j]))
            #beta.star[j]=rj.star*rnorm(1,mean(Beta[non_zeron,j]),1)
          }
          # the normalisation
          beta.star=beta.star/c(sqrt(crossprod(beta.star)))
          if(beta.star[which(r.star==1)[1]]<0){beta.star=-beta.star}
        }# sum(r.star)>1
        nr=sum(r)
        nr.star=sum(r.star)
        
        ### calculate the log posterior 
        s1=s2=0
        XB.star=X%*%beta.star
        w.star=exp(XB.star)/(1+exp(XB.star))
        denbeta.star = beyondWhittle:::betaBasis_k(w.star,K,coarsened=0) 
        
        br.star <- bsplinePsd:::densityMixture(G, denbeta.star)
        
        s1=sum(dnorm(Y,mu_f(A*br.star,baseb),sqrt(8*kc*sigma),log=TRUE))
        s2=sum(dnorm(Y,mu_f(A*br,baseb),sqrt(8*kc*sigma),log=TRUE))
        
        ### the accept or reject for r.star
        lnr=s1-s2+(r[j])*log(q.star/q)+(1-r[j])*log((1-q.star)/(1-q))+(nr-nr.star)*log(pi)/2+log(gamma(nr.star/2))-log(gamma(nr/2))
        if(min(exp(lnr),1)>runif(1)){
          r[j]=rj.star
          beta=beta.star
          br=br.star
          denbeta=denbeta.star
        }
      }# rj.star!=r[j] & sum(r.star)>=1
    } # j ends
    ### outcome store
    R[s,]=r
    
    
    ### step4: beta by Metropolis-within-Gibbs (the rho is updated; rho is the multivariate normal distribution of non zero beta)
    
    non_zerob=which(r==1)
    if(length(non_zerob)>1){
      
      # the beta updating
      beta.star=rep(0,p)
      # the rho updating
      mcmc_beta_r_rho=max(0,mcmc_beta_r_rho+(0.234-mcmc_r_beta_ratio)/sqrt(s))
      Rho[s]=mcmc_beta_r_rho
      
      # the multivariate normal distribution
      beta.star[non_zerob]=as.vector(mvtnorm:::rmvnorm(1,sqrt(2)*mcmc_beta_r_rho*beta[non_zerob],diag(length(non_zerob))))
      beta.star=beta.star/sqrt(c(crossprod(beta.star[non_zerob])))
      if(beta.star[which(beta.star!=0)[1]]<0){beta.star=-beta.star}
      
      
      ### calculate the log posterior 
      s1=s2=0
      XB.star=X%*%beta.star
      w.star=exp(XB.star)/(1+exp(XB.star))
      denbeta.star = beyondWhittle:::betaBasis_k(w.star,K,coarsened=0) 
      
      br.star <- bsplinePsd:::densityMixture(G, denbeta.star)
      
      s1=sum(dnorm(Y,mu_f(A*br.star,baseb),sqrt(8*kc*sigma),log=TRUE))
      s2=sum(dnorm(Y,mu_f(A*br,baseb),sqrt(8*kc*sigma),log=TRUE))
      
      ### the accept or reject for beta.star
      lnr=s1-s2+(t(beta.star-sqrt(2)*mcmc_beta_r_rho*beta)%*%diag(rep(p))%*%(beta.star-sqrt(2)*mcmc_beta_r_rho*beta)-t(beta-sqrt(2)*mcmc_beta_r_rho*beta.star)%*%diag(rep(p))%*%(beta-sqrt(2)*mcmc_beta_r_rho*beta.star))/2
      # updating ratio immediately
      mcmc_r_beta_ratio=min(exp(lnr),1)
      R_beta_ratio[s]=mcmc_r_beta_ratio
      
      if(mcmc_r_beta_ratio>runif(1)){
        beta=beta.star;
        br=br.star
        denbeta=denbeta.star
      }
    }# length(non_zerob)>1
    Beta[s,]=beta
    
    
    ### the A
    A.s=runif(1,max(0,A-d_a),A+d_a)
    #A.s=rtruncnorm(1,0,Inf,A,2.38/sqrt(p))
    s1=s2=0
    
    s1=sum(dnorm(Y,mu_f(A.s*br,baseb),sqrt(8*kc*sigma),log=TRUE))
    s2=sum(dnorm(Y,mu_f(A*br,baseb),sqrt(8*kc*sigma),log=TRUE))
    
    lnr=s1-s2+alpha2*(A-A.s)+(alpha1-1)*(log(A.s)-log(A))
    raa=min(exp(lnr),1);if(raa>runif(1)){A=A.s}
    if(raa<0.15){d_a=d_a/2}else if(raa>0.65){d_a=d_a*2}
    AA[s]=A;Raa[s]=raa
    
    
    
    
    ### step5: V by Metropolis-within-Gibbs (V[l] l is randomly updated)
    
    for(l in sample(1:L,L)){
      
      # random walk of V[]
      dv=l/(l+2*sqrt(n))
      Vl.star=runif(1,V[l]-dv,V[l]+dv)
      if (Vl.star>1){
        Vl.star=Vl.star-1
      }else if (Vl.star<0){
        Vl.star=Vl.star+1
      }
      V.star=V
      V.star[l]=Vl.star
      
      # G is updated
      P.star = bsplinePsd:::pFromV(V.star)
      G.star = bsplinePsd:::mixtureWeight(P.star, Z, K)
      
      
      br.star <- bsplinePsd:::densityMixture(G.star, denbeta)
      
      
      # the log posterior 
      s1=s2=0
      s1=sum(dnorm(Y,mu_f(A*br.star,baseb),sqrt(8*kc*sigma),log=TRUE))
      s2=sum(dnorm(Y,mu_f(A*br,baseb),sqrt(8*kc*sigma),log=TRUE))
      
      ### the accept or reject for V.star
      lnr=s1-s2+(M-1)*(log(1-Vl.star)-log(1-V[l]))
      if(min(exp(lnr),1)>runif(1)){
        V[l]=Vl.star
        G=G.star
        br=br.star
        P=P.star
      }
    }# l ends
    
    
    ### step6: Z by Metropolis-within-Gibbs (Z[l] l is randomly updated)
    
    for(l in sample(1:(L+1),L+1)){
      
      dz=l/(l+2*sqrt(n))
      Zl.star=runif(1,Z[l]-dz,Z[l]+dz)
      if (Zl.star>1){
        Zl.star=Zl.star-1
      }else if (Zl.star<0){
        Zl.star=Zl.star+1
      }
      Z.star=Z
      Z.star[l]=Zl.star
      
      G.star = bsplinePsd:::mixtureWeight(P, Z.star, K)
      
      br.star <- bsplinePsd:::densityMixture(G.star, denbeta)
      
      # the log posterior 
      s1=s2=0
      s1=sum(dnorm(Y,mu_f(A*br.star,baseb),sqrt(8*kc*sigma),log=TRUE))
      s2=sum(dnorm(Y,mu_f(A*br,baseb),sqrt(8*kc*sigma),log=TRUE))
      
      lnr=s1-s2
      if(min(exp(lnr),1)>runif(1)){
        Z[l]=Zl.star
        G=G.star
        br=br.star
      }
    }# l ends
    # with V,Z ready G.star is stored
    GG[s,1:K]=G
    
    
    
    ### step9: sigma& kc (sampled directly)
    
    mu=t(Y)-mu_f(A*br,baseb)
    sye=mu%*%diag(1/kc)%*%t(mu)
    sigma=rinvgamma(1,nu0+3*n/2,s0+sum(kc)+sye/16)
    Sigma[s]=sigma
    
    
    eta2=2/sigma
    for (i in 1:n){
      eta1=((mu[i])^2)/(8*sigma)
      kc[i]=rgig(1,1/2,eta1,eta2)
    } 
    
    
    
    # END: MCMC loop
  }
  
  ###
  # Which iterations to keep
  ###
  
  keep=seq(burnin+1,S,by=thin)
  KK=KK[keep]
  AA=AA[keep]
  R=R[keep,]
  Beta=Beta[keep,]
  GG=GG[keep,]
  
  
  ###
  # the estimate
  ###
  
  R.pred=Beta.pred=NULL
  # variable selection
  R.pred=apply(R,2,median)
  Beta.pred=apply(Beta,2,median)
  
  # Y.est
  Y.est=matrix(nr=n,nc=length(keep))
  w.pred=exp(X%*%Beta.pred)/(1+exp(X%*%Beta.pred))
  
  for (isample in 1:length(keep)){
    denbeta.est = beyondWhittle:::betaBasis_k(w.pred,KK[isample],coarsened=0)   
    br.est = bsplinePsd:::densityMixture(GG[isample,1:KK[isample]], denbeta.est)
    Y.est[,isample]= mu_f(AA[isample]*br.est,baseb)
  }
  
  
  # IAE 
  Y.p05 = apply(Y.est, 1, quantile, probs=0.05)
  Y.p95 = apply(Y.est, 1, quantile, probs=0.95)
  Y.median = apply(Y.est, 1, median)
  
  # Transformed versions of these for uniform CI construction
  madn = apply(Y.est, 1, stats::mad)
  helps = apply(Y.est, 1, uniformmax)
  Cvalue = stats::quantile(helps, 0.9)
  
  # Compute Uniform CIs
  Y.u95 = Y.median + Cvalue * madn
  Y.u05 = Y.median - Cvalue * madn


  IAE=sum(abs(Y.median-YM)/n)
  
  ###
  # the output
  ###
  
  out=list(
    Beta.pred=Beta.pred,
    Y.p05 = Y.p05,
    Y.p95 = Y.p95,
    Y.u05 = Y.u05,
    Y.u95 = Y.u95,
    Y.median = Y.median,
    IAE=IAE,
    k=KK,
    G=GG[,1:max(KK)],
    AA=AA)
  
  out}
#

library(parallel)
# numCores <- detectCores() - 1
numCores = 50
# f=function(i) {
#   return(gibbs_bspline_est_linkfunction(X1, Y1));
# }

# mclapply(1:10, f, mc.cores = numCores)
cl <- makeCluster(numCores)


clusterExport(cl,'gibbs_bspline_est_linkfunction')
clusterExport(cl,'gibbs_bernstein_polynomial_est_linkfunction')
clusterEvalQ(cl, library(mgcv))
clusterEvalQ(cl, library(mvtnorm))
clusterEvalQ(cl, library(truncnorm))
clusterEvalQ(cl, library(GIGrvg))
clusterEvalQ(cl, library(invgamma))
clusterEvalQ(cl, library(beyondWhittle))
clusterEvalQ(cl, library(bsplinePsd))
clusterEvalQ(cl, library(jsonlite))


w=function(i){
  # library(mgcv);library(mvtnorm);library(truncnorm);library(GIGrvg);library(invgamma);library(bsplinePsd)
  
  X=replicate(10,runif(100,0,1))
  rbeta1=c(1,1,1,0,0,0,0,0,0,0)/sqrt(3)
  ww=X%*%rbeta1
  Y=ww*abs(sin(ww))+0.1*rnorm(100)
  YM=ww*abs(sin(ww))
  outa=gibbs_bspline_est_linkfunction(X,Y,S=200000,burnin=100000,thin=10,YM=YM,logstate=FALSE)
  outb=gibbs_bernstein_polynomial_est_linkfunction(X,Y,S=200000,burnin=100000,thin=10,YM=YM,logstate=FALSE)
  info=list(ww=ww,YM=YM)
  out=list(a=outa,b=outb,info=info)
  
  return(out)
}


results = parLapply(cl, 1:100, w)
write(toJSON(results),'./results/all_results_Ne2s100.json')

w2=function(i){
  # library(mgcv);library(mvtnorm);library(truncnorm);library(GIGrvg);library(invgamma);library(bsplinePsd)
  
  X=replicate(10,runif(200,0,1))
  rbeta1=c(1,1,1,0,0,0,0,0,0,0)/sqrt(3)
  ww=X%*%rbeta1
  Y=ww*abs(sin(ww))+0.1*rnorm(200)
  YM=ww*abs(sin(ww))
  outa=gibbs_bspline_est_linkfunction(X,Y,S=200000,burnin=100000,thin=10,YM=YM,logstate=FALSE)
  outb=gibbs_bernstein_polynomial_est_linkfunction(X,Y,S=200000,burnin=100000,thin=10,YM=YM,logstate=FALSE)
  info=list(ww=ww,YM=YM)
  out=list(a=outa,b=outb,info=info)
  
  return(out)
}


results2 = parLapply(cl, 1:100, w2)
write(toJSON(results2),'./results/all_results_Ne2s200.json')

w3=function(i){
  # library(mgcv);library(mvtnorm);library(truncnorm);library(GIGrvg);library(invgamma);library(bsplinePsd)
  
  X=replicate(10,runif(300,0,1))
  rbeta1=c(1,1,1,0,0,0,0,0,0,0)/sqrt(3)
  ww=X%*%rbeta1
  Y=ww*abs(sin(ww))+0.1*rnorm(300)
  YM=ww*abs(sin(ww))
  outa=gibbs_bspline_est_linkfunction(X,Y,S=200000,burnin=100000,thin=10,YM=YM,logstate=FALSE)
  outb=gibbs_bernstein_polynomial_est_linkfunction(X,Y,S=200000,burnin=100000,thin=10,YM=YM,logstate=FALSE)
  info=list(ww=ww,YM=YM)
  out=list(a=outa,b=outb,info=info)
  
  return(out)
}


results3 = parLapply(cl, 1:100, w3)
write(toJSON(results3),'./results/all_results_Ne2s300.json')
stopCluster(cl)


