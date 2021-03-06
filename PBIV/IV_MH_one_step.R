##### IV analysis function with bivariate normal random errors
##### with multiple instruments and multiple observed potential confounders 
##### Author: Xuyang Lu <xylu@ucla.edu>
##### Updated on 2013-08-25

# Y = observed time to event or censoring
# d = censoring status, 1=event, 0=censored
# X = covariate of interest
# G = instruments
# U = observed potential confounders
# G is a matrix of n x kG if multiple instruments are used;
# U is a matrix of n x kU if multiple observed confounders are used; the default is NULL
# where n is sample size, kG is number of instruments, kU is number of observed confounders
# A total of (6+kG+2*kU) parameters to estimate: a0, a1 (length=kG), a2 (length=kU), sig1, b0, b1, b2 (length=kU), sig2, and rho
# m = number of iterations
# wid = a vector of the random walk width for (a0,a1,a2,sig1,b0,b1,b2,sig2,rho)
# init = a vector of the initial values for (a0,a1,a2,sig1,b0,b1,b2,sig2,rho), default value is NULL
# prior_1 = a vector of the first parameter of the priors for (a0,a1,a2,sig1,b0,b1,b2,sig2): 
# mean of the normal priors for a0,a1,a2,b0,b1,b2; shape parameter of the inverse-gamma priors for sig1, sig2
# prior_2 = a vector of the second parameter of the priors for (a0,a1,a2,sig1,b0,b1,b2,sig2): 
# variance of the normal priors for a0,a1,a2,b0,b1,b2; scale parameter of the inverse-gamma priors for sig1, sig2

IV_MH = function(Y, d, X, G, U=NULL, m, wid, init=NULL, prior_1, prior_2) {
  # This function generates MCMC samples and acceptance rate for parameters (a0,a1,a2,sig1,b0,b1,b2,sig2,rho);
  
  if (is.null(U)==TRUE) nocov=TRUE  else nocov=FALSE;
  
  if (nocov==FALSE) {
    
    G=as.matrix(G); 
    U=as.matrix(U);
    kG=ncol(G); 
    kU=ncol(U);
    
    if ((length(X)==length(Y) & length(X)==length(d) & length(X)==nrow(G) & length(X)==nrow(U))== FALSE) 
      stop("Variable lengths differ")
    
    # ignoring missing data
    nomis=na.omit(data.frame(Y,d,X,G,U, row.names=NULL))
    Y=as.vector(nomis[,1])
    d=as.vector(nomis[,2])
    X=as.vector(nomis[,3])
    G=as.matrix(nomis[,4:(3+kG)])
    U=as.matrix(nomis[,(4+kG):(3+kG+kU)])
    n=nrow(nomis)
    
    if (length(wid)!= (6+kG+2*kU))
      stop("Length of width vector 'wid' is incorrect")
    
    if (length(prior_1)!= (5+kG+2*kU) | length(prior_2)!= (5+kG+2*kU))
      stop("Length of prior vector is incorrect")
    
    if (!is.null(init) & length(init)!=(6+kG+2*kU))
      stop("Length of initial values 'init' is incorrect")
    
    # assign space to store simulated samples;
    a0r <- sig1r <- b0r <- b1r <- sig2r <- rhor <- rep(0,m);
    a1r <- matrix(0,ncol=kG,nrow=m)
    a2r<- b2r<- matrix(0,ncol=kU,nrow=m)
    #a0r <- sig1r <- b0r <- b1r <- sig2r <- rhor <- rep(0,500);
    #a1r <- matrix(0,ncol=kG,nrow=500)
    #a2r<- b2r<- matrix(0,ncol=kU,nrow=500)
    
    # record the number of acceptances (to estimate acceptance prob.);
    can_a0 <- can_sig1 <- can_b0 <- can_b1 <- can_sig2 <- can_rho <- 0;
    can_a1 <- rep(0,kG)
    can_a2 <- can_b2 <- rep(0, kU)
    
    # Initial values, also as current state;
    if (!is.null(init) & length(init)==(6+kG+2*kU)) {
      a0 <- a0r[1] <- init[1];
      a1 <- a1r[1,] <- init[2:(1+kG)];
      a2 <- a2r[1,] <- init[(2+kG):(1+kG+kU)];
      sig1 <- sig1r[1] <- init[2+kG+kU];
      b0 <- b0r[1] <- init[3+kG+kU];
      b1 <- b1r[1] <- init[4+kG+kU];
      b2 <- b2r[1,] <- init[(5+kG+kU):(4+kG+2*kU)]
      sig2 <- sig2r[1] <- init[5+kG+2*kU];
      rho <- rhor[1] <- init[6+kG+2*kU];
    }
    if (is.null(init)) {
      # frequentist approach with IV to get reasonable initial values;
      fit1=lm(X~G+U, na.action=na.exclude);
      pred=predict(fit1);
      a0 <- a0r[1] <- summary(fit1)$coef[1,1];
      a1 <- a1r[1,] <- summary(fit1)$coef[2:(1+kG),1];
      a2 <- a2r[1,] <- summary(fit1)$coef[(2+kG):(1+kG+kU),1];
      sig1 <- sig1r[1] <- summary(fit1)$sigma^2;
      library(survival)
      fit2=survreg(Surv(Y,d)~pred+U, dist="gaussian")
      b0 <- b0r[1] <- fit2$coef[1];
      b1 <- b1r[1] <- fit2$coef[2];
      b2 <- b2r[1,] <- fit2$coef[3:(2+kU)];
      sig2 <- sig2r[1] <- fit2$scale^2;
      rho <- rhor[1] <- 0;
    }
    
    q=(Y-(b0+b1*X+U%*%b2)-sqrt(sig2/sig1)*rho*(X-a0-G%*%a1-U%*%a2))/sqrt((1-rho^2)*sig2);
    v=(X-a0-G%*%a1-U%*%a2)/sqrt(sig1)
    
    #i = 2
    for (c in 2:m) {
      
      ##### Metropolis?CHastings algorithm for a0 #####
      w=wid[1];
      mu=prior_1[1]; f=prior_2[1];
      a0n <- a0 + runif(1, -w, w); # candidate
      qn=(Y-(b0+b1*X+U%*%b2)-sqrt(sig2/sig1)*rho*(X-a0n-G%*%a1-U%*%a2))/sqrt((1-rho^2)*sig2);
      vn=(X-a0n-G%*%a1-U%*%a2)/sqrt(sig1);
      logf1d=(q^2-qn^2)/2;
      logf2d=(v^2-vn^2)/2;
      Sr=pnorm(qn, lower.tail=FALSE)/exp(log(pnorm(q, lower.tail=FALSE) ));
      logLd=sum(logf1d*d+log(Sr^(1-d))+logf2d);	# difference of log likelihood 
      logpd=((a0-mu)^2-(a0n-mu)^2)/2/f; # difference of log prior
      if(log(runif(1)) < logLd+logpd) { 
        a0 <- a0r[c] <- a0n; 
        can_a0 <- can_a0 + 1; 
        q <- qn;
        v <- vn;
      } else {
        a0r[c] <- a0;
      }
      
      ##### Metropolis?CHastings algorithm for (vector) a1 #####
      for (j in 1:kG) {
        w=wid[1+j];
        mu=prior_1[1+j]; f=prior_2[1+j];
        a1n <- a1; a1n[j] <- a1[j] + runif(1, -w, w); # candidate
        qn=(Y-(b0+b1*X+U%*%b2)-sqrt(sig2/sig1)*rho*(X-a0-G%*%a1n-U%*%a2))/sqrt((1-rho^2)*sig2);
        vn=(X-a0-G%*%a1n-U%*%a2)/sqrt(sig1);
        logf1d=(q^2-qn^2)/2;
        logf2d=(v^2-vn^2)/2;
        Sr=pnorm(qn, lower.tail=FALSE)/exp(log(pnorm(q, lower.tail=FALSE) ));
        logLd=sum(logf1d*d+log(Sr^(1-d))+logf2d);	# difference of log likelihood 
        logpd=((a1[j]-mu)^2-(a1n[j]-mu)^2)/2/f; # difference of log prior
        if(log(runif(1)) < logLd+logpd) {
          a1[j] <- a1r[c,j] <- a1n[j]; 
          can_a1[j] <- can_a1[j] + 1;
          q <- qn;
          v <- vn;
        } else {
          a1r[c,j] <- a1[j];			
        }
      }
      
      ##### Metropolis?CHastings algorithm for (vector) a2 #####
      for (j in 1:kU) {
        w=wid[1+kG+j];
        mu=prior_1[1+kG+j]; f=prior_2[1+kG+j];
        a2n <- a2; a2n[j] <- a2[j] + runif(1, -w, w); # candidate
        qn=(Y-(b0+b1*X+U%*%b2)-sqrt(sig2/sig1)*rho*(X-a0-G%*%a1-U%*%a2n))/sqrt((1-rho^2)*sig2);
        vn=(X-a0-G%*%a1-U%*%a2n)/sqrt(sig1)
        logf1d=(q^2-qn^2)/2;
        logf2d=(v^2-vn^2)/2;
        Sr=pnorm(qn, lower.tail=FALSE)/exp(log(pnorm(q, lower.tail=FALSE) ));
        logLd=sum(logf1d*d+log(Sr^(1-d))+logf2d);	# difference of log likelihood 
        logpd=((a2[j]-mu)^2-(a2n[j]-mu)^2)/2/f; # difference of log prior
        if(log(runif(1)) < logLd+logpd) {
          a2[j] <- a2r[c,j] <- a2n[j]; 
          can_a2[j] <- can_a2[j] + 1;
          q <- qn;
          v <- vn;
        } else {
          a2r[c,j] <- a2[j];
        }
      }
      
      ##### Metropolis?CHastings algorithm for sig1 #####
      w=wid[2+kG+kU];
      g1=prior_1[2+kG+kU]; g2=prior_2[2+kG+kU];
      sig1n <- runif(1, max(sig1-w, 0), sig1+w); # candidate
      qn=(Y-(b0+b1*X+U%*%b2)-sqrt(sig2/sig1n)*rho*(X-a0-G%*%a1-U%*%a2))/sqrt((1-rho^2)*sig2);
      vn=(X-a0-G%*%a1-U%*%a2)/sqrt(sig1n);
      logf1d=(q^2-qn^2)/2;
      logf2d=(v^2-vn^2+log(sig1)-log(sig1n))/2;
      Sr=pnorm(qn, lower.tail=FALSE)/exp(log(pnorm(q, lower.tail=FALSE) ));
      logLd=sum(logf1d*d+log(Sr^(1-d))+logf2d);	# difference of log likelihood 
      logpd=(g1+1)*(log(sig1)-log(sig1n))+g2*(1/sig1-1/sig1n); # difference of log prior
      loghd=log(sig1+w-max(0,sig1-w))-log(sig1n+w-max(0,sig1n-w)); # difference of log proposal
      if(log(runif(1)) < logLd+logpd+loghd) {
        sig1 <- sig1r[c] <- sig1n; 
        can_sig1 <- can_sig1 + 1;
        q <- qn;
        v <- vn;
      } else {
        sig1r[c] <- sig1;
      }
      
      ##### Metropolis?CHastings algorithm for b0 #####
      w=wid[3+kG+kU];
      mu=prior_1[3+kG+kU]; f=prior_2[3+kG+kU];
      b0n <- b0 + runif(1, -w, w); # candidate
      qn=(Y-(b0n+b1*X+U%*%b2)-sqrt(sig2/sig1)*rho*(X-a0-G%*%a1-U%*%a2))/sqrt((1-rho^2)*sig2);
      logfd=(q^2-qn^2)/2;
      Sr=pnorm(qn, lower.tail=FALSE)/exp(log(pnorm(q, lower.tail=FALSE) ));
      logLd=sum(logfd*d+log(Sr^(1-d)));	# difference of log likelihood 
      logpd=((b0-mu)^2-(b0n-mu)^2)/2/f; # difference of log prior
      if(log(runif(1)) < logLd+logpd) {
        b0 <- b0r[c] <- b0n;
        can_b0 <- can_b0 + 1;
        q <- qn;
      } else {
        b0r[c] <- b0;	
      }
      
      ##### Metropolis?CHastings algorithm for b1 #####
      w=wid[4+kG+kU];
      mu=prior_1[4+kG+kU]; f=prior_2[4+kG+kU];
      b1n <- b1 + runif(1, -w, w); # candidate
      qn=(Y-(b0+b1n*X+U%*%b2)-sqrt(sig2/sig1)*rho*(X-a0-G%*%a1-U%*%a2))/sqrt((1-rho^2)*sig2);
      logfd=(q^2-qn^2)/2;
      Sr=pnorm(qn, lower.tail=FALSE)/exp(log(pnorm(q, lower.tail=FALSE) ));
      logLd=sum(logfd*d+log(Sr^(1-d)));	# difference of log likelihood 
      logpd=((b1-mu)^2-(b1n-mu)^2)/2/f; # difference of log prior
      if(log(runif(1)) < logLd+logpd) {
        b1 <- b1r[c] <- b1n;
        can_b1 <- can_b1 + 1;
        q <- qn;
      } else {
        b1r[c] <- b1;
      }
      
      ##### Metropolis?CHastings algorithm for (vector) b2 #####
      for (j in 1:kU) {
        w=wid[4+kG+kU+j];
        mu=prior_1[4+kG+kU+j]; f=prior_2[4+kG+kU+j];
        b2n <- b2; b2n[j] <- b2[j] + runif(1, -w, w); # candidate
        qn=(Y-(b0+b1*X+U%*%b2n)-sqrt(sig2/sig1)*rho*(X-a0-G%*%a1-U%*%a2))/sqrt((1-rho^2)*sig2);
        logfd=(q^2-qn^2)/2;
        Sr=pnorm(qn, lower.tail=FALSE)/exp(log(pnorm(q, lower.tail=FALSE) ));
        logLd=sum(logfd*d+log(Sr^(1-d)));	# difference of log likelihood 
        logpd=((b2[j]-mu)^2-(b2n[j]-mu)^2)/2/f; # difference of log prior
        if(log(runif(1)) < logLd+logpd) {
          b2[j] <- b2r[c,j] <- b2n[j]; 
          can_b2[j] <- can_b2[j] + 1;
          q <- qn;
        } else {
          b2r[c,j] <- b2[j];
        }
      }
      
      ##### Metropolis?CHastings algorithm for sig2 #####
      w=wid[5+kG+2*kU];
      g1=prior_1[5+kG+2*kU]; g2=prior_2[5+kG+2*kU];
      sig2n <- runif(1, max(sig2-w, 0), sig2+w); # candidate
      qn=(Y-(b0+b1*X+U%*%b2)-sqrt(sig2n/sig1)*rho*(X-a0-G%*%a1-U%*%a2))/sqrt((1-rho^2)*sig2n);
      logfd=(q^2-qn^2+log(sig2)-log(sig2n))/2;
      Sr=pnorm(qn, lower.tail=FALSE)/exp(log(pnorm(q, lower.tail=FALSE) ));
      logLd=sum(logfd*d+log(Sr^(1-d)));	# difference of log likelihood  
      logpd=(g1+1)*(log(sig2)-log(sig2n))+g2*(1/sig2-1/sig2n); # difference of log prior
      loghd=log(min(2*w,sig2+w))-log(min(2*w,sig2n+w)); # difference of log proposal
      if(log(runif(1)) < logLd+logpd+loghd) {
        sig2 <- sig2r[c] <- sig2n; 
        can_sig2 <- can_sig2 + 1;
        q <- qn;
      } else {
        sig2r[c] <- sig2;
      }
      
      ##### Metropolis?CHastings algorithm for rho #####
      w=wid[6+kG+2*kU];
      rhon <- runif(1, max(rho-w, -1), min(rho+w,1)); # candidate
      qn=(Y-(b0+b1*X+U%*%b2)-sqrt(sig2/sig1)*rhon*(X-a0-G%*%a1-U%*%a2))/sqrt((1-rhon^2)*sig2);
      logfd=(q^2-qn^2+log(1-rho^2)-log(1-rhon^2))/2;
      Sr=pnorm(qn, lower.tail=FALSE)/exp(log(pnorm(q, lower.tail=FALSE) ));
      logLd=sum(logfd*d+log(Sr^(1-d)));	# difference of log likelihood 
      loghd=log(min(1,rho+w)-max(-1,rho-w))-log(min(rhon+w,1)-max(-1,rhon-w)); # difference of log proposal
      if(log(runif(1)) < logLd+loghd) {
        rho <- rhor[c] <- rhon; 
        can_rho <- can_rho + 1;
        q <- qn;
      } else {
        rhor[c] <- rho;   
      }
      
      #i = i + 1
      if (c %% 100 == 0) {
        print(c)
        print("Acceptance rate: a0 a1 a2 sig1 b0 b1 b2 sig2 rho")
        print(round(c(can_a0, can_a1, can_a2, can_sig1, 
                      can_b0, can_b1, can_b2, can_sig2, can_rho)/(c-1), 4))
        #res <- as.data.frame(list(a0r=a0r,a1r=a1r,a2r=a2r,sig1r=sig1r,b0r=b0r,b1r=b1r,b2r=b2r,sig2r=sig2r,rhor=rhor))
        #write_delim(res, paste("0323temp_", round(c / 500), sep = ""))
        #i = 1
      }
    }
    
    res <- list(a0r=a0r,a1r=a1r,a2r=a2r,sig1r=sig1r,b0r=b0r,b1r=b1r,b2r=b2r,sig2r=sig2r,rhor=rhor,
                ap=c(can_a0, can_a1, can_a2, can_sig1, can_b0, can_b1, can_b2, can_sig2, can_rho)/(m-1))
    return(res)
    
  }
  
  
  if (nocov==TRUE) {
    
    G=as.matrix(G); 
    kG=ncol(G);
    
    if ((length(X)==length(Y) & length(X)==length(d) & length(X)==nrow(G))== FALSE) 
      stop("Variable lengths differ")
    
    # ignoring missing data
    nomis=na.omit(data.frame(Y,d,X,G, row.names=NULL))
    Y=as.vector(nomis[,1])
    d=as.vector(nomis[,2])
    X=as.vector(nomis[,3])
    G=as.matrix(nomis[,4:(3+kG)])
    n=nrow(nomis)
    
    if (length(wid)!= (6+kG))
      stop("Length of width vector 'wid' is incorrect")
    
    if (length(prior_1)!= (5+kG) | length(prior_2)!= (5+kG))
      stop("Length of prior vector is incorrect")
    
    if (!is.null(init) & length(init)!=(6+kG))
      stop("Length of initial values 'init' is incorrect")
    
    # assign space to store simulated samples;
    a0r <- sig1r <- b0r <- b1r <- sig2r <- rhor <- rep(0,m);
    a1r <- matrix(0,ncol=kG,nrow=m)
    
    # record the number of candidates (to estimate acceptance prob);
    can_a0 <- can_sig1 <- can_b0 <- can_b1 <- can_sig2 <- can_rho <- 0;
    can_a1 <- rep(0,kG)
    
    # Initial values, also as current state;
    if (!is.null(init) & length(init)==(6+kG)) {
      a0 <- a0r[1] <- init[1];
      a1 <- a1r[1,] <- init[2:(1+kG)];
      sig1 <- sig1r[1] <- init[2+kG];
      b0 <- b0r[1] <- init[3+kG];
      b1 <- b1r[1] <- init[4+kG];
      sig2 <- sig2r[1]  <- init[5+kG];
      rho <- rhor[1] <- init[6+kG];
    }
    if (is.null(init)) {
      # frequentist approach with IV to get reasonable initial values;
      fit1=lm(X~G, na.action=na.exclude);
      pred=predict(fit1);
      a0 <- a0r[1] <- summary(fit1)$coef[1,1];
      a1 <- a1r[1,] <- summary(fit1)$coef[2:(1+kG),1];
      sig1 <- sig1r[1] <- summary(fit1)$sigma^2;
      library(survival)
      fit2=survreg(Surv(Y,d)~pred, dist="gaussian")
      b0 <- b0r[1] <- fit2$coef[1];
      b1 <- b1r[1] <- fit2$coef[2];
      sig2 <- sig2r[1] <- fit2$scale^2;
      rho <- rhor[1] <- 0;
    }
    
    q=(Y-(b0+b1*X)-sqrt(sig2/sig1)*rho*(X-a0-G%*%a1))/sqrt((1-rho^2)*sig2);
    v=(X-a0-G%*%a1)/sqrt(sig1)
    
    for (i in 2:m) {
      
      ##### Metropolis?CHastings algorithm for a0 #####
      w=wid[1];
      mu=prior_1[1]; f=prior_2[1];
      a0n <- a0 + runif(1, -w, w); # candidate
      qn=(Y-(b0+b1*X)-sqrt(sig2/sig1)*rho*(X-a0n-G%*%a1))/sqrt((1-rho^2)*sig2);
      vn=(X-a0n-G%*%a1)/sqrt(sig1);
      logf1d=(q^2-qn^2)/2;
      logf2d=(v^2-vn^2)/2;
      Sr=pnorm(qn, lower.tail=FALSE)/pnorm(q, lower.tail=FALSE);
      logLd=sum(logf1d*d+log(Sr^(1-d))+logf2d);	# difference of log likelihood 
      logpd=((a0-mu)^2-(a0n-mu)^2)/2/f; # difference of log prior
      if(log(runif(1)) < logLd+logpd) {
        a0 <- a0r[i] <- a0n; 
        can_a0 <- can_a0 + 1; 
        q <- qn;
        v <- vn;
      } else {
        a0r[i] <- a0;
      }
      
      ##### Metropolis?CHastings algorithm for (vector) a1 #####
      for (j in 1:kG) {
        w=wid[1+j];
        mu=prior_1[1+j]; f=prior_2[1+j];
        a1n <- a1; a1n[j] <- a1[j] + runif(1, -w, w); # candidate
        qn=(Y-(b0+b1*X)-sqrt(sig2/sig1)*rho*(X-a0-G%*%a1n))/sqrt((1-rho^2)*sig2);
        vn=(X-a0-G%*%a1n)/sqrt(sig1);
        logf1d=(q^2-qn^2)/2;
        logf2d=(v^2-vn^2)/2;
        Sr=pnorm(qn, lower.tail=FALSE)/pnorm(q, lower.tail=FALSE);
        logLd=sum(logf1d*d+log(Sr^(1-d))+logf2d);	# difference of log likelihood 
        logpd=((a1[j]-mu)^2-(a1n[j]-mu)^2)/2/f; # difference of log prior
        if(log(runif(1)) < logLd+logpd) {
          a1[j] <- a1r[c,j] <- a1n[j]; 
          can_a1[j] <- can_a1[j] + 1;
          q <- qn;
          v <- vn;
        } else {
          a1r[c,j] <- a1[j];			
        }
      }
      
      ##### Metropolis?CHastings algorithm for sig1 #####
      w=wid[2+kG];
      g1=prior_1[2+kG]; g2=prior_2[2+kG];
      sig1n <- runif(1, max(sig1-w, 0), sig1+w); # candidate
      qn=(Y-(b0+b1*X)-sqrt(sig2/sig1n)*rho*(X-a0-G%*%a1))/sqrt((1-rho^2)*sig2);
      vn=(X-a0-G%*%a1)/sqrt(sig1n);
      logf1d=(q^2-qn^2)/2;
      logf2d=(v^2-vn^2+log(sig1)-log(sig1n))/2;
      Sr=pnorm(qn, lower.tail=FALSE)/pnorm(q, lower.tail=FALSE);
      logLd=sum(logf1d*d+log(Sr^(1-d))+logf2d);	# difference of log likelihood 
      logpd=(g1+1)*(log(sig1)-log(sig1n))+g2*(1/sig1-1/sig1n); # difference of log prior
      loghd=log(sig1+w-max(0,sig1-w))-log(sig1n+w-max(0,sig1n-w)); # difference of log proposal
      if(log(runif(1)) < logLd+logpd+loghd) {
        sig1 <- sig1r[i] <- sig1n; 
        can_sig1 <- can_sig1 + 1;
        q <- qn;
        v <- vn;
      } else {
        sig1r[i] <- sig1;
      }
      
      ##### Metropolis?CHastings algorithm for b0 #####
      w=wid[3+kG];
      mu=prior_1[3+kG]; f=prior_2[3+kG];
      b0n <- b0 + runif(1, -w, w); # candidate
      qn=(Y-(b0n+b1*X)-sqrt(sig2/sig1)*rho*(X-a0-G%*%a1))/sqrt((1-rho^2)*sig2);
      logfd=(q^2-qn^2)/2;
      Sr=pnorm(qn, lower.tail=FALSE)/pnorm(q, lower.tail=FALSE);
      logLd=sum(logfd*d+log(Sr^(1-d)));	# difference of log likelihood 
      logpd=((b0-mu)^2-(b0n-mu)^2)/2/f; # difference of log prior
      if(log(runif(1)) < logLd+logpd) {
        b0 <- b0r[i] <- b0n;
        can_b0 <- can_b0 + 1;
        q <- qn;
      } else {
        b0r[i] <- b0;	
      }
      
      ##### Metropolis?CHastings algorithm for b1 #####
      w=wid[4+kG];
      mu=prior_1[4+kG]; f=prior_2[4+kG];
      b1n <- b1 + runif(1, -w, w); # candidate
      qn=(Y-(b0+b1n*X)-sqrt(sig2/sig1)*rho*(X-a0-G%*%a1))/sqrt((1-rho^2)*sig2);
      logfd=(q^2-qn^2)/2;
      Sr=pnorm(qn, lower.tail=FALSE)/pnorm(q, lower.tail=FALSE);
      logLd=sum(logfd*d+log(Sr^(1-d)));	# difference of log likelihood 
      logpd=((b1-mu)^2-(b1n-mu)^2)/2/f; # difference of log prior
      if(log(runif(1)) < logLd+logpd) {
        b1 <- b1r[i] <- b1n;
        can_b1 <- can_b1 + 1;
        q <- qn;
      } else {
        b1r[i] <- b1;
      }
      
      ##### Metropolis?CHastings algorithm for sig2 #####
      w=wid[5+kG];
      g1=prior_1[5+kG]; g2=prior_2[5+kG];
      sig2n <- runif(1, max(sig2-w, 0), sig2+w); # candidate
      qn=(Y-(b0+b1*X)-sqrt(sig2n/sig1)*rho*(X-a0-G%*%a1))/sqrt((1-rho^2)*sig2n);
      logfd=(q^2-qn^2+log(sig2)-log(sig2n))/2;
      Sr=pnorm(qn, lower.tail=FALSE)/pnorm(q, lower.tail=FALSE);
      logLd=sum(logfd*d+log(Sr^(1-d)));	# difference of log likelihood  
      logpd=(g1+1)*(log(sig2)-log(sig2n))+g2*(1/sig2-1/sig2n); # difference of log prior
      loghd=log(min(2*w,sig2+w))-log(min(2*w,sig2n+w)); # difference of log proposal
      if(log(runif(1)) < logLd+logpd+loghd) {
        sig2 <- sig2r[i] <- sig2n; 
        can_sig2 <- can_sig2 + 1;
        q <- qn;
      } else {
        sig2r[i] <- sig2;
      }
      
      ##### Metropolis?CHastings algorithm for rho #####
      w=wid[6+kG];
      rhon <- runif(1, max(rho-w, -1), min(rho+w,1)); # candidate
      qn=(Y-(b0+b1*X)-sqrt(sig2/sig1)*rhon*(X-a0-G%*%a1))/sqrt((1-rhon^2)*sig2);
      logfd=(q^2-qn^2+log(1-rho^2)-log(1-rhon^2))/2;
      Sr=pnorm(qn, lower.tail=FALSE)/pnorm(q, lower.tail=FALSE);
      logLd=sum(logfd*d+log(Sr^(1-d)));	# difference of log likelihood 
      loghd=log(min(1,rho+w)-max(-1,rho-w))-log(min(rhon+w,1)-max(-1,rhon-w)); # difference of log proposal
      if(log(runif(1)) < logLd+loghd) {
        rho <- rhor[i] <- rhon; 
        can_rho <- can_rho + 1;
        q <- qn;
      } else {
        rhor[i] <- rho;   
      }
      
      
      if (i %% 100 == 0) {
        print("Acceptance rate:")
        print(c(can_a0, can_a1, can_sig1, can_b0, can_b1, can_sig2, can_rho)/(m-1))
      }
    }
    
    res <- list(a0r=a0r,a1r=a1r,sig1r=sig1r,b0r=b0r,b1r=b1r,sig2r=sig2r,rhor=rhor,
                ap=c(can_a0, can_a1, can_sig1, can_b0, can_b1, can_sig2, can_rho)/(m-1))
  }
  
  res
}


