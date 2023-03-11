KSDet = function (result.prelim, genotype, M = 5, region.pos, alpha = 0.05, bigmemory = T, leveraging = T, LD.filter = 0.75) 
{
  chr = 1
  null.model <- F
  result.summary.single <- c()
  variant.info <- c()
  G.now <- c()
  G.before <- c()
  G.after <- c()
  start.index <- 1
  while (start.index < length(region.pos)) {
    start <- region.pos[start.index]
    end <- region.pos[start.index + 1] - 1
     
    start.index <- start.index + 1
    print("extracting genotype data")
    temp.pos <- floor(unique(c(seq(start, end, by = min(25000, ceiling((end - start)/5))), end)))
    G <- genotype[, start:end]
    
    print("processing genotype data")
    if (length(G) == 0) {
      msg <- "Number of variants in the specified range is 0"
      warning(msg, call. = F)
      next
    }
    else {
      if (ncol(G) == 1) {
        msg <- "Number of variants in the specified range is 1"
        warning(msg, call. = F)
        next
      }
    }
    # G <- G[, match(unique(colnames(G)), colnames(G)), drop = F]
  
    MAF <- colMeans(G)/2
    MAC <- colSums(G)
    MAF[MAF > 0.5] <- 1 - MAF[MAF > 0.5]
    MAC[MAF > 0.5] <- nrow(G) * 2 - MAC[MAF > 0.5]
     
    pos <- as.numeric(gsub("^.*\\:", "", colnames(G)))
    p.single <- as.matrix(Get.p(G, result.prelim))
    temp.variant.info <- cbind(colnames(G), p.single, MAF, MAC)
    colnames(temp.variant.info) <- c("pos", "pvalue", "MAF", 
                                     "MAC")
    if (length(LD.filter) != 0) {
      sparse.fit <- sparse.cor(G)
      cor.X <- sparse.fit$cor
      cov.X <- sparse.fit$cov
      Sigma.distance = as.dist(1 - abs(cor.X))
      if (ncol(G) > 1) {
        fit = hclust(Sigma.distance, method = "complete")
        corr_max = LD.filter
        clusters = cutree(fit, h = 1 - corr_max)
      }
      else {
        clusters <- 1
      }
      temp.index <- sample(length(p.single))
      temp.index <- temp.index[match(unique(clusters), 
                                     clusters[temp.index])]
      if (length(temp.index) <= 1) {
        msg <- "Number of variants after LD filtering in the specified range is <=1"
        warning(msg, call. = F)
        next
      }
      G <- G[, temp.index, drop = F]
      p.single <- p.single[temp.index, , drop = F]
      temp.variant.info <- cbind(temp.variant.info, paste0(start, 
                                                           "-", end, "-", clusters))
      colnames(temp.variant.info)[ncol(temp.variant.info)] <- "cluster"
    }
    variant.info <- rbind(variant.info, temp.variant.info)
    pos <- as.numeric(gsub("^.*\\:", "", colnames(G)))
    G <- G[, order(pos), drop = F]
    p.single <- p.single[order(pos), , drop = F]
    MAF <- colMeans(G)/2
    G <- as.matrix(G)
    #G[, MAF > 0.5 & !is.na(MAF)] <- 2 - G[, MAF > 0.5 & !is.na(MAF)]
    MAF <- colMeans(G)/2
    MAC <- colSums(G)
    G <- Matrix(G, sparse = T)
    pos <- as.numeric(gsub("^.*\\:", "", colnames(G)))
    print("generating knockoffs")
    gc()
    if (leveraging == T) {
      G_k <- create.MK(G, pos, M = M, corr_max = 0.75, 
                       maxN.neighbor = Inf, maxBP.neighbor = 0.1 * 10^6, 
                       thres.ultrarare = 0, bigmemory = bigmemory, 
                       R2.thres = 0.75)
    }
    else {
      G_k <- create.MK(G, pos, M = M, corr_max = 0.75, 
                       maxN.neighbor = Inf, maxBP.neighbor = 0.1 * 10^6, 
                       thres.ultrarare = 0, bigmemory = bigmemory, 
                       n.AL = nrow(G), R2.thres = 0.75)
    }
    print("knockoff analysis")
    p.single_k <- c()
    for (k in 1:M) {
      temp.p <- c()
      for (i in 1:ceiling(ncol(G)/1000)) {
        temp.X <- G_k[[k]][, (1 + (i - 1) * 1000):min(ncol(G), 
                                                      i * 1000), drop = F]
        temp.p <- rbind(temp.p, Get.p(temp.X, result.prelim = result.prelim))
      }
      p.single_k <- cbind(p.single_k, temp.p)
      gc()
    }
    
    p.common <- p.single
    p.common_k <- p.single_k
    W <- (-log10(p.common) - apply(-log10(p.common_k), 
                                   1, median)) * (-log10(p.common) >= apply(-log10(p.common_k), 
                                                                            1, max))
    W[is.na(W)] <- 0
    MK.stat <- MK.statistic(-log10(p.common), -log10(p.common_k), 
                            method = "median")
    temp.summary.single <- cbind(chr, pos, pos, pos, pos, MK.stat, W, p.common, p.common_k, MAF)
    colnames(temp.summary.single) <- c("chr", "start", 
                                       "end", "actual_start", "actual_end", "kappa", 
                                       "tau", "W_KS", "P_KS", paste0("P_KS_k", 1:M), 
                                       "MAF")
     
    result.summary.single <- rbind(result.summary.single, temp.summary.single)
  }
  
  ressum = KS.single.summary(result.summary.single, 5)
  ressum = as.data.frame(ressum)
  
  KSDet_res = which(ressum$Qvalue < alpha)
  SPADet_res = which(ressum$P_KS < alpha/ncol(genotype))
  
  return(list(KSDet_res = KSDet_res, SPADet_res = SPADet_res))
}

Get_Liu_PVal.MOD.Lambda<-function(Q.all, lambda, log.p=FALSE){
  param<-Get_Liu_Params_Mod_Lambda(lambda)
  Q.Norm<-(Q.all - param$muQ)/param$sigmaQ
  Q.Norm1<-Q.Norm * param$sigmaX + param$muX
  p.value<- pchisq(Q.Norm1,  df = param$l,ncp=param$d, lower.tail=FALSE, log.p=log.p)
  return(p.value)
}

Get_Liu_Params_Mod_Lambda<-function(lambda){
  ## Helper function for getting the parameters for the null approximation
  
  c1<-rep(0,4)
  for(i in 1:4){
    c1[i]<-sum(lambda^i)
  }
  
  muQ<-c1[1]
  sigmaQ<-sqrt(2 *c1[2])
  s1 = c1[3] / c1[2]^(3/2)
  s2 = c1[4] / c1[2]^2
  
  beta1<-sqrt(8)*s1
  beta2<-12*s2
  type1<-0
  
  #print(c(s1^2,s2))
  if(s1^2 > s2){
    a = 1/(s1 - sqrt(s1^2 - s2))
    d = s1 *a^3 - a^2
    l = a^2 - 2*d
  } else {
    type1<-1
    l = 1/s2
    a = sqrt(l)
    d = 0
  }
  muX <-l+d
  sigmaX<-sqrt(2) *a
  
  re<-list(l=l,d=d,muQ=muQ,muX=muX,sigmaQ=sigmaQ,sigmaX=sigmaX)
  return(re)
}

Get.p.SKAT<-function(score,K,window.matrix,weight,result.prelim){
  
  mu<-result.prelim$nullglm$fitted.values;Y.res<-result.prelim$Y-mu
  X0<-result.prelim$X0;outcome<-result.prelim$out_type
  
  Q<-as.vector(t(score^2)%*%(weight*window.matrix)^2)
  K.temp<-weight*t(weight*K)
  
  p<-rep(NA,length(Q))
  for(i in 1:length(Q)){
    #print(i)
    temp<-K.temp[window.matrix[,i]!=0,window.matrix[,i]!=0]
    if(sum(temp^2)==0){p[i]<-NA;next}
    
    lambda=eigen(temp,symmetric=T,only.values=T)$values
    if(sum(is.na(lambda))!=0){p[i]<-NA;next}
    
    #temp.p<-SKAT_davies(Q[i],lambda,acc=10^(-6))$Qq
    temp.p<-davies(Q[i],lambda,acc=10^(-6))$Qq
    
    if(temp.p > 1 || temp.p <= 0 ){
      temp.p<-Get_Liu_PVal.MOD.Lambda(Q[i],lambda)
    }
    p[i]<-temp.p
  }
  
  return(as.matrix(p))
}



#percentage notation
percent <- function(x, digits = 3, format = "f", ...) {
  paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
}

sparse.cor <- function(x){
  n <- nrow(x)
  cMeans <- colMeans(x)
  covmat <- (as.matrix(crossprod(x)) - n*tcrossprod(cMeans))/(n-1)
  sdvec <- sqrt(diag(covmat))
  cormat <- covmat/tcrossprod(sdvec)
  list(cov=covmat,cor=cormat)
}

sparse.cov.cross <- function(x,y){
  n <- nrow(x)
  cMeans.x <- colMeans(x);cMeans.y <- colMeans(y)
  covmat <- (as.matrix(crossprod(x,y)) - n*tcrossprod(cMeans.x,cMeans.y))/(n-1)
  list(cov=covmat)
}




max.nth<-function(x,n){return(sort(x,partial=length(x)-(n-1))[length(x)-(n-1)])}

Get.p.base<-function(X,result.prelim){
  #X<-Matrix(X)
  mu<-result.prelim$nullglm$fitted.values;Y.res<-result.prelim$Y-mu
  outcome<-result.prelim$out_type
  if(outcome=='D'){v<-mu*(1-mu)}else{v<-rep(as.numeric(var(Y.res)),nrow(X))}
  A<-(t(X)%*%Y.res)^2
  B<-colSums(v*X^2)
  C<-t(X)%*%(v*result.prelim$X0)%*%result.prelim$inv.X0
  D<-t(t(result.prelim$X0)%*%as.matrix(v*X))
  p<-pchisq(as.numeric(A/(B-rowSums(C*D))),df=1,lower.tail=F)
  #p<-pchisq(as.numeric((t(X)%*%Y.res)^2/(apply(X*(v*X),2,sum)-apply(t(X)%*%(v*result.prelim$X0)%*%result.prelim$inv.X0*t(t(result.prelim$X0)%*%as.matrix(v*X)),1,sum))),df=1,lower.tail=F)
  #p[is.na(p)]<-NA
  return(as.matrix(p))
}

Get.p<-function(X,result.prelim){
  #X<-as.matrix(X)
  mu<-result.prelim$nullglm$fitted.values;Y.res<-result.prelim$Y-mu
  outcome<-result.prelim$out_type
  if(outcome=='D'){
    p<-ScoreTest_SPA(t(X),result.prelim$Y,result.prelim$X,method=c("fastSPA"),minmac=-Inf)$p.value
  }else{
    v<-rep(as.numeric(var(Y.res)),nrow(X))
    A<-(t(X)%*%Y.res)^2
    B<-colSums(v*X^2)
    C<-t(X)%*%(v*result.prelim$X0)%*%result.prelim$inv.X0
    D<-t(t(result.prelim$X0)%*%as.matrix(v*X))
    p<-pchisq(as.numeric(A/(B-rowSums(C*D))),df=1,lower.tail=F)
    #p<-pchisq(as.numeric((t(X)%*%Y.res)^2/(apply(X*(v*X),2,sum)-apply(t(X)%*%(v*result.prelim$X0)%*%result.prelim$inv.X0*t(t(result.prelim$X0)%*%as.matrix(v*X)),1,sum))),df=1,lower.tail=F)
  }
  return(as.matrix(p))
}

Get.Z<-function(X,result.prelim){
  #X<-Matrix(X)
  mu<-result.prelim$nullglm$fitted.values;Y.res<-result.prelim$Y-mu
  sd.X<-apply(X,2,sd)
  Z<-t(apply(X,2,scale))%*%as.matrix(scale(Y.res))/sqrt(length(Y.res))
  Z[sd.X==0,]<-0
  return(as.matrix(Z))
}

MK.statistic<-function (T_0,T_k,method='median'){
  T_0<-as.matrix(T_0);T_k<-as.matrix(T_k)
  T.temp<-cbind(T_0,T_k)
  T.temp[is.na(T.temp)]<-0
  
  which.max.alt<-function(x){
    temp.index<-which(x==max(x))
    if(length(temp.index)!=1){return(temp.index[2])}else{return(temp.index[1])}
  }
  kappa<-apply(T.temp,1,which.max.alt)-1
  
  if(method=='max'){tau<-apply(T.temp,1,max)-apply(T.temp,1,max.nth,n=2)}
  if(method=='median'){
    Get.OtherMedian<-function(x){median(x[-which.max(x)])}
    tau<-apply(T.temp,1,max)-apply(T.temp,1,Get.OtherMedian)
  }
  return(cbind(kappa,tau))
}

MK.threshold.byStat<-function (kappa,tau,M,fdr = 0.1,Rej.Bound=10000){
  b<-order(tau,decreasing=T)
  c_0<-kappa[b]==0
  ratio<-c();temp_0<-0
  for(i in 1:length(b)){
    #if(i==1){temp_0=c_0[i]}
    temp_0<-temp_0+c_0[i]
    temp_1<-i-temp_0
    temp_ratio<-(1/M+1/M*temp_1)/max(1,temp_0)
    ratio<-c(ratio,temp_ratio)
    if(i>Rej.Bound){break}
  }
  ok<-which(ratio<=fdr)
  if(length(ok)>0){
    #ok<-ok[which(ok-ok[1]:(ok[1]+length(ok)-1)<=0)]
    return(tau[b][ok[length(ok)]])
  }else{return(Inf)}
}

MK.threshold<-function (T_0,T_k, fdr = 0.1,method='median',Rej.Bound=10000){
  stat<-MK.statistic(T_0,T_k,method=method)
  kappa<-stat[,1];tau<-stat[,2]
  t<-MK.threshold.byStat(kappa,tau,M=ncol(T_k),fdr=fdr,Rej.Bound=Rej.Bound)
  return(t)
}

MK.q.byStat<-function (kappa,tau,M,Rej.Bound=10000){
  b<-order(tau,decreasing=T)
  c_0<-kappa[b]==0
  #calculate ratios for top Rej.Bound tau values
  ratio<-c();temp_0<-0
  for(i in 1:length(b)){
    #if(i==1){temp_0=c_0[i]}
    temp_0<-temp_0+c_0[i]
    temp_1<-i-temp_0
    temp_ratio<-(1/M+1/M*temp_1)/max(1,temp_0)
    ratio<-c(ratio,temp_ratio)
    if(i>Rej.Bound){break}
  }
  #calculate q values for top Rej.Bound values
  q<-rep(1,length(tau));index_bound<-max(which(tau[b]>0))
  for(i in 1:length(b)){
    temp.index<-i:min(length(b),Rej.Bound,index_bound)
    if(length(temp.index)==0){next}
    q[b[i]]<-min(ratio[temp.index])*c_0[i]+1-c_0[i]
    if(i>Rej.Bound){break}
  }
  return(q)
}

Get.cauchy<-function(p){
  p[p>0.99]<-0.99
  is.small<-(p<1e-16) & !is.na(p)
  is.regular<-(p>=1e-16) & !is.na(p)
  temp<-rep(NA,length(p))
  temp[is.small]<-1/p[is.small]/pi
  temp[is.regular]<-as.numeric(tan((0.5-p[is.regular])*pi))
  
  cct.stat<-mean(temp,na.rm=T)
  if(is.na(cct.stat)){return(NA)}
  if(cct.stat>1e+15){return((1/cct.stat)/pi)}else{
    return(1-pcauchy(cct.stat))
  }
}

Get.cauchy.scan<-function(p,window.matrix){
  p[p>0.99]<-0.99
  is.small<-(p<1e-16) & !is.na(p)
  temp<-rep(0,length(p))
  temp[is.small]<-1/p[is.small]/pi
  temp[!is.small]<-as.numeric(tan((0.5-p[!is.small])*pi))
  #window.matrix.MAC10<-(MAC>=10)*window.matrix0
  
  cct.stat<-as.numeric(t(temp)%*%window.matrix/apply(window.matrix,2,sum))
  #cct.stat<-as.numeric(t(temp)%*%window.matrix.MAC10/apply(window.matrix.MAC10,2,sum))
  is.large<-cct.stat>1e+15 & !is.na(cct.stat)
  is.regular<-cct.stat<=1e+15 & !is.na(cct.stat)
  pval<-rep(NA,length(cct.stat))
  pval[is.large]<-(1/cct.stat[is.large])/pi
  pval[is.regular]<-1-pcauchy(cct.stat[is.regular])
  return(pval)
}

Get.p.moment<-function(Q,re.Q){ #Q a A*q matrix of test statistics, re.Q a B*q matrix of resampled test statistics
  re.mean<-apply(re.Q,2,mean)
  re.variance<-apply(re.Q,2,var)
  re.kurtosis<-apply((t(re.Q)-re.mean)^4,1,mean)/re.variance^2-3
  re.df<-(re.kurtosis>0)*12/re.kurtosis+(re.kurtosis<=0)*100000
  re.p<-t(pchisq((t(Q)-re.mean)*sqrt(2*re.df)/sqrt(re.variance)+re.df,re.df,lower.tail=F))
  #re.p[re.p==1]<-0.99
  return(re.p)
}


Impute<-function(Z, impute.method){
  p<-dim(Z)[2]
  if(impute.method =="random"){
    for(i in 1:p){
      IDX<-which(is.na(Z[,i]))
      if(length(IDX) > 0){
        maf1<-mean(Z[-IDX,i])/2
        Z[IDX,i]<-rbinom(length(IDX),2,maf1)
      }
    }
  } else if(impute.method =="fixed"){
    for(i in 1:p){
      IDX<-which(is.na(Z[,i]))
      if(length(IDX) > 0){
        maf1<-mean(Z[-IDX,i])/2
        Z[IDX,i]<-2 * maf1
      }
    }
  } else if(impute.method =="bestguess") {
    for(i in 1:p){
      IDX<-which(is.na(Z[,i]))
      if(length(IDX) > 0){
        maf1<-mean(Z[-IDX,i])/2
        Z[IDX,i]<-round(2 * maf1)
      }
    }
  } else {
    stop("Error: Imputation method shoud be \"fixed\", \"random\" or \"bestguess\" ")
  }
  return(as.matrix(Z))
}


MK.q.byStat<-function (kappa,tau,M,Rej.Bound=10000){
  b<-order(tau,decreasing=T)
  c_0<-kappa[b]==0
  #calculate ratios for top Rej.Bound tau values
  ratio<-c();temp_0<-0
  for(i in 1:length(b)){
    #if(i==1){temp_0=c_0[i]}
    temp_0<-temp_0+c_0[i]
    temp_1<-i-temp_0
    temp_ratio<-(1/M+1/M*temp_1)/max(1,temp_0)
    ratio<-c(ratio,temp_ratio)
    if(i>Rej.Bound){break}
  }
  #calculate q values for top Rej.Bound values
  q<-rep(1,length(tau));index_bound<-max(which(tau[b]>0))
  for(i in 1:length(b)){
    temp.index<-i:min(length(b),Rej.Bound,index_bound)
    if(length(temp.index)==0){next}
    q[b[i]]<-min(ratio[temp.index])*c_0[i]+1-c_0[i]
    if(i>Rej.Bound){break}
  }
  return(q)
}

KS.single.summary<-function(result.single,M){
  
  temp<-result.single[,colnames(result.single)]
  #colnames(temp)<-colnames(result.window)
  
  result<-temp
  result<-result[order(result[,2]),]
  result<-result[order(result[,1]),]
  
  q<-MK.q.byStat(result[,'kappa'],result[,'tau'],M=5)
  result.summary<-cbind(result[,1:5],q,result[,-(1:5)])
  colnames(result.summary)[6]<-'Qvalue'
  
  return(result.summary)
}