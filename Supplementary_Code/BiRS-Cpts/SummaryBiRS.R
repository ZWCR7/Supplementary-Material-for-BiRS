BiRS_Summary = function(UnVec, SneX, foldlen, trunc = 6, alpha = 0.05, ReMax = 10, COMPU = 'None', Pvm = 'Emp')
{
  p = length(UnVec)
  Signal = rep(0, p)
  split = ceiling(p/foldlen)
  
  S.split = list()
  
  if (split == 1)
  {
    S.split = c(S.split, list(1:p))
  }
  if (split != 1)
  {
    for (i in 1:(split - 1))
    {
      S.split = c(S.split, list(((i-1)*foldlen+1):(i*foldlen)))
    }
    S.split = c(S.split, list(((split-1)*foldlen+1):p))
  }
  
  bmax_prm = quantile(apply(SneX, 1, max), prob = 1 - alpha)
   
  if ((COMPU == 'Point') && (Pvm == 'Asy')) UFC = UnVec
  
  #rm(totalre)
  
  Unprm = rep(0, split)
  for (sp in 1:split)
  {
    Unprm[sp] = max(UnVec[S.split[[sp]]])
  }
  
  reject_ind = rep(0, split)
  Slist = list()
  for (sp in 1:split)
  {
    if (Unprm[sp] > bmax_prm) 
    {
      reject_ind[sp] = 1
      ps = length(S.split[[sp]]); ms = ceiling(ps/2)
      Slist = c(Slist, list(S.split[[sp]][1:ms]))
      Slist = c(Slist, list(S.split[[sp]][(ms + 1):ps]))
    }
  }
  
  if (sum(reject_ind) == 0) return("there is no evidence that there exist signal region")
  
  Signal = BiSearchDCF(UnVec, SneX, Slist, Signal, trunc, MB, alpha)
  
  ind = which(Signal == 1)
  UnVec[ind] = 0
  SneX[ ,ind] = 0
  loop.res = Rejec(UnVec, SneX, alpha)
  loop.ind = loop.res$rej
  loop.bound = loop.res$bd
  
  rm(loop.res)
  
  index = 0
  while ((loop.ind == T) & (index <= ReMax))
  {
    Unprm.loop = rep(0, split)
    for (sp in 1:split)
    {
      Unprm.loop[sp] = max(UnVec[S.split[[sp]]])
    }
    
    Slist.loop = list()
    for (sp in 1:split)
    {
      if (Unprm.loop[sp] > loop.bound) 
      {
        p.loop = length(S.split[[sp]]); m.loop = ceiling(p.loop/2)
        Slist.loop = c(Slist.loop, list(S.split[[sp]][1:m.loop]))
        Slist.loop = c(Slist.loop, list(S.split[[sp]][(m.loop + 1):p.loop]))
      }
    }
    
    Signal = BiSearchDCF(UnVec, SneX, Slist.loop, Signal, trunc, MB, alpha)
    
    ind = which(Signal == 1)
    UnVec[ind] = 0
    SneX[ ,ind] = 0
    loop.res = Rejec(UnVec, SneX, alpha)
    loop.ind = loop.res$rej
    loop.bound = loop.res$bd
    
    rm(loop.res)
    
    index = index + 1
  }
  
  indF = which(Signal == 1)
  nF = length(indF)
  startind = indF[1]; endind = NULL
  for (i in 1:(nF - 1))
  {
    if ((indF[i + 1] - indF[i]) > 1)
    {
      endind = c(endind, indF[i])
      startind = c(startind, indF[i + 1])
    }
  }
  endind = c(endind, indF[nF])
  
  Stat = rep(0, length(startind))
  Pv = rep(1, length(startind))
  
  PvPoint = rep(1, nF)
  StPoint = rep(0, nF)
  
  if (COMPU == 'Region')
  {
    if (is.null(Y))
    {
      for (j in 1:length(startind))
      {
        resj = XFTest(as.matrix(X[, startind[j]:endind[j]]), Y = NULL, MB, alpha)
        Stat[j] = resj$Un
        
        if (Pvm == 'Emp') Pv[j] = resj$Pvalue
      }
    }
    if (!is.null(Y))
    {
      for (j in 1:length(startind))
      {
        resj = XFTest(as.matrix(X[, startind[j]:endind[j]]), as.matrix(Y[, startind[j]:endind[j]]), MB, alpha)
        Stat[j] = resj$Un
        
        if (Pvm == 'Emp') Pv[j] = resj$Pvalue
      }
    }
  }
  
  if (COMPU == 'Point')
  {
    if (Pvm == 'Emp')
    {
      if (is.null(Y))
      {
        for (j in 1:nF)
        {
          resj = XFTest(as.matrix(X[ ,indF[j]]), Y = NULL, MB, alpha)
          StPoint[j] = resj$Un; PvPoint[j] = resj$Pvalue
        }
      }
      if (!is.null(Y))
      {
        for (j in 1:nF)
        {
          resj = XFTest(as.matrix(X[ ,indF[j]]), as.matrix(Y[ ,indF[j]]), MB, alpha)
          StPoint[j] = resj$Un; PvPoint[j] = resj$Pvalue
        }
      }
    }
    
    if (Pvm == 'Asy')
    {
      if (is.null(Y))
      {
        for (j in 1:nF)
        {
          StPoint[j] = UFC[indF[j]]
          PvPoint[j] = pnorm(UFC[indF[j]], 0, sqrt(sd(X[ ,indF[j]])^2), lower.tail = F)
        }
      }
      if (!is.null(Y))
      {
        for (j in 1:nF)
        {
          n = nrow(X); m = nrow(Y)
          StPoint[j] = UFC[indF[j]]
          PvPoint[j] = pnorm(UFC[indF[j]], 0, sqrt(sd(X[ ,indF[j]])^2 + n/m*sd(Y[ ,indF[j]])^2), lower.tail = F)
        }
      }
    }
  }
  
  
  BSDCF_res = data.frame(startind, endind, Stat, Pv)
  return(list(BSDCF_res = BSDCF_res, StPoint = StPoint, PvPoint = PvPoint, COMPU = COMPU, index = index))
}