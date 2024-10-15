#' The \code{SigDetDCF}: function for DCF binary search based signal detection control FWER without dimension reduction
#' @param X: sample matrix
#' @param Y: sample matrix for second population, alternatively
#' @param foldlen: The length be splitted, default 4096, if dimension of X or Y less than 4096, suggest the dimension of X
#' @param decay: decay rate for seeded intervals
#' @param trunc: The truncation parameter for the smallest signal region
#' @param MB: multiplier bootstrap size, default 1000
#' @param alpha: the size of test, default 0.05
#' @param ReMax: the maximum number for re-searching
#' @param COMPU: the type for computing p value, c('None', 'Region', 'Point'), default 'None' for not computing; 'Region' for computing
#' region based p value and 'Point' for point based p value.
#' @param Pvm: the type of p value, c('Emp', 'Asy'), default 'Emp' for empirical p value, alternative 'Asy' for asymptotic p value 
#' (only for point based p value)
#' @return \code{BSDCF_res}: the detection results, which include the specific signal region and the corresponding statistics and p-value
#' @return \code{index}: the number of research procedure

source('~/ZW/BiRS-Simulation/SBSDCF.R')

SeedDetDCF = function(X, Y = NULL, foldlen = 4096, decay = 2, trunc = 3, MB = 1000, alpha = 0.05, ReMax = 10, COMPU = 'None', Pvm = 'Emp')
{
  p = ncol(X)
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
  
  if (is.null(Y))
  {
    totalre = XFTest(X, Y = NULL, MB, alpha)
  }
  else
  {
    totalre = XFTest(X, Y, MB, alpha)
  }
  
  UnVec = totalre$UnVec
  SneX = totalre$SneX
  bmax_prm = totalre$bound
  
  if ((COMPU == 'Point') && (Pvm == 'Asy')) UFC = UnVec
  
  rm(totalre)
  
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
      ps = length(S.split[[sp]])
      
      int_length = ps/decay
      
      n_int = ceiling(round(ps/int_length, 14))*2-1		# sometimes very slight numerical inaccuracies
      
      st = S.split[[sp]][1]; et = S.split[[sp]][ps]
      new_candidate = cbind(floor(seq(st, et-int_length, length.out = (n_int))), ceiling(seq(st+int_length-1, et, length.out = (n_int))))
      
      for (j in 1:n_int)
      {
        Slist = c(Slist, list(new_candidate[j, 1]:new_candidate[j, 2]))
      }
    }
  }
  
  if (sum(reject_ind) == 0) return("there is no evidence that there exist signal region")
  
  Signal = SeedSearchDCF(UnVec, SneX, Slist, Signal, decay, trunc)
  
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
        p.loop = length(S.split[[sp]])
        
        int_length.loop = p.loop/decay
        n_int.loop = ceiling(round(p.loop/int_length.loop, 14))*2-1
        
        
        st.loop = S.split[[sp]][1]; et.loop = S.split[[sp]][p.loop]
        new_candidate.loop = cbind(floor(seq(st.loop, et.loop-int_length.loop, length.out = (n_int.loop))), ceiling(seq(st.loop+int_length.loop-1, et.loop, length.out = (n_int.loop))))
        
        for (j in 1:n_int.loop)
        {
          Slist.loop = c(Slist.loop, list(new_candidate.loop[j, 1]:new_candidate.loop[j, 2]))
        }
      }
    }
    
    Signal = SeedSearchDCF(UnVec, SneX, Slist.loop, Signal, decay, trunc)
    
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
  
  
  SBDCF_res = data.frame(startind, endind, Stat, Pv)
  return(list(SBDCF_res = SBDCF_res, StPoint = StPoint, PvPoint = PvPoint, COMPU = COMPU, index = index))
}