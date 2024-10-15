p = 0
for (block in 1:40)
{
  load(paste0('MAFS/block', block, '_maf.RData'))
  pb = length(MAF)
  
  p = p + pb
}

Anno = matrix(rnorm(10*p), 10, p)
rank = colRanks(t(Anno), preserveShape = TRUE)
PHRED = -10*log10(1-rank/dim(rank)[1])

set.seed(2)
dense = 4
locuind = sort(sample(1:40, size = dense, replace = F))

causal_locu = rep(0, 40); causal_locu[locuind] = 1
gamma1 = log(0.01); gamma2 = log(4)

MAFYMat = matrix(0, 7, p)
deltav = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7)

startH = rep(0, 9)

Hnum = 0
for (hh in 1:9)
{
  st.hh = (hh - 1)*4 + 1
  et.hh = hh*4
  
  for (block in st.hh:et.hh)
  {
    load(paste0('MAFS/block', block, '_maf.RData'))
    Hnum = Hnum + length(MAF)
  }
  
  startH[hh] = Hnum
}

startH = c(0, startH)

true_mat = matrix(0, 7, p)
for (j in 1:7)
{
  delta = deltav[j]
  
  true_causal = rep(0, p)
  MAFX = c(); MAFY = c(); startb = 0
  
  set.seed(1)
  for (block in 1:40)
  {
    load(paste0('MAFS/block', block, '_maf.RData'))
    pb = length(MAF)
    
    MAFC = MAF + 0.01
    
    MAFX = c(MAFX, MAFC)
    
    if (causal_locu[block] == 0)
    {
      MAFY = c(MAFY, MAFC)
    }
    else
    {
      length.window = floor(sample(c(1, 1.5, 2), 1)/50*pb)
      
      print(length.window)
      window.st = sample(1:(pb - length.window - 2), 1)
      window.et = window.st + length.window - 1
      
      true_causal[(startb + window.st):(startb + window.et)] = 1
      
      anno_sel = sample(1:10, 5, replace = F)
      anno_block =  Anno[anno_sel, (startb + window.st):(startb + window.et)]
      
      causal_prob = exp(gamma1 + gamma2*apply(anno_block, 2, sum))/(1 + exp(gamma1 + gamma2*apply(anno_block, 2, sum)))
      
      causal_index = rep(0, pb)
      sign = rep(0, pb)
      
      causal_index[window.st:window.et] = rbinom(length.window, 1, causal_prob)
      sign[window.st:window.et] = sample(c(-1, 1), length.window, replace = T)
      
      MAFY = c(MAFY, MAFC + delta*sign*causal_index*MAFC)
    }
    
    startb = startb + pb
  }
  
  MAFYMat[j, ] = MAFY
  true_mat[j, ] = true_causal
}

alpha = 0.05

BiRSDCF_Det = list()
BiRSCL_Det = list()
QSCAN_Det = list()
STAARO_Det = list()
STAARS_Det = list()
STAARB_Det = list()
KnockoffScreen_Det = list()
SSSS_Det = list()

for (seed in 1:1000)
{
  BSDCFcandidate = NULL
  block_index = NULL
  quantileMax = rep(-9999, 1000); quantileMin = rep(-9999, 1000)
  
  startind = c(); endind = c()
  for (hh in 1:10)
  {
    load(paste0('Genetic_Size_Res/BiRS-DCF_delta = ', 0, '_hh = ', hh, '_seed = ', seed, '.RData'))
    
    quantileMax = pmax(quantileMax, reBSD$quantileMax)
    quantileMin = pmax(quantileMin, reBSD$quantileMin)
    
    if (!is.null(reBSD$BSDCF))
    {
      BSDCFb = reBSD$BSDCF
      BSDCFb[, 1] = BSDCFb[, 1] + startH[hh]
      BSDCFb[, 2] = BSDCFb[, 2] + startH[hh]
      
      BSDCFb = cbind(BSDCFb, rep(hh, nrow(BSDCFb)))
      
      BSDCFcandidate = rbind(BSDCFcandidate, BSDCFb)
      block_index = c(block_index, hh) 
    }
  }
  
  thresMax = quantile(quantileMax, 1 - alpha)
  thresMin = quantile(quantileMin, 1 - alpha)
  
  Tn = -10000
  
  if (!is.null(BSDCFcandidate))
  {
    Tn = max(BSDCFcandidate[, 3])
  }
  
  BSDCF_res = NULL
  for (j in 1:length(block_index))
  {
    indexj = which(BSDCFcandidate[, 4] == block_index[j])
    Tnj = max(BSDCFcandidate[indexj, 3])
    
    #BSDCF_res = rbind(BSDCF_res, BSDCFcandidate[indexj, ])
    
    if (Tnj > thresMax)
    {
      BSDCF_res = rbind(BSDCF_res, BSDCFcandidate[indexj, ])
    }
  }
  
  if (!is.null(BSDCF_res))
  {
    index_thres = which(BSDCF_res[, 3] > thresMin)
    BSDCF_res = BSDCF_res[index_thres, ]
    
    startind = BSDCF_res$startind; endind = BSDCF_res$endind
  }
  
  BiRSDCF_Det = c(BiRSDCF_Det, list(data.frame(startind, endind)))
  
  ##############################################################################################################
  startind = c(); endind = c()
  
  for (hh in 1:10)
  {
    load(paste0('Genetic_Size_Res/BiRS-CL_delta = ', 0, '_hh = ', hh, '_seed = ', seed, '.RData'))
    
    if (!is.character(reBSC))
    {
      startind = c(startind, reBSC$BSCL_res$startind + startH[hh])
      endind = c(endind, reBSC$BSCL_res$endind + startH[hh])
    }
  }
  
  BiRSCL_Det = c(BiRSCL_Det, list(data.frame(startind, endind)))
  
  ##############################################################################################################
  startind = c(); endind = c(); Qn = c(); quantileMax = rep(-9999, 1000)
  for (hh in 1:10)
  {
    load(paste0('Genetic_Size_Res/QSCAN_delta = ', 0, '_hh = ', hh, '_seed = ', seed, '.RData'))
    
    quantileMax = pmax(quantileMax, reSCQ$SCAN_thres_boot)
    if (sum(reSCQ$SCAN_res) > 1)
    {
      if (is.null(ncol(reSCQ$SCAN_res)))
      {
        startind = c(startind, reSCQ$SCAN_res[2] + startH[hh])
        endind = c(endind, reSCQ$SCAN_res[3] + startH[hh])
        Qn = c(Qn, reSCQ$SCAN_res[1])
      }
      else
      {
        startind = c(startind, reSCQ$SCAN_res[, 2] + startH[hh])
        endind = c(endind, reSCQ$SCAN_res[, 3] + startH[hh])
        Qn = c(Qn, reSCQ$SCAN_res[, 1])
      }
    }
  }
  
  thresMax = quantile(quantileMax, 1 - alpha)
  if (!is.null(Qn))
  {
    index_sel = which(Qn > thresMax)
    
    if (length(index_sel > 0))
    {
      startind = startind[index_sel]
      endind = endind[index_sel]
    }
    else
    {
      startind = c(); endind = c()
    }
  }
  
  QSCAN_Det = c(QSCAN_Det, list(data.frame(startind, endind)))
  
  ##############################################################################################################
  startind = c(); endind = c(); Qn = c(); quantileMax = rep(-9999, 2000)
  for (hh in 1:10)
  {
    load(paste0('Genetic_Size_Res/SCANG-STAAR_delta = ', 0, '_hh = ', hh, '_seed = ', seed, '.RData'))
    
    quantileMax = pmax(quantileMax, reSCT$SCANG_O_thres_boot)
    if (sum(reSCT$SCANG_O_res) > 1)
    {
      if (is.null(ncol(reSCT$SCANG_O_res)))
      {
        startind = c(startind, reSCT$SCANG_O_res[2] + startH[hh])
        endind = c(endind, reSCT$SCANG_O_res[3] + startH[hh])
        Qn = c(Qn, reSCT$SCANG_O_res[1])
      }
      else
      {
        startind = c(startind, reSCT$SCANG_O_res[, 2] + startH[hh])
        endind = c(endind, reSCT$SCANG_O_res[, 3] + startH[hh])
        Qn = c(Qn, reSCT$SCANG_O_res[, 1])
      }
    }
  }
  
  thresMax = quantile(quantileMax, 1 - alpha)
  if (!is.null(Qn))
  {
    index_sel = which(Qn > thresMax)
    
    if (length(index_sel > 0))
    {
      startind = startind[index_sel]
      endind = endind[index_sel]
    }
    else
    {
      startind = c(); endind = c()
    }
  }
  
  STAARO_Det = c(STAARO_Det, list(data.frame(startind, endind)))
  
  ##############################################################################################################
  startind = c(); endind = c(); Qn = c(); quantileMax = rep(-9999, 2000)
  for (hh in 1:10)
  {
    load(paste0('Genetic_Size_Res/SCANG-STAAR_delta = ', 0, '_hh = ', hh, '_seed = ', seed, '.RData'))
    
    quantileMax = pmax(quantileMax, reSCT$SCANG_S_thres_boot)
    if (sum(reSCT$SCANG_S_res) > 1)
    {
      if (is.null(ncol(reSCT$SCANG_S_res)))
      {
        startind = c(startind, reSCT$SCANG_S_res[2] + startH[hh])
        endind = c(endind, reSCT$SCANG_S_res[3] + startH[hh])
        Qn = c(Qn, reSCT$SCANG_S_res[1])
      }
      else
      {
        startind = c(startind, reSCT$SCANG_S_res[, 2] + startH[hh])
        endind = c(endind, reSCT$SCANG_S_res[, 3] + startH[hh])
        Qn = c(Qn, reSCT$SCANG_S_res[, 1])
      }
    }
  }
  
  thresMax = quantile(quantileMax, 1 - alpha)
  if (!is.null(Qn))
  {
    index_sel = which(Qn > thresMax)
    
    if (length(index_sel > 0))
    {
      startind = startind[index_sel]
      endind = endind[index_sel]
    }
    else
    {
      startind = c(); endind = c()
    }
  }
  
  STAARS_Det = c(STAARS_Det, list(data.frame(startind, endind)))
  
  ##############################################################################################################
  startind = c(); endind = c(); Qn = c(); quantileMax = rep(-9999, 2000)
  for (hh in 1:10)
  {
    load(paste0('Genetic_Size_Res/SCANG-STAAR_delta = ', 0, '_hh = ', hh, '_seed = ', seed, '.RData'))
    
    quantileMax = pmax(quantileMax, reSCT$SCANG_B_thres_boot)
    if (sum(reSCT$SCANG_B_res) > 1)
    {
      if (is.null(ncol(reSCT$SCANG_B_res)))
      {
        startind = c(startind, reSCT$SCANG_B_res[2] + startH[hh])
        endind = c(endind, reSCT$SCANG_B_res[3] + startH[hh])
        Qn = c(Qn, reSCT$SCANG_B_res[1])
      }
      else
      {
        startind = c(startind, reSCT$SCANG_B_res[, 2] + startH[hh])
        endind = c(endind, reSCT$SCANG_B_res[, 3] + startH[hh])
        Qn = c(Qn, reSCT$SCANG_B_res[, 1])
      }
    }
  }
  
  thresMax = quantile(quantileMax, 1 - alpha)
  if (!is.null(Qn))
  {
    index_sel = which(Qn > thresMax)
    
    if (length(index_sel > 0))
    {
      startind = startind[index_sel]
      endind = endind[index_sel]
    }
    else
    {
      startind = c(); endind = c()
    }
  }
  
  STAARB_Det = c(STAARB_Det, list(data.frame(startind, endind)))
  
  ##################################################################################################################
  startind = c(); endind = c()
  
  for (hh in 1:10)
  {
    load(paste0('Genetic_Size_Res/4S_delta = ', 0, '_hh = ', hh, '_seed = ', seed, '.RData'))
    
    if (nrow(reSSS$SSSS_Res) != 0)
    {
      startind = c(startind, reSSS$SSSS_Res$startind + startH[hh])
      endind = c(endind, reSSS$SSSS_Res$endind + startH[hh])
    }
  }
  
  SSSS_Det = c(SSSS_Det, list(data.frame(startind, endind)))
  
}

save(list = c('BiRSDCF_Det'), file = paste0('Genetic_Size_Res/BiRS-DCF_delta = ', 0, '.RData'))
save(list = c('BiRSCL_Det'), file = paste0('Genetic_Size_Res/BiRS-CL_delta = ', 0, '.RData'))
save(list = c('QSCAN_Det'), file = paste0('Genetic_Size_Res/QSCAN_delta = ', 0, '.RData'))
save(list = c('STAARO_Det'), file = paste0('Genetic_Size_Res/STAARO_delta = ', 0, '.RData'))
save(list = c('STAARS_Det'), file = paste0('Genetic_Size_Res/STAARS_delta = ', 0, '.RData'))
save(list = c('STAARB_Det'), file = paste0('Genetic_Size_Res/STAARB_delta = ', 0, '.RData'))
save(list = c('SSSS_Det'), file = paste0('Genetic_Size_Res/4S_delta = ', 0, '.RData'))

sizeBiRSDCF = 0; sizeBiRSCL = 0; sizeQSCAN = 0; sizeSTAARO = 0; sizeSTAARB = 0; sizeSTAARS = 0; sizeSSSS = 0
for (i in 1:1000)
{
  if (nrow(BiRSDCF_Det[[i]]) > 0) sizeBiRSDCF = sizeBiRSDCF + 1
  if (nrow(BiRSCL_Det[[i]]) > 0) sizeBiRSCL = sizeBiRSCL + 1
  if (nrow(QSCAN_Det[[i]]) > 0) sizeQSCAN = sizeQSCAN + 1
  if (nrow(STAARO_Det[[i]]) > 0) sizeSTAARO = sizeSTAARO + 1
  if (nrow(STAARB_Det[[i]]) > 0) sizeSTAARB = sizeSTAARO + 1
  if (nrow(STAARS_Det[[i]]) > 0) sizeSTAARS = sizeSTAARO + 1
  if (nrow(SSSS_Det[[i]]) > 0) sizeSSSS = sizeSSSS + 1
}







