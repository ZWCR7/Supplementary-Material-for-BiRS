##Significant Chromosome
alpha = 0.05
quantileMax = rep(-9999, 1000); quantileMin = rep(-9999, 1000)
UnVec = rep(0, 22)

for (chr in 1:22)
{
  BSDCFcandidate = NULL; block_index= NULL #BPcandidate = list()
  
  print(paste0('chr', chr, '---start'))
  load(paste0('Blocks-index-50kb/Block-Chr', chr, '.RData'))
  partnum = nrow(Block_Index)
  
  UnChr = rep(0, partnum)
  start_index = 0; #cluster = 1
  for (part in 1:partnum)
  {
    load(paste0('results', chr, '/part', part, '.RData'))
    #load(paste0('mafs', chr, '/part', part, '.RData'))
    
    if (!is.null(RES$BiRSRES))
    {
      quantileMax = pmax(quantileMax, RES$BiRSRES$quantileMax)
      quantileMin = pmax(quantileMin, RES$BiRSRES$quantileMin)
      
      UnChr[part] = max(abs(RES$BiRSRES$UnStat))
    }
  }
  
  # SneX[chr] = quantile(quantileMax, 0.)
  UnVec[chr] = max(UnChr)
  
}

thresMax = quantile(quantileMax, 1 - alpha)
chrsel = which(UnVec > thresMax)
##########################################################################################

##Significant region
alpha = 0.05
for (chr in chrsel)
{
  BSDCFcandidate = NULL; block_index= NULL #BPcandidate = list()
  
  print(paste0('chr', chr, '---start'))
  load(paste0('Blocks-index-50kb/Block-Chr', chr, '.RData'))
  partnum = nrow(Block_Index)
  
  #Spart = matrix(0, 1000, partnum)
  quantileMax = rep(-9999, 1000); quantileMin = rep(-9999, 1000)
  
  start_index = 0; #cluster = 1
  for (part in 1:partnum)
  {
    load(paste0('results', chr, '/part', part, '.RData'))
    #load(paste0('mafs', chr, '/part', part, '.RData'))
    
    if (!is.null(RES$BiRSRES))
    {
      quantileMax = pmax(quantileMax, RES$BiRSRES$quantileMax)
      quantileMin = pmax(quantileMin, RES$BiRSRES$quantileMin)
      
      if (!is.null(RES$BiRSRES$BSDCF))
      {
        BSDCFb = RES$BiRSRES$BSDCF
        
        startbp = NULL; endbp = NULL; bplist = NULL
        for (r in 1:nrow(BSDCFb))
        {
          #bps = RES$bp[BSDCFb[r, 1]:BSDCFb[r, 2]]
          #index_bps = match(bps, bp); mafs = MAFs[index_bps]
          
          #mafs = unname(mafs); clusters = rep(cluster, length(bps))
          
          #bplist = rbind(bplist, data.frame(bps, mafs, clusters))
          
          startbp = c(startbp, as.numeric(strsplit(RES$bp[BSDCFb[r, 1]], ':')[[1]][2]))
          endbp = c(endbp, as.numeric(strsplit(RES$bp[BSDCFb[r, 2]], ':')[[1]][2]))
          
          #cluster = cluster + 1
        }
        
        BSDCFb[, 1] = BSDCFb[, 1] + start_index
        BSDCFb[, 2] = BSDCFb[, 2] + start_index
        
        BSDCFb = cbind(BSDCFb, rep(part, nrow(BSDCFb)))
        BSDCFb = cbind(BSDCFb, startbp, endbp)
        
        #BPcandidate = c(BPcandidate, list(bplist))
        BSDCFcandidate = rbind(BSDCFcandidate, BSDCFb)
        block_index = c(block_index, part)
      }
      
      start_index = start_index + RES$nsnp
    }
  }
  
  # SneX[chr] = quantile(quantileMax, 0.)
  # UnVec[chr] = max(Stat)
  thresMax = quantile(quantileMax, 1 - alpha)
  thresMin = quantile(quantileMin, 1 - alpha)

  Tn = -10000

  if (!is.null(BSDCFcandidate))
  {
    Tn = max(BSDCFcandidate[, 3])
  }
  
  BSDCF_res = NULL; BP_res = NULL
  for (j in 1:length(block_index))
  {
    indexj = which(BSDCFcandidate[, 4] == block_index[j])
    Tnj = max(BSDCFcandidate[indexj, 3])

    if (Tnj > thresMax)
    {
      BSDCF_res = rbind(BSDCF_res, BSDCFcandidate[indexj, ])
      #BP_res = rbind(BP_res, BPcandidate[[j]])
    }
  }

  if (!is.null(BSDCF_res))
  {
    index_thres = which(BSDCF_res[, 3] > thresMin)
    BSDCF_res = BSDCF_res[index_thres, ]
    
    #index_thres1 = unique(BP_res$clusters)[index_thres]
    #index_thres2 = match(BP_res$clusters, index_thres1)
    #index_thres3 = index_thres2[-which(is.na(index_thres2))]
    
    #BP_res = BP_res[-which(is.na(index_thres2)), ]
    
    startind = BSDCF_res$startind; endind = BSDCF_res$endind
  }

  BSDCFR = BSDCF_res[, c(5, 6, 3)]
  write.csv(BSDCFR, file = paste0('chr', chr, '_WGS.csv'))
  #write.csv(BP_res, file = paste0('chr', chr, '_WGS_variant_info.csv'))
  #BiRSDCF_Det = c(BiRSDCF_Det, list(data.frame(startind, endind)))
}

#save(list = c('FrateBSD', 'TrateBSD', 'powerBSD', 'DifTiBSD'), file = 'BSD_Describe.RData')

 