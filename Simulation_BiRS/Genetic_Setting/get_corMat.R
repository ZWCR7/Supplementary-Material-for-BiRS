library(BiRS)
library(QSCAN)
library(genio)
library(Matrix)
library(mvnfast)
library(foreach)
library(seqminer)
library(doParallel)
library(KnockoffScreen)

BIMList = read_bim('genotype.bim')

chr = 1
BP = BIMList$pos
RSID = BIMList$id

posmin = min(BP); posmax = max(BP)

ind.start = 1; Block.bed = NULL; Block.ind = NULL; 
pos.end = 0

while(pos.end < posmax)
{
  pos.start = BP[ind.start]
  pos.trunc = pos.start + 50*1000
  
  ind.end = max(which(BP <= pos.trunc))
  pos.end = BP[ind.end]
  
  Block.bed = rbind(Block.bed, c(chr, pos.start, pos.end))
  Block.ind = rbind(Block.ind, c(chr, ind.start, ind.end))
  
  print(ind.end - ind.start)
  ind.start = ind.end + 1
}

# ind.start = 1; Block.bed.small = NULL; Block.ind.small = NULL; 
# pos.end = 0
# 
# while(pos.end < posmax)
# {
#   pos.start = BP[ind.start]
#   pos.trunc = pos.start + 1000
#   
#   ind.end = max(which(BP <= pos.trunc))
#   pos.end = BP[ind.end]
#   
#   Block.bed.small = rbind(Block.bed.small, c(chr, pos.start, pos.end))
#   Block.ind.small = rbind(Block.ind.small, c(chr, ind.start, ind.end))
#   
#   #print(ind.end - ind.start)
#   ind.start = ind.end + 1
# }

FAMList = read_fam('genotype.fam')
n = length(FAMList$fam)

source('Impute_Func.R')
Get_MAF = function(i)
{
  Gi =  -1*(readPlinkToMatrixByIndex(plinkFilePrefix = 'genotype', sampleIndex = 1:n, markerIndex = Block.ind[i, 2]:Block.ind[i, 3]) - 2)
  Gi[Gi < 0 | Gi > 2] = NA
  if (sum(is.na(Gi)) > 0)
  {
    Gi = apply(Gi, 2, Impute_v, impute.method = 'bestguess')
  }
  
  EleVar = apply(Gi, 2, var)
  index_include = which(EleVar != 0)
  
  Gi = Gi[, index_include]
  MAF = colMeans(Gi)/2
  #corg = cor(Gi, y = NULL, use = 'everything', method = 'pearson')
  
  part_num = ceiling(ncol(Gi)/100)
  for (part in 1:part_num)
  {
    st = (part - 1)*100 + 1
    et = min(part*100, ncol(Gi))
    Gip = Gi[, st:et]
    
    corg = cor(Gip, y = NULL)
    
    save(list = c('corg'), file = paste0('corMatrix/block', i, '_part_', part, '_cor.RData'))
  }
  
  save(list = c('MAF'), file = paste0('MAFS/block', i, '_maf.RData'))
  
  return(part_num)
}

cl = makeCluster(40)
registerDoParallel(cl)

aa = Sys.time()
NUMS = foreach(i = 1:40, .packages = c("seqminer", "Matrix", "genio")) %dopar% Get_MAF(i)
bb = Sys.time()

stopImplicitCluster()
stopCluster(cl)

save(list = 'NUMS', file = 'part numbers.RData')
 



# set.seed(20)
# 
# num_region = 5
# #locuind = c(2, 4, 6, 7, 12)
# locuind = sort(sample(1:nrow(Block.bed), size = num_region, replace = F))
# 
# p = length(RSID)
# FAMList = read_fam('Gene_Sample/ukb_imp_chr1_simu.fam')
# n = length(FAMList$fam); FinalID = rep(0, n)
# for (i in 1:n)
# {
#   FinalID[i] = paste0(FAMList$fam[i], '_', FAMList$id[i])
# }
# 
# 
# X1 = rnorm(n, 0, 1); X2 = sample(c(0, 1), size = n, replace = T)
# X = cbind(X1, X2)
# 
# gamma = c(0.5, 0.5)
# 
# ##Continuous Type for Block BiRS
# #cbase = 1e-04
# cv = c(0.2, 0.4, 0.6, 0.8, 1.0)
# for (setc in 1:length(cv))
# {
#   source('BlockBiRS.R')
#   source('BlockQSCAN.R')
#   
#   c = cv[setc]
#   print(paste0('Simulate Signal Level c = ', c, '--------------------Start'))
#   
#   GenoFile = 'Gene_Sample/ukb_imp_chr1_simu.vcf.gz'
#   
#   beta = rep(0, p)
#   set.seed(2)
#   Ycore = rep(0, n)
#   index_causal = c()
#   for (i in 1:num_region)
#   {
#     range = paste0(Block.bed[locuind[i], 1], ":", Block.bed[locuind[i], 2], "-", Block.bed[locuind[i], 3])
#     Gi = Matrix(t(readVCFToMatrixByRange(GenoFile, range, annoType='')[[1]]), sparse = T)
#     Gi[Gi < 0 | Gi > 2] = NA
#     Gi = apply(Gi, 2, Impute_v, impute.method = 'bestguess')
#     
#     MAFi = colMeans(Gi)/2; EleVari = apply(Gi, 2, var)
#     #print(sum(2*MAFi < 0.001))
#     MAFi[which(MAFi == 0)] = 0.0001
#     colBP = as.numeric(gsub("^.*\\:","",colnames(Gi)))
#     
#     betai = rep(0, length(colBP))
#     
#     causal.posstart = sample(colBP, 1); causal.posend = causal.posstart + 10*1000
#     causal.indstart = which(colBP == causal.posstart); causal.indend = max(which(colBP < causal.posend))
#     
#     signalpos = causal.indstart:causal.indend
#     
#     signalpos1 = sample(signalpos, size = ceiling(0.1*length(signalpos)), replace = F)
#     signalpos2 = setdiff(signalpos, signalpos1)
#     signsel = sample(c(-1, 1), size = length(signalpos1), prob = c(0.5, 0.5), replace = T)
#     
#     print(MAFi[signalpos1])
#     
#     # betai[signalpos1] = 1/sqrt(2*MAFi[signalpos1]*(1 - MAFi[signalpos1]))*signsel
#     # betai1 = betai[signalpos1]
#     # 
#     # a = sqrt(c/t(betai1^2 %*% EleVari[signalpos1]))
#     # betai = a[1, 1]*betai
#     
#     betai[signalpos1] = c*abs(log10(MAFi[signalpos1]))
#     
#     #betai[signalpos2] = runif(length(signalpos2), -cbase, cbase)
#     betai[signalpos2] = 0
#     
#     origin.startind = which(BP == min(colBP)); origin.endind = which(BP == max(colBP))
#     index_causal = c(index_causal, origin.startind + signalpos - 1)
#     
#     print(length(signalpos))
#     
#     Ycore = Ycore + Gi %*% betai 
#     beta[origin.startind:origin.endind] = betai
#   }
#   
#   Ycore = X %*% gamma + Ycore
#   nsimu = 100
#   
#   SimuC_BBiRS = function(s)
#   {
#     set.seed(s)
#     inv_logi = function(x) exp(x)/(1 + exp(x))
#     
#     eta = Ycore
#     pi = inv_logi(eta)
#     Y = rbinom(n, size = 1, prob = pi)
#     
#     #Y = Ycore + epsilon
#     phenotype = as.vector(Y)
#     
#     XTrans = scale(X, center = T, scale = F)
#     
#     res_BBiRS = Block_BiRS(phenotype = phenotype, X = XTrans, GenoFile = GenoFile,
#                            Block.bed = Block.bed, method = 'Distribute', IfThres = T, IfScale = F, family = 'binomial', trunc = 6, MB = 1000, alpha = 0.05, bigmemory = T)
#     
#     return(list(res_BBiRS = res_BBiRS))
#   }
#   
#   print('Block BiRS Detection-----------------Start')
#   aa = Sys.time()
#   cl = makeCluster(50)
#   registerDoParallel(cl)
#   
#   RES_BBiRS = foreach(s = 1:nsimu, .packages = c("BiRS", "MASS", "QSCAN", "Matrix", "seqminer", "genio", "KnockoffScreen")) %dopar% SimuC_BBiRS(s)
#   
#   stopImplicitCluster()
#   stopCluster(cl)
#   bb = Sys.time()
#   print(paste0('Block BiRS Detection---------------End, Time = ', bb - aa))
#   #####################################################################################################################################################################
#   
#   SimuC_KS = function(s)
#   {
#     set.seed(s)
#     inv_logi = function(x) exp(x)/(1 + exp(x))
#     
#     eta = Ycore
#     pi = inv_logi(eta)
#     Y = rbinom(n, size = 1, prob = pi)
#     phenotype = as.vector(Y)
#     
#     result.prelim = KS.prelim(phenotype, X = X, id = FinalID, out_type = "D")
#     
#     result.fit = KS.chr(result.prelim = result.prelim, seq.filename = GenoFile, window.bed = window.bed, region.pos = NULL,
#                         tested.pos = BP, excluded.pos = NULL, M = 5, thres.single = 0.01, thres.ultrarare = 1, thres.missing = 0.05,
#                         midout.dir = NULL, temp.dir = NULL, jobtitle = NULL, Gsub.id = NULL, impute.method = "bestguess", bigmemory = T,
#                         leveraging = T, LD.filter = 0.75)
#     
#     result.summary = as.data.frame(KS.summary(result.fit$result.single, result.fit$result.window, 5))
#     indsel = which(result.summary$Qvalue < 0.05)
#     result.significant = result.summary[indsel, ]
#     
#     startind = rep(0, length(indsel)); endind = rep(0, length(indsel))
#     if (nrow(result.significant) != 0)
#     {
#       for (l in 1:nrow(result.significant))
#       {
#         startind[l] = which(BP == result.significant$actual_start[l])
#         endind[l] = which(BP == result.significant$actual_end[l])
#       }
#       
#       res_KS = data.frame(startind, endind)
#     }
#     else
#     {
#       res_KS = 'There is no significant regions.'
#     }
#     
#     return(list(res_KS = res_KS))
#   }
#   
#   print('KnockoffScreen Detection-----------------Start')
#   aa = Sys.time()
#   cl = makeCluster(50)
#   registerDoParallel(cl)
#   
#   RES_KS = foreach(s = 1:nsimu, .packages = c("BiRS", "MASS", "QSCAN", "Matrix", "seqminer", "genio", "KnockoffScreen")) %dopar% SimuC_KS(s)
#   
#   stopImplicitCluster()
#   stopCluster(cl)
#   bb = Sys.time()
#   print(paste0('KnockoffScreen Detection---------------End, Time = ', bb - aa))
#   #################################################################################################################################################################
#   SimuC_QSCAN = function(s)
#   {
#     set.seed(s)
#     inv_logi = function(x) exp(x)/(1 + exp(x))
#     
#     eta = Ycore
#     pi = inv_logi(eta)
#     Y = rbinom(n, size = 1, prob = pi)
#     phenotype = as.vector(Y)
#     
#     res_QSCAN = Block_QSCAN(phenotype = phenotype, X = X, GenoFile = GenoFile, Block.bed = Block.bed, family = 'binomial', Lmax = 80, Lmin = 40)
#     
#     return(list(res_QSCAN = res_QSCAN))
#   }
#   
#   print('QSCAN Detection-----------------Start')
#   aa = Sys.time()
#   cl = makeCluster(50)
#   registerDoParallel(cl)
#   
#   RES_QSCAN = foreach(s = 1:nsimu, .packages = c("BiRS", "MASS", "QSCAN", "Matrix", "seqminer", "genio", "KnockoffScreen")) %dopar% SimuC_QSCAN(s)
#   
#   stopImplicitCluster()
#   stopCluster(cl)
#   bb = Sys.time()
#   print(paste0('QSCAN Detection---------------End, Time = ', bb - aa))
#   
#   
#   save(list = c('RES_BBiRS', 'RES_KS', 'RES_QSCAN', 'beta', 'index_causal'), file = paste0('Gene_RES/Signal_Level_WGS_Lin_Signal_Binarys c = ', c, '.RData'))
#   
#   print(paste0('Simulate Signal Level c = ', c, '--------------------End'))
#   
#   #rm(list = ls()[-c(which(ls() == 'cv'), which(ls() == 'setc'))])
# }
# 
# rm(list = ls())
# 
