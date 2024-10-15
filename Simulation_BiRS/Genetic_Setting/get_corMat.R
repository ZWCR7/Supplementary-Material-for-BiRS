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
