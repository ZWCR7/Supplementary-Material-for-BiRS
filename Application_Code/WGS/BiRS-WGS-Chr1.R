CHR = 1

library(BiRS)
library(genio)
library(Matrix)
library(seqminer)
library(foreach)
library(doParallel)

source('BlockBiRS_WGS.R')
load('Breast_QC_Sample.RData')
load(paste0('Block-Chr', CHR, '.RData'))

partnum = nrow(Block_Index)

GenoPrefix = paste0('/mnt/project/Bulk/Previous WGS releases/GATK and GraphTyper WGS/GraphTyper population level WGS variants, PLINK format [200k release]/ukb24305_c', CHR, '_b0_v1')

Parallel_BBiRS = function(part)
{
  RES = Block_BiRS(phenotype = phenoBreast, GenoPrefix = GenoPrefix, Block.ind = Block_Index[part, ], subsample = sample_index, trunc = 4, alpha = 0.05)
  save(list = c('RES'), file = paste0('results', CHR, '/part', part, '.RData'))
}

TT = ceiling(partnum/50)
for (tt in 1:TT)
{
  aa = Sys.time()
  start_tt = 50*(tt - 1) + 1
  end_tt = min(50*tt, partnum)

  cl = makeCluster(50)
  registerDoParallel(cl)

  RES = foreach(part = start_tt:end_tt, .packages = c('genio', 'BiRS', 'seqminer', 'Matrix')) %dopar% Parallel_BBiRS(part)

  stopImplicitCluster()
  stopCluster(cl)
  bb = Sys.time()

  print(paste0('Big part', tt, 'time = ', bb - aa))
}

system(paste0('tar -czvf results', CHR, '.tar.gz results', CHR))
system(paste0('dx upload results', CHR, '.tar.gz'))
