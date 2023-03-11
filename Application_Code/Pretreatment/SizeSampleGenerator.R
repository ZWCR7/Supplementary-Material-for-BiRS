library(data.table)
library(foreach)
library(Rcpp)
library(BiRS)
library(doParallel)
library(VariantAnnotation)

##################################################################################################################################

Disease = "C50 Malignant neoplasm of breast"
Field = '40006'
Sex = "Female"

lenchr = rep(0, 22)
load('FinalVariant.RData')
for (chr in 1:22)
{
  lenchr[chr] = length(which(CHR == chr))
}

rm(BP, SNP, CHR)

chrratio = sum(lenchr)/lenchr
partnum = c(13, 12, 10, 10, 9, 11, 9, 8, 7, 8, 8, 8, 6, 5, 5, 6, 5, 5, 5, 4, 3, 3)

partend = cumsum(partnum)
partst = partend - partnum + 1

load('GLMSampleQC/Binary_Sample_field40006_disease_C50 Malignant neoplasm of breast_covariate_phenotype.RData')
rm(covariate)

for (chr in 1:22)
{
  CaseSamID = which(phenotype == 1); ControlSamID = which(phenotype == 0)
  ncase = length(CaseSamID); ncontrol = length(ControlSamID)
  
  CaseSel = sample(CaseSamID, size = floor(ncase/chrratio[chr]), replace = F)
  ControlSel = sample(ControlSamID, size = floor(ncontrol/chrratio[chr]), replace = F)
  
  rm(CaseSamID); rm(ControlSamID)
  
  TotalSel = c(CaseSel, ControlSel)
  
  SizeGenotype = NULL
  for (part in partst[chr]:partend[chr])
  {
    name_disease = paste0('GLMSampleQC/Binary_Sample_', Sex, '_part_', part ,'.RData')
    load(file = name_disease)
    
    SizeGenotype = cbind(SizeGenotype, genotype[TotalSel, ])
    rm(genotype); gc()
  }
  
  ncase = floor(ncase/chrratio[chr]); ncontrol = floor(ncontrol/chrratio[chr])
  
  save(list = c('SizeGenotype', 'ncase', 'ncontrol'), file =  paste0('SizeControl Sample/Binary_Sample_field', Field, '_disease_', Disease,
                                                                     '_chr', chr,'.RData'))
  
  print(chr)
}