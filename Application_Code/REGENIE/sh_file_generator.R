library(genio)

load('FinalVariant.RData')

for (chr in 1:22)
{
  snp = SNP[which(CHR == chr)]

  snpchr = read_bim(paste0('ukb22418/chr', chr, '/ukb22418_c', chr, '_b0_v2.bim'))

  snprm = setdiff(snpchr$id, snp)

  write.table(snprm, file = paste0('ukb22418/chr', chr, '/snplist_rm', chr, '.txt'), row.names = F, col.names = F, quote = F, eol = '\n')

  print(paste0('chr', chr, '--------------------finished'))
}

for (chr in 1:22)
{
  sink(paste0('ukb22418/regenie_fit_chr', chr, '.sh'))
  cat("#!/bin/bash")
  cat("\n")
  cat(paste0('regenie --step 1 --bed chr', chr, '/ukb22418_c', chr, '_b0_v2 --exclude chr', chr, '/snplist_rm', chr, '.txt --covarFile covariate.txt --phenoFile phenotype.txt --bsize 1000 --bt --out chr', chr, '/fit_bin_out'))
  sink()
}


for (chr in 1:22)
{
  sink(paste0('ukb22418/regenie_test_chr', chr, '.sh'))
  cat("#!/bin/bash")
  cat("\n")
  cat(paste0('regenie --step 2 --bed chr', chr, '/ukb22418_c', chr, '_b0_v2 --exclude chr', chr, '/snplist_rm', chr, '.txt --covarFile covariate.txt --phenoFile phenotype.txt --bsize 1000 --bt --firth --approx --pred chr', chr, '/fit_bin_out_pred.list --out chr', chr, '/test_bin_out_firth'))
  sink()
}

library(data.table)

CHR = NULL
SNP = NULL
LogPV = NULL
BP = NULL

for (chr in 1:22)
{
  chrres = read.table(paste0('ukb22418/chr', chr, '/test_bin_out_firth_phenotype.regenie'), header = T)

  CHR = c(CHR, chrres$CHROM); SNP = c(SNP, chrres$ID); LogPV = c(LogPV, chrres$LOG10P); BP = c(BP, chrres$GENPOS)
}

P = 10^(-LogPV)

regenie_res = data.frame(CHR, BP, P, SNP)
write.csv(regenie_res, file = "C50Breast_RegenieGWAS.csv")