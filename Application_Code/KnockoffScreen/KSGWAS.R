load('GLMSampleQC/Binary_Sample_field40006_disease_C50 Malignant neoplasm of breast_covariate_phenotype.RData')
#load('FInput_Female.RData')
load('KSID.RData')

n = length(KSID)
FinalID = rep(0, n)
for (i in 1:n)
{
  FinalID[i] = paste0(KSID[i], '_', KSID[i])
}

#rm(FInputid); rm(FInputin)

load('FinalVariant.RData')

KSDet = function(chr)
{
  seq.filename = paste0('ukb22418/chr', chr, '/chr', chr, '.vcf.gz')
  
  poschr = BP[which(CHR == chr)]
  posmin = min(poschr); posmax = max(poschr)
  
  #poschr = poschr[which(poschr >= posmin & poschr <= posmax)]
  
  window.bed = c()
  num_BP = length(poschr)
  
  indgap = seq(1, num_BP, by = 4)
  num_win = length(indgap)
  
  posbegin = poschr[indgap]
  posend = c(poschr[(indgap - 1)[-1]], posmax)
  window.bed = rbind(window.bed, cbind(chr, posbegin, posend))
  r = nrow(window.bed)
  
  window.bed<-window.bed[order(as.numeric(window.bed[,2])),]
  
  X = covariate[, -1]
  result.prelim = KS.prelim(phenotype, X = X, id = FinalID, out_type = "D")
  
  # midout.dir = paste0('KSResults/chr', chr, '/')
  # temp.dir = paste0('KSResults/chr', chr, '/')
  # jobtitle = paste0('GWASKSChr', chr)
  
  fit = KS.chr(result.prelim = result.prelim, seq.filename = seq.filename, window.bed = window.bed, region.pos = NULL,  
               tested.pos = poschr, excluded.pos = NULL, M = 5, thres.single = 0.01, thres.ultrarare = 25, thres.missing = 0.05, 
               midout.dir = NULL, temp.dir = NULL, jobtitle = NULL, Gsub.id = NULL, impute.method = "fixed", bigmemory = T, 
               leveraging = T, LD.filter = 0.75)
  
  
  save(fit, file = paste0('KSResults/KSfitchr', chr, '.RData'))
  
  return(fit)
}

#res11 = KSDet(11)

cl = makeCluster(22)
registerDoParallel(cl)

aa = Sys.time()
res = foreach(chr = 1:22, .packages = c("VariantAnnotation", "KnockoffScreen")) %dopar% KSDet(chr)
bb = Sys.time()

stopImplicitCluster()
stopCluster(cl)

result.single = NULL; result.window = NULL
for (chr in 1:22)
{
  load(paste0('KSResults/KSfitchr', chr, '.RData'))
  
  result.single = rbind(result.single, fit$result.single)
  result.window = rbind(result.window, fit$result.window)
}

result.summary = as.data.frame(KS.summary(result.single, result.window, 5))

indsel = which(result.summary$Qvalue < 0.05)
result.significant = result.summary[indsel, ]

KnockSGWAS = result.significant[, c(1, 2, 6, 10)]

load('FinalVariant.RData')
snp = NULL
for (chr in 1:22)
{
  poschr = BP[which(CHR == chr)]
  idchr = SNP[which(CHR == chr)]
  
  indchrsel = which(KnockSGWAS$chr == chr)
  poschrsel = KnockSGWAS$start[indchrsel]
  
  indpossel = match(poschrsel, poschr)
  
  snp = c(snp, idchr[indpossel])
}

bp = KnockSGWAS$start
chr = KnockSGWAS$chr
Qv = KnockSGWAS$Qvalue
Pv = KnockSGWAS$P_KS
KSGWASResults = data.frame(chr, snp, bp, Qv, Pv)

write.csv(KSGWASResults, file = 'C50Breast_KSGWAS.csv')

















