load('SCANBegid.RData')

QScanRNT = function(X, phenotype, family)
{
  samplesize = length(phenotype)
  
  if(family == "gaussian")
  {
    lmnull = lm(phenotype ~ -1 + X)
    sigma = summary(lmnull)$sigma
    
    fam = 0
    working = rep(1, samplesize)
    
    mu0 = lmnull$fitted
  }
  
  if(family != "gaussian")
  {
    # fit global null model
    glmnull = glm(phenotype ~ -1 + X, family = family)
    sigma = sqrt(summary(glmnull)$dispersion)
    
    fam = 1
    working = glmnull$weights
    
    mu0 = glmnull$fitted
  }
  
  return(list(working = working, fam = fam, sigma = sigma, mu0 = mu0))
}

QScanGen = function(genotype, phenotype, X, working, sigma, fam, mu0, Lmax, Lmin, begid, steplength = 1, times = 1000, alpha = 0.05, f = 0)
{
  maf = colMeans(genotype)/2
  
  samplesize = dim(genotype)[1]
  
  ## SCANG-O
  res = c()
  resmost = c()
  
  weights = dbeta(maf, 1, 1)
  #subnum <- floor(dim(genotype)[2]/folds)
  
  set.seed(19880615)
  threstemp = Q_SCAN_Thres(genotype, X, working, sigma, fam, times, Lmax, Lmin, weights)
  
  #begid <- subnum*(i-1)+1
  
  ##### SCANG-O
  emL20 = threstemp
  
  th0 = quantile(emL20, 1 - alpha)
  
  restemp = Q_SCAN_Search(genotype, X, working, sigma, fam, phenotype, mu0, th0, Lmax, Lmin, begid, f, weights)
  res = rbind(res, restemp$res)
  resmost = rbind(resmost, restemp$resmost)
  
  rm(mu0)
  gc()
  
  return(list(emL20 = emL20, res = res, resmost = resmost))
}

field = '40006'

load('GLMSampleQC/Binary_Sample_field40006_disease_C50 Malignant neoplasm of breast_covariate_phenotype.RData')

resRNT = QScanRNT(covariate, phenotype, family = 'binomial')

working = resRNT$working; sigma = resRNT$sigma; fam = resRNT$fam; mu0 = resRNT$mu0
rm(resRNT)

Sex = 'Female'
ScanPara = function(part)
{
  load(paste0('GLMSampleQC/Binary_Sample_', Sex, '_part_', part, '.RData'))
  begid = begidv[part]
  
  genotype[which(is.na(genotype))] = 0
  
  RESChrPart = QScanGen(genotype, phenotype, covariate, working, sigma, fam, mu0, Lmax = 20, Lmin = 4, begid)
  
  return(RESChrPart)
}

cl = makeCluster(10)
registerDoParallel(cl)

RES = foreach(par = 1:160, .packages = c("QSCAN")) %dopar% ScanPara(par)

stopImplicitCluster()
stopCluster(cl)

res = c()
resmost = c()
L20 = matrix(0, 1000, 160)

for (par in 1:160)
{
  res = rbind(res, RES[[par]]$res)
  resmost = rbind(resmost, RES[[par]]$resmost)
  L20[, par] = RES[[par]]$emL20
}

emL20 = apply(L20, 1, max)
rm(L20)

Organize = function(emL20, res, resmost, f = 0, alpha = 0.05)
{
  th0 = quantile(emL20, 1 - alpha)
  
  res <- res[res[,1]>th0,]
  if(length(res)==0)
  {
    res <- c(0,0,0,1)
  }
  if(length(res)>4)
  {
    res <- regionfilter(res,f)
    if(length(res)==4)
    {
      res[4] <- mean(emL20>res[1])
    }else
    {
      res[,4] <- apply(res,1, function(z) mean(emL20>z[1]))
    }
  }
  
  mostnum <- which.max(resmost[,1])
  resmost <- resmost[mostnum,]
  resmost[4] <- mean(emL20>resmost[1])
  
  Lst <- list(SCAN_res=res,SCAN_top1=resmost,SCAN_thres=th0,SCAN_thres_boot=emL20)
  
  return(Lst)
}

ScanRES = Organize(emL20, res, resmost)
ScanInd = ScanRES$SCAN_res

num = 0
for (i in 1:nrow(ScanInd))
{
  num = num + ScanInd[i, 3] - ScanInd[i, 2] + 1
}

save(ScanRES, file = 'Scan_Results_disease_C50 Malignant neoplasm of breast_Sex_Female.RData')
###################################################################################################
rm(list = ls())

load('Scan_Results_disease_C50 Malignant neoplasm of breast_Sex_Female.RData')

load('FinalVariant.RData')

rsidsel = NULL; possel = NULL; chrsel = NULL; cluster = NULL
for (i in 1:nrow(ScanRES$SCAN_res))
{
  ind1 = ScanRES$SCAN_res[i, 2]; ind2 = ScanRES$SCAN_res[i, 3]
  len = ind2 - ind1 + 1
  cluster = c(cluster, rep(i, len))
  rsidsel = c(rsidsel, SNP[ind1:ind2])
  possel = c(possel, BP[ind1:ind2])
  chrsel = c(chrsel, CHR[ind1:ind2])
}

write.csv(data.frame(rsidsel, possel, chrsel, cluster), file = 'C50Breast_Scan.csv')
