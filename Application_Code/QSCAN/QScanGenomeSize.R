for (chr in 1:22)
{
  load(paste0('SizeControl Sample/Binary_Sample_field40006_disease_C50 Malignant neoplasm of breast_chr', chr, '.RData'))
  sampleid_ns = rownames(SizeGenotype)
  sampleid_ls = strsplit(sampleid_ns, split = "_")
  
  sampleid = rep(0, length(sampleid_ns))  
  for (i in 1:length(sampleid_ns))
  {
    sampleid[i] = sampleid_ls[[i]][1]
  }
  
  rm(SizeGenotype, sampleid_ls, sampleid_ns)
  gc()
  
  load('GLMSampleQC/Binary_Sample_field40006_disease_C50 Malignant neoplasm of breast_covariate_phenotype.RData')
  load('FInput_Female.RData')
  
  chooseind = rep(0, length(sampleid))
  for (i in 1:length(sampleid))
  {
    chooseind[i] = which(FInputid == sampleid[i])
  }
  
  covariateSize = covariate[chooseind, ]
  
  save(covariateSize, file = paste0('SizeControl Sample/Binary_Sample_field40006_disease_C50 Malignant neoplasm of breast_covariate_chr', chr, '.RData'))
  
  rm(list = ls()[-c(which(ls() == "chr"))])
  gc()
}


for (chr in 1:22)
{
  load(paste0('SizeControl Sample/Binary_Sample_field40006_disease_C50 Malignant neoplasm of breast_chr', chr, '.RData'))
  load(paste0('SizeControl Sample/Binary_Sample_field40006_disease_C50 Malignant neoplasm of breast_covariate_chr', chr, '.RData'))
  
  QScanGenSize = function(genotype, phenotype, X, family, Lmax, Lmin, steplength = 1, times = 2000, alpha = 0.05, f = 0)
  {
    maf <- colMeans(genotype)/2
    ## crop the sequence into sub sequence
    folds <- floor(dim(genotype)[2]/4000)
    
    samplesize <- dim(genotype)[1]
    # rank normal transformation
    if(family=="gaussian")
    {
      lmnull <- lm(phenotype~-1+X)
      sigma <- summary(lmnull)$sigma
      
      fam <- 0
      working <- rep(1,samplesize)
      
      mu0 <- lmnull$fitted
    }
    if(family!="gaussian")
    {
      # fit global null model
      glmnull <- glm(phenotype~-1+X,family=family)
      sigma <- sqrt(summary(glmnull)$dispersion)
      
      fam <- 1
      working <- glmnull$weights
      
      mu0 <- glmnull$fitted
    }
    
    
    ## SCANG-O
    L20 <- matrix(0,folds,times)
    res <- c()
    resmost <- c()
    
    weights <- dbeta(maf,1,1)
    subnum <- floor(dim(genotype)[2]/folds)
    
    for(i in 1:folds)
    {
      
      if(i<folds)
      {
        genotypesub <- genotype[,(subnum*(i-1)+1):(i*subnum+Lmax)]
        weightssub <- weights[(subnum*(i-1)+1):(i*subnum+Lmax)]
      }
      if(i==folds)
      {
        genotypesub <- genotype[,(subnum*(i-1)+1):dim(genotype)[2]]
        weightssub <- weights[(subnum*(i-1)+1):dim(genotype)[2]]
      }
      
      
      
      set.seed(19880615)
      threstemp <- Q_SCAN_Thres(genotypesub,X,working,sigma,fam,times,Lmax,Lmin,weightssub)
      
      begid <- subnum*(i-1)+1
      
      ##### SCANG-O
      L20[i,] <- threstemp
      
      emL20 <- apply(L20,2,max)
      th0 <- quantile(emL20,1-alpha)
      
      restemp <- Q_SCAN_Search(genotypesub,X,working,sigma,fam,phenotype,mu0,th0,Lmax,Lmin,begid,f,weightssub)
      res <- rbind(res,restemp$res)
      resmost <- rbind(resmost,restemp$resmost)
      
      
    }
    
    rm(L20)
    gc()
    
    rm(mu0)
    gc()
    
    ## SCANG-O
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
  
  SizeGenotype[which(is.na(SizeGenotype))] = 0
  TotalN = nrow(SizeGenotype) 
  
  SimuSizeScan = function(s)
  {
    set.seed(s)
    
    PMind = sample(1:TotalN, size = TotalN, replace = F)
    phenotype = rep(0, TotalN)
    
    phenotype[PMind[1:ncase]] = 1; phenotype[PMind[-(1:ncase)]] = 0
    
    SizeRES = QScanGenSize(SizeGenotype, phenotype, covariateSize, family = 'binomial', Lmax = 20, Lmin = 4)
    
    return(SizeRES)
  }
  
  nsimu = 1000
  
  cl = makeCluster(20)
  registerDoParallel(cl)
  
  RES = foreach(s = 1:nsimu, .packages = c("QSCAN")) %dopar% SimuSizeScan(s)
  
  stopImplicitCluster()
  stopCluster(cl)
  
  reject = rep(0, 1000)
  for (i in 1:1000)
  {
    if (sum(RES[[i]]$SCAN_res) != 1) reject[i] = 1
  }
  
  save(reject, file = paste0('SizeControl Sample/Results_Size_Scan_disease_C50 Malignant neoplasm of breast_chr', chr, '.RData'))
  
  rm(list = ls()[-c(which(ls() == "chr"))])
  gc()
}

Rej = rep(0, 22)
for (chr in 1:22)
{
  load(paste0('SizeControl Sample/Results_Size_Scan_disease_C50 Malignant neoplasm of breast_chr', chr, '.RData'))
  
  Rej[chr] = sum(reject)
}