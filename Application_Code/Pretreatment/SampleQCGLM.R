bd <- read.table("/home/zhw98/Biobank/phenotype.tab", header=TRUE, sep="\t")
lvl.0009 <- c(0,1)
lbl.0009 <- c("Female","Male")
bd$f.31.0.0 <- ordered(bd$f.31.0.0, levels=lvl.0009, labels=lbl.0009)
lvl.1001 <- c(-3,-1,1,2,3,4,5,6,1001,1002,1003,2001,2002,2003,2004,3001,3002,3003,3004,4001,4002,4003)
lbl.1001 <- c("Prefer not to answer","Do not know","White","Mixed","Asian or Asian British","Black or Black British","Chinese","Other ethnic group","British","Irish","Any other white background","White and Black Caribbean","White and Black African","White and Asian","Any other mixed background","Indian","Pakistani","Bangladeshi","Any other Asian background","Caribbean","African","Any other Black background")
bd$f.21000.0.0 <- ordered(bd$f.21000.0.0, levels=lvl.1001, labels=lbl.1001)
bd$f.21000.1.0 <- ordered(bd$f.21000.1.0, levels=lvl.1001, labels=lbl.1001)
bd$f.21000.2.0 <- ordered(bd$f.21000.2.0, levels=lvl.1001, labels=lbl.1001)
bd$f.22001.0.0 <- ordered(bd$f.22001.0.0, levels=lvl.0009, labels=lbl.0009)
lvl.1002 <- c(1)
lbl.1002 <- c("Caucasian")
bd$f.22006.0.0 <- ordered(bd$f.22006.0.0, levels=lvl.1002, labels=lbl.1002)
lvl.0001 <- c(1)
lbl.0001 <- c("Yes")
bd$f.22027.0.0 <- ordered(bd$f.22027.0.0, levels=lvl.0001, labels=lbl.0001)

load('CTtraits.RData')
load('covariates.RData')

bd = data.frame(bd, bdtraits, covariate)
rm('bdtraits')
rm('covariate')
gc()

Outlier = bd$f.22027.0.0
excluindOut = which(Outlier == "Yes")

ReSex = bd$f.31.0.0; GeSex = bd$f.22001.0.0
Coresindex = (ReSex == GeSex)
excluindSex = which(Coresindex == FALSE)

ReEthMat = matrix(c(bd$f.21000.0.0, bd$f.21000.1.0, bd$f.21000.2.0), 3, 502412, byrow = T)
EthFunc = function(x) 
{
  return(sum((x == "British"), na.rm = T) + sum((x == "Irish"), na.rm = T) + sum((x == "White"), na.rm = T)
         + sum((x == "Any other white background"), na.rm = T))
}
excluindEthRe = which(apply(ReEthMat, 2, EthFunc) == 0)

PC1Eth = bd$f.22009.0.1; PC2Eth = bd$f.22009.0.2
minbound = -0.0144823 - 5*10.5773; maxbound = -0.0144823 + 5*10.5773
excluindEthPC1 = which((PC1Eth < minbound) | (PC1Eth > maxbound)) 
excluindEthPC2 = which((PC2Eth < minbound) | (PC2Eth > maxbound))
excluindEthPC = unique(c(excluindEthPC1, excluindEthPC2))
excluindEth = unique(c(excluindEthRe, excluindEthPC))
rm(excluindEthPC, excluindEthPC1, excluindEthPC2, excluindEthRe)

MissingRate = bd$f.22005.0.0
excluindMR = which(MissingRate > 0.05)

Age = bd$f.21022.0.0
excluindAge = which(Age > 50)

#AssesCenter = bd$f.54.0.0
#excluindAC = which((AssesCenter != 11010) & (AssesCenter != 11016) & 
#(AssesCenter != 11001) & (AssesCenter != 11017))

exclutotal = unique(c(excluindEth, excluindMR, excluindOut, excluindSex))
bdmid = bd[-exclutotal, ]

#traitssel = read.table('SampleQC.txt')
#eclutraits = NULL

inclusel = which((bdmid$f.31.0.0 == 'Male'))
InputData = bdmid[inclusel, ]

# CaseGroup = bdmid[inclusel1, ]
# ControlGroup = bdmid[inclusel0, ]

rm(list = ls()[-c(which(ls() == 'InputData'))])
gc()


field = '40006'
node = 500
Sex = 'Female'
  
if (field == "20002") 
{
  colind = 3:138
  codemap = fread('coding6.tsv', data.table = F)
}

if (field == "41202, 41204") 
{
  colind = 204:470
  codemap = fread('coding19.tsv', data.table = F)
}
if (field == "40006") 
{
  colind = 186:203
  codemap = fread('coding19.tsv', data.table = F)
}

Inputfield = InputData[, colind]

coding1 = codemap$coding[which(codemap$node_id == node)]
coding2 = codemap$coding[which(codemap$parent_id == node)]
coding = c(coding1, coding2)

FindID = function(x)
{
  IN = 0
  
  for (j in 1:length(coding))
  {
    if (sum(x == coding[j], na.rm = T) > 0)
    {
      IN = IN + 1
    }
  }
  
  if (IN > 0) 
  {
    return(1)
  }
  else
  {
    return(0)
  }
  
}

phenotype = apply(Inputfield, 1, FindID)
Sex = as.numeric(InputData$f.31.0.0)
AssesCenter = as.numeric(InputData$f.54.0.0)
Age = as.numeric(InputData$f.21022.0.0)
PCA = cbind(InputData$f.22009.0.1, InputData$f.22009.0.2, InputData$f.22009.0.3, InputData$f.22009.0.4, InputData$f.22009.0.5, InputData$f.22009.0.6, InputData$f.22009.0.7, 
            InputData$f.22009.0.8, InputData$f.22009.0.9, InputData$f.22009.0.10, InputData$f.22009.0.11, InputData$f.22009.0.12, InputData$f.22009.0.13, 
            InputData$f.22009.0.14, InputData$f.22009.0.15, InputData$f.22009.0.16, InputData$f.22009.0.17, InputData$f.22009.0.18, InputData$f.22009.0.19, 
            InputData$f.22009.0.20)

covariate = cbind(Sex, AssesCenter, Age, PCA)
InputID = InputData$f.eid

##########################################################################################################################
library(data.table)

readsample = fread(file = 'ukb22418/chr6/chr6.vcf', skip = 6, nrows = 1, data.table = F, header = F)
readsample = readsample[1, -(1:9)]
samplevec = as.vector(readsample, mode = "character")

samplelist = strsplit(samplevec, "")
sampleid = rep(0, length(samplevec))

for (i in 1:length(samplelist))
{
  charid = as.numeric(samplelist[[i]][1:7])
  sampleid[i] = (1e+06)*charid[1] + (1e+05)*charid[2] + (1e+04)*charid[3] + (1e+03)*charid[4] + (1e+02)*charid[5] + (1e+01)*charid[6] + charid[7]
}
rm(readsample); rm(samplelist); rm(samplevec)

indna = which(is.na(sampleid))
sampleid[indna] = '-1'

FInputin = NULL
FInputid = NULL
for (j in 1:length(InputID))
{
  indfindj = which(sampleid == InputID[j])

  if (!isEmpty(indfindj))
  {
    FInputin = c(FInputin, indfindj)
    FInputid = c(FInputid, sampleid[indfindj])
  }
}

save(list = c("FInputin", "FInputid"), file = 'FInput_Female.RData')

load('FInput_Female.RData')
Fphcoin = rep(0, length(FInputid))
for (j in 1:length(FInputid))
{
  Fphcoin[j] = which(InputID == FInputid[j])
}

save(Fphcoin, file = 'Fphcoin_Female.RData')

phenotype = phenotype[Fphcoin]
covariate = covariate[Fphcoin, ]

disease = codemap$meaning[which(codemap$node_id == node)]

save(list = c('phenotype', 'covariate'), file = paste0('GLMSampleQC/Binary_Sample_field', field, '_disease_', disease,
                                                       '_covariate_phenotype', '.RData'))

##############################################################################################################################
load('FInput_Female.RData')
Sex = 'Female'

vcfconvert = function(vcffile, chr, pos1, pos2, rsid, InputID = NULL)
{
  rngs = GRanges(seqnames = paste0(chr), IRanges(pos1, pos2))
  param = ScanVcfParam(which = rngs)

  vcf = readVcf(file = vcffile, paste0(chr), param)

  chrdata = assay(vcf)

  NAME = rownames(chrdata); indin = NULL
  for (i in 1:length(NAME))
  {
    if (sum(rsid == NAME[i]) > 0) indin = c(indin, i)
  }

  chrdata = chrdata[indin, ]

  if (is.null(nrow(chrdata))) chrdata = as.matrix(t(chrdata))

  trans = function(x)
  {
    if (x == "./.") return(NA)
    if (x == "0/0") return(0)
    if (x == "0/1") return(1)
    if (x == "1/0") return(1)
    if (x == "1/1") return(2)
  }

  if (nrow(chrdata) > 1)
  {
    chrnum = apply(chrdata, c(1, 2), trans)
    rm(chrdata); gc()

    InputSam = t(chrnum[, InputID])
    rm(chrnum);gc()

    return(InputSam)
  }

  if (nrow(chrdata) == 1)
  {
    chrnum = apply(chrdata, c(1, 2), trans)
    rm(chrdata); gc()

    InputSam = chrnum[, InputID]
    rm(chrnum);gc()

    return(InputSam)
  }

  if (nrow(chrdata) == 0)
  {
    return(NULL)
  }
}

partind = 14
for (chr in 2:22)
{
  aa = Sys.time()

  load(paste0('QCrsidChr/QCrsidC', chr, '.RData'))
  L = length(QCposl)
  K = ceiling(L/200)

  Generator = function(i)
  {
    ind1 = 200*(i - 1) + 1
    ind2 = min(200*i, L)

    pos1 = QCposl[ind1]; pos2 = QCposl[ind2]

    vcffile = paste0('ukb22418/chr', chr, '/chr', chr, '.vcf.gz')
    genotype = vcfconvert(vcffile, chr, pos1, pos2, QCrsidl, FInputin)
    gc()

    return(genotype)
  }

  cl = makeCluster(25)
  registerDoParallel(cl)

  aa = Sys.time()
  res = foreach(i = 1:K, .packages = c("data.table", "VariantAnnotation")) %dopar% Generator(i)
  bb = Sys.time()

  stopImplicitCluster()
  stopCluster(cl)


  genotype = NULL; part = 1
  for (l in 1:K)
  {
    if (!is.null(res[[l]]))
    {
      genotype = cbind(genotype, res[[l]])
    }

    if ((ncol(genotype) >= 4000) | (l == K))
    {
      save(list = c('genotype'), file = paste0('GLMSampleQC/Binary_Sample_', Sex, '_part_', partind, '.RData'))
      genotype = NULL
      partind = partind + 1
    }

    print(paste0('part', part, 'finished'))
    gc()

    part = part + 1

  }

  bb = Sys.time()

  print(paste0('chr', chr, 'finished----------time:', bb - aa))
  gc()

}

##########################################################################################################################

