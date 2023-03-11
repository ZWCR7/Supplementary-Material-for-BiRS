########################################################################################################################################
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

# AssesCenter = bd$f.54.0.0
# excluindAC = which((AssesCenter != 11010) & (AssesCenter != 11016) & 
#                      (AssesCenter != 11001) & (AssesCenter != 11017))

exclutotal = unique(c(excluindEth, excluindMR, excluindOut, excluindSex))
bdmid = bd[-exclutotal, ]

#traitssel = read.table('SampleQC.txt')
#eclutraits = NULL

inclusel0 = which((bdmid$f.31.0.0 == 'Female'))
inclusel1 = which((bdmid$f.31.0.0 == 'Female'))

CaseGroup = bdmid[inclusel1, ]
ControlGroup = bdmid[inclusel0, ]

rm(list = ls()[-c(which(ls() == 'CaseGroup'), which(ls() == "ControlGroup"))])
gc()

###################################################################################################################
BinaryQC = function(CaseData, ControlData, field, node, sex)
{
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

  Casefield = CaseData[, colind]
  Controlfield = ControlData[, colind]

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

  CaseInd = apply(Casefield, 1, FindID)
  ControlInd = apply(Controlfield, 1, FindID)

  CaseID = CaseData$f.eid[which(CaseInd == 1)]
  ControlID = ControlData$f.eid[which(ControlInd == 0)]


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

  FCasein = NULL
  FControlin = NULL
  
  for (j in 1:length(CaseID))
  {
    indfindj = which(sampleid == CaseID[j])
    
    if (!isEmpty(indfindj))
    {
      FCasein = c(FCasein, indfindj)
    }
  }
  
  for (k in 1:length(ControlID))
  {
    indfindk = which(sampleid == ControlID[k])
    
    if (!isEmpty(indfindk))
    {
      FControlin = c(FControlin, indfindk)
    }
  }
  
  # for (i in 1:length(sampleid))
  # {
  #   if (sum(CaseID == sampleid[i]) == 1)
  #   {
  #     FCaseid = c(FCaseid, sampleid[i])
  #     FCasein = c(FCasein, i)
  #   }
  # 
  #   if (sum(ControlID == sampleid[i]) == 1)
  #   {
  #     FControlid = c(FControlid, sampleid[i])
  #     FControlin = c(FControlin, i)
  #   }
  # }

  disease = codemap$meaning[which(codemap$node_id == node)]

  save(list = c('FCasein', 'FControlin'), file = paste0('CaseControl SampleQC ', sex, '/Binary_Sample_field', field, '_disease_', disease, '.RData'))
}
##########################################################################################################################

BinaryQC(CaseGroup, ControlGroup, 40006, 500, 'Female')
