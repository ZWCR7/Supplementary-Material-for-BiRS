SNP = names(UnVec)
rm(SneX); rm(UnVec); gc()

DATAPos = NULL
for (chr in 1:22)
{
  load(paste0('QCrsidChr/QCrsidC', chr, '.RData'))
  
  Chrl = rep(chr, length(QCposl))
  
  DATAPos = rbind(DATAPos, data.frame(QCrsidl, QCposl, Chrl))
}

BP = NULL; CHR = NULL
for (i in 1:length(SNP))
{
  indi = which(DATAPos$QCrsidl == SNP[i])
  BP = c(BP, DATAPos$QCposl[indi])
  CHR = c(CHR, DATAPos$Chrl[indi])
}

save(list = c('BP', 'SNP', 'CHR'), file = 'FinalVariant.RData')
#######################################################################################
lenvec = rep(0, 160)

for (i in 1:160)
{
  load(paste0('GLMSampleQC/Binary_Sample_Male_part_', i, '.RData'))
  lenvec[i] = ncol(genotype)
  
  rm(genotype); gc()
}

begid = c(1, cumsum(lenvec)[-160] + 1)

save(begid, file = 'SCANBegid.RData')
