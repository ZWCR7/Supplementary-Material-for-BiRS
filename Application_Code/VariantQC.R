library(data.table)
QCUKB = fread('ukb_snp_qc.txt', data.table = F)

for (l in 1:22)
{
  QCChrl = QCUKB[which(QCUKB$chromosome == l), ]
  
  QCChrl = QCChrl[which(QCChrl$array == 2), ]
  
  for (i in 1:95)
  {
    if (i < 10)
    {
      batchqc = paste0('Batch_b00', i, '_qc')
    }
    else
    {
      batchqc = paste0('Batch_b0', i, "_qc")
    }
    
    indxi = which(names(QCChrl) == batchqc)
    
    QCChrl = QCChrl[which(QCChrl[ ,indxi] == 1), ]
  }
  
  for (j in 1:11)
  {
    ukbqc = paste0('UKBiLEVEAX_b', j, "_qc")
    
    indxj = which(names(QCChrl) == ukbqc)
    
    QCChrl = QCChrl[which(QCChrl[ ,indxj] == 1), ]
  }
  
  QCrsidl = QCChrl$rs_id
  
  biml = read_bim(file = paste0('ukb22418/chr', l, '/ukb22418_c', l, '_b0_v2.bim'))
  QCposl = rep(0, length(QCrsidl))
  
  for (k in 1:length(QCrsidl))
  {
    QCposl[k] = biml$pos[which(biml$id == QCrsidl[k])]
  }
  
  save(list = c('QCrsidl', 'QCposl'), file = paste0('QCrsidChr/QCrsidC', l, '.RData'))
}









