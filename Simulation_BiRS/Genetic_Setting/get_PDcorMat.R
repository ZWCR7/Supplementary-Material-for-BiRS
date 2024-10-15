load('part numbers.RData')

for (block in 1:40)
{
  #print(paste0('extract block-', block))
  
  numb = NUMS[[block]]
  
  for (part in 1:numb)
  {
    load(paste0('corMatrix/block', block, '_part_', part, '_cor.RData'))
    pp = ncol(corg)
    
    LD = nearPD(corg, corr = T, maxit = 100)$mat
    
    LD = as.matrix(LD)
    colnames(LD) = rownames(LD)
    
    save(list = c('LD'), file = paste0('corMatrixPD/block', block, '_part_', part, '_cor.RData'))
  }
  
}