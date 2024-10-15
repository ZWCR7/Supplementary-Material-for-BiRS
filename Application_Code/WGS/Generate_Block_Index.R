library(genio)

for (chr in 1:22)
{
  BP = read_bim(paste0('/mnt/project/Bulk/Previous WGS releases/GATK and GraphTyper WGS/GraphTyper population level WGS variants, PLINK format [200k release]/ukb24305_c', chr, '_b0_v1.bim'))$pos
  
  Block.seq = unique(seq(1, length(BP), 5000), length(BP))
  start = Block.seq[-length(Block.seq)]
  end = c(start[-1] - 1, length(BP))
  CHR = rep(chr, length(start))
  
  Block_Index = cbind(CHR, start, end)
  
  save(Block_Index, file = paste0('Blocks-10kb/Block-Chr', chr, '.RData'))
}
