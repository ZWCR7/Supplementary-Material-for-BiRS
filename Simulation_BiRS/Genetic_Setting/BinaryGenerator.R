library(data.table)
library(genio)

file.pos = read.table('ukb_simu.pos-1', header = T)
pos = file.pos$CHROM_POS
chr.pos = paste0('1:', file.pos$CHROM_POS)

rm(file.pos); gc()

generate_binary = function(part)
{
  if (part <= 10)
  {
    filename = paste0('ukb_simu_txt/haplotype0', part - 1)
  }
  
  if ((part > 10) & (part <= 90))
  {
    filename = paste0('ukb_simu_txt/haplotype', part - 1)
  }
  
  if (part > 90)
  {
    filename = paste0('ukb_simu_txt/haplotype', 9000 + part - 91)
  }
  
  file.hap = fread(file = filename, header = F)
  gc()
  
  N = length(file.hap$V3)
  haplo_all = c()
  
  for (i in c(1:N))
  {
    haplo_i = strsplit(file.hap$V3[i], " ")[[1]]
    haplo_i = as.numeric(haplo_i)
    haplo_i[which(haplo_i == 2)] = 0
    haplo_all = rbind(haplo_all, haplo_i)
  }
  
  rm(file.hap); gc()
  
  geno.transposed = t(haplo_all[c(TRUE, FALSE), ] + haplo_all[c(FALSE, TRUE), ])
  
  partname = paste0('ukb_simu_binary/genotype', part)
  write_plink(partname, geno.transposed)
  
  ##calculate the genotype frequency
  AF = apply(t(geno.transposed), 2, mean)/2
  rm(geno.transposed); gc()
  
  ##Edit the bim file
  bimfile = fread(file = paste0(partname, '.bim'))
  
  M = dim(bimfile)[1]
  
  bimfile$V2 = chr.pos
  bimfile$V4 = pos
  bimfile$V5[which(AF > 0.5)] = '1'
  bimfile$V6[which(AF > 0.5)] = '2'
  write.table(bimfile, file = paste0(partname, '.bim'), quote = F, sep = '\t', col.names = F, row.names = F)
  
  return(NULL)
}

for (set in 1:20)
{
  aa = Sys.time() 
  set1 = 50*(set - 1) + 1
  set2 = 50*set
  
  cl = makeCluster(50)
  registerDoParallel(cl)
  
  APP = foreach(part = set1:set2, .packages = c("data.table", "genio")) %dopar% generate_binary(part)
  
  stopImplicitCluster()
  stopCluster(cl)
  bb = Sys.time()
  
  print(paste0('set', set, '---time:', bb - aa))
  
}


##Edit fam file
for (part in 1:1000)
{
  partname = paste0('ukb_simu_binary/genotype', part)
  famfile = fread(file = paste0(partname, '.fam'))
  id1 = 50*(part - 1) + 1; id2 = 50*part
  
  famfile$V1 = id1:id2 
  famfile$V2 = id1:id2
  
  write.table(famfile, file = paste0(partname, '.fam'), quote = F, sep = '\t', col.names = F, row.names = F)
}













