N=91536
for (seg in c(1:100)){
  print(seg)
  indx.SNP.start=(seg-1)*100+1
  indx.SNP.end=indx.SNP.start+99
  names.ukbdata=paste0('ukb_british_100000sample_10000SNP_maf0.05')
  UKBdata <- names.ukbdata
  geno<- readPlinkToMatrixByIndex(UKBdata, 1:N, indx.SNP.start:indx.SNP.end)
  cor.g=cor(geno,y=NULL,use="everything",method='pearson')
  cor.g[1:5,1:5]
  write.table(cor.g,file=paste0('cor100loci/cor_loci',seg,'.txt'),quote=F,sep='\t',col.names = F,row.names = F)
}
