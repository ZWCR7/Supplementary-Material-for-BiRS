library(genio)
library(Matrix)
library(seqminer)

load('Binary_Sample_field40006_breast_cancer_covariate_phenotype.RData')
load('KSID.RData')
GenoPrefix = '/mnt/project/Bulk/Previous WGS releases/GATK and GraphTyper WGS/GraphTyper population level WGS variants, PLINK format [200k release]/ukb24305_c19_b0_v1'

wgs_samples = read_fam(paste0(GenoPrefix, '.fam'))
wgs_sample_ids = wgs_samples$id

index1 = match(wgs_sample_ids, KSID)
sample_index = which(!is.na(index1))

sample_id = wgs_sample_ids[sample_index]

index2 = index1[sample_index]
phenoBreast = phenotype[index2]

save(list = c('phenoBreast', 'sample_id', 'sample_index'), file = 'Breast_QC_Sample.RData')

system('dx upload Breast_QC_Sample.RData')