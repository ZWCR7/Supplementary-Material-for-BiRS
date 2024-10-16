BiRSDCF_T = matrix(0, 10, 25)
BiRSCL_T = matrix(0, 10, 25)
QSCAN_T = matrix(0, 10, 25)
STAAR_T = matrix(0, 10, 25)
KnockoffScreen_T = matrix(0, 10, 25)
SSSS_T = matrix(0, 10, 25)
SCANDCF_T = matrix(0, 10, 25)
DCF_T = matrix(0, 10, 25)

delta = 0.7

for (seed in 1:25)
{
  for (hh in 1:10)
  {
    print(paste0('seed = ', seed, '_hh = ', hh))
    load(paste0('Genetic_RES/BiRS-DCF_delta = ', delta, '_hh = ', hh, '_seed = ', seed, '.RData'))
    #load(paste0('Genetic_RES/BiRS-CL_delta = ', delta, '_hh = ', hh, '_seed = ', seed, '.RData'))
    load(paste0('Genetic_RES/QSCAN_delta = ', delta, '_hh = ', hh, '_seed = ', seed, '.RData'))
    load(paste0('Genetic_RES/SCAN-DCF_delta = ', delta, '_hh = ', hh, '_seed = ', seed, '.RData'))
    load(paste0('Genetic_RES/KnockoffScreen_delta = ', delta, '_hh = ', hh, '_seed = ', seed, '.RData'))
    load(paste0('Genetic_RES/SCANG-STAAR_delta = ', delta, '_hh = ', hh, '_seed = ', seed, '.RData'))
    load(paste0('Genetic_RES/4S_delta = ', delta, '_hh = ', hh, '_seed = ', seed, '.RData'))
    load(paste0('Genetic_RES/DCF_delta = ', delta, '_hh = ', hh, '_seed = ', seed, '.RData'))
    
    BiRSDCF_T[hh, seed] = as.numeric(diffBSD, units = 'secs')
    #BiRSCL_T[hh, seed] = as.numeric(diffBSC, units = 'secs')
    QSCAN_T[hh, seed] = as.numeric(diffSCQ, units = 'secs')
    SCANDCF_T[hh, seed] = as.numeric(diffSCD, units = 'secs')
    KnockoffScreen_T[hh, seed] = as.numeric(diffKSD, units = 'secs')
    SSSS_T[hh, seed] = as.numeric(diffSSS, units = 'secs')
    STAAR_T[hh, seed] = as.numeric(diffSCT, units = 'secs')
    DCF_T[hh, seed] = as.numeric(diffDCF, units = 'secs')
  }
  
}

print(sum(apply(BiRSDCF_T, 1, mean)))
print(sum(apply(SCANDCF_T, 1, mean)))
print(sum(apply(DCF_T, 1, mean)))
print(sum(apply(QSCAN_T, 1, mean)))
print(sum(apply(KnockoffScreen_T, 1, mean)))
print(sum(apply(SSSS_T, 1, mean)))
print(sum(apply(STAAR_T, 1, mean)))
