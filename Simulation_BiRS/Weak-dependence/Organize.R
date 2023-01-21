rm(list = ls())

load('Normal_ESBA/BSD_Describe.RData')
load('Normal_ESBA/BSC_Describe.RData')
#load('Normal_ESBA/SCD_Describe.RData')
#load('Normal_ESBA/SCC_Describe.RData')
load('Normal_ESBA/SCQ_Describe.RData')
#load('Normal_ESBA/SCM_Describe.RData')

powerBSD = matrix(powerBSD/1000, 5, 4)
powerBSC = matrix(powerBSC/1000, 5, 4)
#powerSCD = matrix(powerSCD/1000, 5, 4)
#powerSCC = matrix(powerSCC/1000, 5, 4)
powerSCQ = matrix(powerSCQ/1000, 5, 4)
#powerSCM = matrix(powerSCM/100, 5, 4)

TRBSD = matrix(apply(TrateBSD, 2, mean), 5, 4)
TRBSC = matrix(apply(TrateBSC, 2, mean), 5, 4)
#TRSCD = matrix(apply(TrateSCD, 2, mean), 5, 4)
#TRSCC = matrix(apply(TrateSCC, 2, mean), 5, 4)
TRSCQ = matrix(apply(TrateSCQ, 2, mean), 5, 4)
#TRSCM = matrix(apply(TrateSCM, 2, mean), 5, 4)

FRBSD = matrix(apply(FrateBSD, 2, mean), 5, 4)
FRBSC = matrix(apply(FrateBSC, 2, mean), 5, 4)
#FRSCD = matrix(apply(FrateSCD, 2, mean), 5, 4)
#FRSCC = matrix(apply(FrateSCC, 2, mean), 5, 4)
FRSCQ = matrix(apply(FrateSCQ, 2, mean), 5, 4)
#FRSCM = matrix(apply(FrateSCM, 2, mean), 5, 4)

DTBSD = matrix(apply(DifTiBSD, 2, mean), 5, 4)
DTBSC = matrix(apply(DifTiBSC, 2, mean), 5, 4)
#DTSCD = matrix(apply(DifTiSCD, 2, mean), 5, 4)
#DTSCC = matrix(apply(DifTiSCC, 2, mean), 5, 4)
DTSCQ = matrix(apply(DifTiSCQ, 2, mean), 5, 4)
#DTSCM = matrix(apply(DifTiSCM, 2, mean), 5, 4)

BSD = NULL
for (j in 1:5)
{
  BSD = cbind(BSD, powerBSD[j, ], TRBSD[j, ], FRBSD[j, ])
}
BSD = as.data.frame(BSD)

BSC = NULL
for (j in 1:5)
{
  BSC = cbind(BSC, powerBSC[j, ], TRBSC[j, ], FRBSC[j, ])
}
BSC = as.data.frame(BSC)

# SCD = NULL
# for (j in 1:5)
# {
#   SCD = cbind(SCD, powerSCD[j, ], TRSCD[j, ], FRSCD[j, ])
# }
# SCD = as.data.frame(SCD)
# 
# SCC = NULL
# for (j in 1:5)
# {
#   SCC = cbind(SCC, powerSCC[j, ], TRSCC[j, ], FRSCC[j, ])
# }
# SCC = as.data.frame(SCC)

SCQ = NULL
for (j in 1:5)
{
  SCQ = cbind(SCQ, powerSCQ[j, ], TRSCQ[j, ], FRSCQ[j, ])
}
SCQ = as.data.frame(SCQ)

# DT = NULL
# for (j in 1:4)
# {
#   DT = cbind(DT, DTBSD[ ,j], DTBSC[ ,j], DTSCD[ ,j], DTSCC[ ,j], 
#              DTSCQ[ ,j])
# }
# DT = as.data.frame(DT)


delta = c(0.20, 0.25, 0.30, 0.35, 0.40)
for (j in 1:4)
{
  
  namePj = paste0('Figures/Normal_ESBA_power_beta', j, '.jpg')
  nameDj = paste0('Figures/Normal_ESBA_Diftime_beta', j, '.jpg')
  nameTj = paste0('Figures/Normal_ESBA_TDR_beta', j, '.jpg')
  nameFj = paste0('Figures/Normal_ESBA_FDR_beta', j, '.jpg')
  
  mainPj = paste0('Unequal covariance, beta = ', j, ' power')
  mainTj = paste0('Unequal covariance, beta = ', j, ' TPR')
  mainFj = paste0('Unequal covariance, beta = ', j, ' FDR')
  
  jpeg(filename = namePj, width = 1200, height = 900, quality = 120)
  plot(delta, powerBSD[ ,j], type = 'b', pch = 1, lty = 1, lwd = 4, col = 1, main = mainPj, 
       ylab = "", xlab = "delta", cex.lab = 3, cex.axis = 2.8, ylim = c(0, 1), cex.main = 4)
  lines(delta, powerBSC[ ,j], type = "b", pch = 2, lty = 1, col = 2, lwd = 4)
  lines(delta, powerSCQ[ ,j], type = "b", pch = 3, lty = 1, col = 3, lwd = 4)
  legend("bottomright", c("BSD", "BSC", "SCQ"), pch = c(1, 2, 3), lty = 1, box.col = "grey", cex = 4, col = c(1, 2, 3), lwd = 4)
  dev.off()
  
  jpeg(filename = nameTj, width = 1200, height = 900, quality = 120)
  plot(delta, TRBSD[ ,j], type = 'b', pch = 1, lty = 1, lwd = 4, col = 1, main = mainTj, 
       ylab = "", xlab = "delta", cex.lab = 3, cex.axis = 2.8, ylim = c(0, 1), cex.main = 4)
  lines(delta, TRBSC[ ,j], type = "b", pch = 2, lty = 1, col = 2, lwd = 4)
  lines(delta, TRSCQ[ ,j], type = "b", pch = 3, lty = 1, col = 3, lwd = 4)
  # lines(delta, TRSCD[ ,j], type = "b", pch = 4, lty = 1, col = 4, lwd = 2)
  # lines(delta, TRSCC[ ,j], type = "b", pch = 20, lty = 1, col = 6, lwd = 2)
  # legend("bottomright", c("BSD", "BSC", "SCQ", "SCD", "SCC"), pch = c(1, 2, 3, 4, 20), lty = 1, box.col = "grey", cex = 2.5, 
  #        col = c(1, 2, 3, 4, 6), lwd = 2)
  legend("topleft", c("BSD", "BSC", "SCQ"), pch = c(1, 2, 3), lty = 1, box.col = "grey", cex = 4, col = c(1, 2, 3), lwd = 4)
  dev.off()
  
  jpeg(filename = nameFj, width = 1200, height = 900, quality = 120)
  plot(delta, FRBSD[ ,j], type = 'b', pch = 1, lty = 1, lwd = 4, col = 1, main = mainFj, 
       ylab = "", xlab = "delta", cex.lab = 3, cex.axis = 2.8, ylim = c(0, 0.3), cex.main = 4)
  lines(delta, FRBSC[ ,j], type = "b", pch = 2, lty = 1, col = 2, lwd = 4)
  lines(delta, FRSCQ[ ,j], type = "b", pch = 3, lty = 1, col = 3, lwd = 4)
  #lines(delta, rep(0.05, length(delta)), type = 'l', lty = 2, col = 4)
  #lines(delta, rep(0.06, length(delta)), type = 'l', lty = 2, col = 4)
  # lines(delta, FRSCD[ ,j], type = "b", pch = 4, lty = 1, col = 4, lwd = 2)
  # lines(delta, FRSCC[ ,j], type = "b", pch = 20, lty = 1, col = 6, lwd = 2)
  # legend("topleft", c("BSD", "BSC", "SCQ", "SCD", "SCC"), pch = c(1, 2, 3, 4, 20), lty = 1, box.col = "grey", cex = 2.5, 
  #        col = c(1, 2, 3, 4, 6), lwd = 2)
  legend("topleft", c("BSD", "BSC", "SCQ"), pch = c(1, 2, 3), lty = 1, box.col = "grey", cex = 4, col = c(1, 2, 3), lwd = 4)
  dev.off()
}

######################################################################################################################################
rm(list = ls())

load('Normal_NSBA/BSD_Describe.RData')
load('Normal_NSBA/BSC_Describe.RData')
#load('Normal_NSBA/SCD_Describe.RData')
#load('Normal_NSBA/SCC_Describe.RData')
load('Normal_NSBA/SCQ_Describe.RData')
#load('Normal_NSBA/SCM_Describe.RData')

powerBSD = matrix(powerBSD/1000, 5, 4)
powerBSC = matrix(powerBSC/1000, 5, 4)
#powerSCD = matrix(powerSCD/1000, 5, 4)
#powerSCC = matrix(powerSCC/1000, 5, 4)
powerSCQ = matrix(powerSCQ/1000, 5, 4)
#powerSCM = matrix(powerSCM/100, 5, 4)

TRBSD = matrix(apply(TrateBSD, 2, mean), 5, 4)
TRBSC = matrix(apply(TrateBSC, 2, mean), 5, 4)
#TRSCD = matrix(apply(TrateSCD, 2, mean), 5, 4)
#TRSCC = matrix(apply(TrateSCC, 2, mean), 5, 4)
TRSCQ = matrix(apply(TrateSCQ, 2, mean), 5, 4)
#TRSCM = matrix(apply(TrateSCM, 2, mean), 5, 4)

FRBSD = matrix(apply(FrateBSD, 2, mean), 5, 4)
FRBSC = matrix(apply(FrateBSC, 2, mean), 5, 4)
#FRSCD = matrix(apply(FrateSCD, 2, mean), 5, 4)
#FRSCC = matrix(apply(FrateSCC, 2, mean), 5, 4)
FRSCQ = matrix(apply(FrateSCQ, 2, mean), 5, 4)
#FRSCM = matrix(apply(FrateSCM, 2, mean), 5, 4)

DTBSD = matrix(apply(DifTiBSD, 2, mean), 5, 4)
DTBSC = matrix(apply(DifTiBSC, 2, mean), 5, 4)
#DTSCD = matrix(apply(DifTiSCD, 2, mean), 5, 4)
#DTSCC = matrix(apply(DifTiSCC, 2, mean), 5, 4)
DTSCQ = matrix(apply(DifTiSCQ, 2, mean), 5, 4)
#DTSCM = matrix(apply(DifTiSCM, 2, mean), 5, 4)

BSD = NULL
for (j in 1:5)
{
  BSD = cbind(BSD, powerBSD[j, ], TRBSD[j, ], FRBSD[j, ])
}
BSD = as.data.frame(BSD)

BSC = NULL
for (j in 1:5)
{
  BSC = cbind(BSC, powerBSC[j, ], TRBSC[j, ], FRBSC[j, ])
}
BSC = as.data.frame(BSC)

# SCD = NULL
# for (j in 1:5)
# {
#   SCD = cbind(SCD, powerSCD[j, ], TRSCD[j, ], FRSCD[j, ])
# }
# SCD = as.data.frame(SCD)
# 
# SCC = NULL
# for (j in 1:5)
# {
#   SCC = cbind(SCC, powerSCC[j, ], TRSCC[j, ], FRSCC[j, ])
# }
# SCC = as.data.frame(SCC)

SCQ = NULL
for (j in 1:5)
{
  SCQ = cbind(SCQ, powerSCQ[j, ], TRSCQ[j, ], FRSCQ[j, ])
}
SCQ = as.data.frame(SCQ)

# DT = NULL
# for (j in 1:4)
# {
#   DT = cbind(DT, DTBSD[ ,j], DTBSC[ ,j], DTSCD[ ,j], DTSCC[ ,j], 
#              DTSCQ[ ,j])
# }
# DT = as.data.frame(DT)


delta = c(0.20, 0.25, 0.30, 0.35, 0.40)
for (j in 1:4)
{
  
  namePj = paste0('Figures/Normal_NSBA_power_beta', j, '.jpg')
  nameDj = paste0('Figures/Normal_NSBA_Diftime_beta', j, '.jpg')
  nameTj = paste0('Figures/Normal_NSBA_TDR_beta', j, '.jpg')
  nameFj = paste0('Figures/Normal_NSBA_FDR_beta', j, '.jpg')
  
  mainPj = paste0('Unequal covariance, beta = ', j, ' power')
  mainTj = paste0('Unequal covariance, beta = ', j, ' TPR')
  mainFj = paste0('Unequal covariance, beta = ', j, ' FDR')
  
  jpeg(filename = namePj, width = 1200, height = 900, quality = 120)
  plot(delta, powerBSD[ ,j], type = 'b', pch = 1, lty = 1, lwd = 4, col = 1, main = mainPj, 
       ylab = "", xlab = "delta", cex.lab = 3, cex.axis = 2.8, ylim = c(0, 1), cex.main = 4)
  lines(delta, powerBSC[ ,j], type = "b", pch = 2, lty = 1, col = 2, lwd = 4)
  lines(delta, powerSCQ[ ,j], type = "b", pch = 3, lty = 1, col = 3, lwd = 4)
  legend("bottomright", c("BSD", "BSC", "SCQ"), pch = c(1, 2, 3), lty = 1, box.col = "grey", cex = 4, col = c(1, 2, 3), lwd = 4)
  dev.off()
  
  jpeg(filename = nameTj, width = 1200, height = 900, quality = 120)
  plot(delta, TRBSD[ ,j], type = 'b', pch = 1, lty = 1, lwd = 4, col = 1, main = mainTj, 
       ylab = "", xlab = "delta", cex.lab = 3, cex.axis = 2.8, ylim = c(0, 1), cex.main = 4)
  lines(delta, TRBSC[ ,j], type = "b", pch = 2, lty = 1, col = 2, lwd = 4)
  lines(delta, TRSCQ[ ,j], type = "b", pch = 3, lty = 1, col = 3, lwd = 4)
  # lines(delta, TRSCD[ ,j], type = "b", pch = 4, lty = 1, col = 4, lwd = 2)
  # lines(delta, TRSCC[ ,j], type = "b", pch = 20, lty = 1, col = 6, lwd = 2)
  # legend("bottomright", c("BSD", "BSC", "SCQ", "SCD", "SCC"), pch = c(1, 2, 3, 4, 20), lty = 1, box.col = "grey", cex = 2.5, 
  #        col = c(1, 2, 3, 4, 6), lwd = 2)
  legend("topleft", c("BSD", "BSC", "SCQ"), pch = c(1, 2, 3), lty = 1, box.col = "grey", cex = 4, col = c(1, 2, 3), lwd = 4)
  dev.off()
  
  jpeg(filename = nameFj, width = 1200, height = 900, quality = 120)
  plot(delta, FRBSD[ ,j], type = 'b', pch = 1, lty = 1, lwd = 4, col = 1, main = mainFj, 
       ylab = "", xlab = "delta", cex.lab = 3, cex.axis = 2.8, ylim = c(0, 0.3), cex.main = 4)
  lines(delta, FRBSC[ ,j], type = "b", pch = 2, lty = 1, col = 2, lwd = 4)
  lines(delta, FRSCQ[ ,j], type = "b", pch = 3, lty = 1, col = 3, lwd = 4)
  #lines(delta, rep(0.05, length(delta)), type = 'l', lty = 2, col = 4)
  #lines(delta, rep(0.06, length(delta)), type = 'l', lty = 2, col = 4)
  # lines(delta, FRSCD[ ,j], type = "b", pch = 4, lty = 1, col = 4, lwd = 2)
  # lines(delta, FRSCC[ ,j], type = "b", pch = 20, lty = 1, col = 6, lwd = 2)
  # legend("topleft", c("BSD", "BSC", "SCQ", "SCD", "SCC"), pch = c(1, 2, 3, 4, 20), lty = 1, box.col = "grey", cex = 2.5, 
  #        col = c(1, 2, 3, 4, 6), lwd = 2)
  legend("topleft", c("BSD", "BSC", "SCQ"), pch = c(1, 2, 3), lty = 1, box.col = "grey", cex = 4, col = c(1, 2, 3), lwd = 4)
  dev.off()
}