rm(list = ls())

load('Normal_ESBA/BSD_Describe.RData')
load('Normal_ESBA/BSC_Describe.RData')
load('Normal_ESBA/SCQ_Describe.RData')
load('Normal_ESBA/KSD_Describe.RData')

powerBSD = matrix(powerBSD/500, 5, 4)
powerBSC = matrix(powerBSC/500, 5, 4)
powerSCQ = matrix(powerSCQ/500, 5, 4)
powerSCM = matrix(powerSCM/100, 5, 4)

TRBSD = matrix(apply(TrateBSD, 2, mean), 5, 4)
TRBSC = matrix(apply(TrateBSC, 2, mean), 5, 4)
TRSCQ = matrix(apply(TrateSCQ, 2, mean), 5, 4)
TRKSD = matrix(apply(TrateKSD, 2, mean), 5, 4)

FRBSD = matrix(apply(FrateBSD, 2, mean), 5, 4)
FRBSC = matrix(apply(FrateBSC, 2, mean), 5, 4)
FRSCQ = matrix(apply(FrateSCQ, 2, mean), 5, 4)
FRKSD = matrix(apply(FrateKSD, 2, mean), 5, 4)

DTBSD = matrix(apply(DifTiBSD, 2, mean), 5, 4)
DTBSC = matrix(apply(DifTiBSC, 2, mean), 5, 4)
DTSCQ = matrix(apply(DifTiSCQ, 2, mean), 5, 4)
DTKSD = matrix(apply(DifTiKSD, 2, mean), 5, 4)

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

SCQ = NULL
for (j in 1:5)
{
  SCQ = cbind(SCQ, powerSCQ[j, ], TRSCQ[j, ], FRSCQ[j, ])
}
SCQ = as.data.frame(SCQ)

KSD = NULL
for (j in 1:5)
{
  KSD = cbind(KSD, powerKSD[j, ], TRKSD[j, ], FRKSD[j, ])
}
KSD = as.data.frame(KSD)


delta = c(0.20, 0.25, 0.30, 0.35, 0.40)
for (j in 1:4)
{
  
  namePj = paste0('Figures/Normal_ESBA_power_beta', j, '.jpg')
  nameDj = paste0('Figures/Normal_ESBA_Diftime_beta', j, '.jpg')
  nameTj = paste0('Figures/Normal_ESBA_TDR_beta', j, '.jpg')
  nameFj = paste0('Figures/Normal_ESBA_FDR_beta', j, '.jpg')
  
  mainPj = paste0('Equal covariance, beta = ', j, ' power')
  mainTj = paste0('Equal covariance, beta = ', j, ' TPR')
  mainFj = paste0('Equal covariance, beta = ', j, ' FDR')
  
  jpeg(filename = namePj, width = 1200, height = 900, quality = 120)
  plot(delta, powerBSD[ ,j], type = 'b', pch = 1, lty = 1, lwd = 4, col = 1, main = mainPj, 
       ylab = "", xlab = "delta", cex.lab = 3, cex.axis = 2.8, ylim = c(0, 1), cex.main = 4)
  lines(delta, powerBSC[ ,j], type = "b", pch = 2, lty = 1, col = 2, lwd = 4)
  lines(delta, powerSCQ[ ,j], type = "b", pch = 3, lty = 1, col = 3, lwd = 4)
  lines(delta, powerKSD[ ,j], type = "b", pch = 4, lty = 1, col = 4, lwd = 4)
  legend("bottomright", c("BSD", "BSC", "SCQ", "KSD"), pch = c(1, 2, 3, 4), lty = 1, box.col = "grey", cex = 3.5, col = c(1, 2, 3, 4), lwd = 4)
  dev.off()
  
  jpeg(filename = nameTj, width = 1200, height = 900, quality = 120)
  plot(delta, TRBSD[ ,j], type = 'b', pch = 1, lty = 1, lwd = 4, col = 1, main = mainTj, 
       ylab = "", xlab = "delta", cex.lab = 3, cex.axis = 2.8, ylim = c(0, 1), cex.main = 4)
  lines(delta, TRBSC[ ,j], type = "b", pch = 2, lty = 1, col = 2, lwd = 4)
  lines(delta, TRSCQ[ ,j], type = "b", pch = 3, lty = 1, col = 3, lwd = 4)
  lines(delta, TRKSD[ ,j], type = "b", pch = 4, lty = 1, col = 4, lwd = 4)
  legend("topleft", c("BSD", "BSC", "SCQ", "KSD"), pch = c(1, 2, 3, 4), lty = 1, box.col = "grey", cex = 3.5, col = c(1, 2, 3, 4), lwd = 4)
  dev.off()
  
  jpeg(filename = nameFj, width = 1200, height = 900, quality = 120)
  plot(delta, FRBSD[ ,j], type = 'b', pch = 1, lty = 1, lwd = 4, col = 1, main = mainFj, 
       ylab = "", xlab = "delta", cex.lab = 3, cex.axis = 2.8, ylim = c(0, 0.3), cex.main = 4)
  lines(delta, FRBSC[ ,j], type = "b", pch = 2, lty = 1, col = 2, lwd = 4)
  lines(delta, FRSCQ[ ,j], type = "b", pch = 3, lty = 1, col = 3, lwd = 4)
  lines(delta, FRKSD[ ,j], type = "b", pch = 4, lty = 1, col = 4, lwd = 4)
  legend("topleft", c("BSD", "BSC", "SCQ", "KSD"), pch = c(1, 2, 3, 4), lty = 1, box.col = "grey", cex = 3.5, col = c(1, 2, 3, 4), lwd = 4)
  dev.off()
}

######################################################################################################################################
rm(list = ls())

load('Normal_NSBA/BSD_Describe.RData')
load('Normal_NSBA/BSC_Describe.RData')
load('Normal_NSBA/SCQ_Describe.RData')
load('Normal_NSBA/KSD_Describe.RData')

powerBSD = matrix(powerBSD/500, 5, 4)
powerBSC = matrix(powerBSC/500, 5, 4)
powerSCQ = matrix(powerSCQ/500, 5, 4)
powerSCM = matrix(powerSCM/100, 5, 4)

TRBSD = matrix(apply(TrateBSD, 2, mean), 5, 4)
TRBSC = matrix(apply(TrateBSC, 2, mean), 5, 4)
TRSCQ = matrix(apply(TrateSCQ, 2, mean), 5, 4)
TRKSD = matrix(apply(TrateKSD, 2, mean), 5, 4)

FRBSD = matrix(apply(FrateBSD, 2, mean), 5, 4)
FRBSC = matrix(apply(FrateBSC, 2, mean), 5, 4)
FRSCQ = matrix(apply(FrateSCQ, 2, mean), 5, 4)
FRKSD = matrix(apply(FrateKSD, 2, mean), 5, 4)

DTBSD = matrix(apply(DifTiBSD, 2, mean), 5, 4)
DTBSC = matrix(apply(DifTiBSC, 2, mean), 5, 4)
DTSCQ = matrix(apply(DifTiSCQ, 2, mean), 5, 4)
DTKSD = matrix(apply(DifTiKSD, 2, mean), 5, 4)

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

SCQ = NULL
for (j in 1:5)
{
  SCQ = cbind(SCQ, powerSCQ[j, ], TRSCQ[j, ], FRSCQ[j, ])
}
SCQ = as.data.frame(SCQ)

KSD = NULL
for (j in 1:5)
{
  KSD = cbind(KSD, powerKSD[j, ], TRKSD[j, ], FRKSD[j, ])
}
KSD = as.data.frame(KSD)


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
  lines(delta, powerKSD[ ,j], type = "b", pch = 4, lty = 1, col = 4, lwd = 4)
  legend("bottomright", c("BSD", "BSC", "SCQ", "KSD"), pch = c(1, 2, 3, 4), lty = 1, box.col = "grey", cex = 3.5, col = c(1, 2, 3, 4), lwd = 4)
  dev.off()
  
  jpeg(filename = nameTj, width = 1200, height = 900, quality = 120)
  plot(delta, TRBSD[ ,j], type = 'b', pch = 1, lty = 1, lwd = 4, col = 1, main = mainTj, 
       ylab = "", xlab = "delta", cex.lab = 3, cex.axis = 2.8, ylim = c(0, 1), cex.main = 4)
  lines(delta, TRBSC[ ,j], type = "b", pch = 2, lty = 1, col = 2, lwd = 4)
  lines(delta, TRSCQ[ ,j], type = "b", pch = 3, lty = 1, col = 3, lwd = 4)
  lines(delta, TRKSD[ ,j], type = "b", pch = 4, lty = 1, col = 4, lwd = 4)
  legend("topleft", c("BSD", "BSC", "SCQ", "KSD"), pch = c(1, 2, 3, 4), lty = 1, box.col = "grey", cex = 3.5, col = c(1, 2, 3, 4), lwd = 4)
  dev.off()
  
  jpeg(filename = nameFj, width = 1200, height = 900, quality = 120)
  plot(delta, FRBSD[ ,j], type = 'b', pch = 1, lty = 1, lwd = 4, col = 1, main = mainFj, 
       ylab = "", xlab = "delta", cex.lab = 3, cex.axis = 2.8, ylim = c(0, 0.3), cex.main = 4)
  lines(delta, FRBSC[ ,j], type = "b", pch = 2, lty = 1, col = 2, lwd = 4)
  lines(delta, FRSCQ[ ,j], type = "b", pch = 3, lty = 1, col = 3, lwd = 4)
  lines(delta, FRKSD[ ,j], type = "b", pch = 4, lty = 1, col = 4, lwd = 4)
  legend("topleft", c("BSD", "BSC", "SCQ", "KSD"), pch = c(1, 2, 3, 4), lty = 1, box.col = "grey", cex = 3.5, col = c(1, 2, 3, 4), lwd = 4)
  dev.off()
}