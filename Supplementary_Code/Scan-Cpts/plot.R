rm(list = ls())

load('Review-Scan-Cpts/MES-BSD_Describe.RData')
load('Review-Scan-Cpts/MES-SBS_Describe.RData')
load('Review-Scan-Cpts/MES-WBS_Describe.RData')
load('Review-Scan-Cpts/MES-SCQ_Describe.RData')

powerBSD = matrix(powerBSD/100, 5, 4)
powerSBS = matrix(powerSBS/100, 5, 4)
powerWBS = matrix(powerWBS/100, 5, 4)

TRBSD = matrix(apply(TrateBSD, 2, mean), 5, 4)
TRSBS = matrix(apply(TrateSBS, 2, mean), 5, 4)
TRWBS = matrix(apply(TrateWBS, 2, mean), 5, 4)
TRSCQ = matrix(apply(TrateSCQ, 2, mean), 5, 4)
 
FRBSD = matrix(apply(FrateBSD, 2, mean), 5, 4)
FRSBS = matrix(apply(FrateSBS, 2, mean), 5, 4)
FRWBS = matrix(apply(FrateWBS, 2, mean), 5, 4)
FRSCQ = matrix(apply(FrateSCQ, 2, mean), 5, 4)
 
delta = c(0.20, 0.25, 0.30, 0.35, 0.40)
for (j in 1:4)
{
  nameTj = paste0('Figures/MES-TPR_beta', j, '.jpg')
  nameFj = paste0('Figures/MES-FDR_beta', j, '.jpg')
  
  # nameCTj = paste0('Review-Figures/MES-CP-TPR_beta', j, '.jpg')
  # nameCFj = paste0('Review-Figures/MES-CP-FDR_beta', j, '.jpg')
  
  mainTj = paste0('MES, beta = ', j, ' TPR')
  mainFj = paste0('MES, beta = ', j, ' FDR')
  
  jpeg(filename = nameTj, width = 1200, height = 900, quality = 100)
  plot(delta, TRBSD[ ,j], type = 'b', pch = 1, lty = 1, lwd = 4, col = 1, main = mainTj, 
       ylab = "", xlab = "delta", cex.lab = 3, cex.axis = 2.8, ylim = c(0, 1), cex.main = 4)
  lines(delta, TRSCQ[ ,j], type = "b", pch = 2, lty = 1, col = 2, lwd = 4)
  lines(delta, TRSBS[ ,j], type = "b", pch = 3, lty = 1, col = 3, lwd = 4)
  lines(delta, TRWBS[ ,j], type = "b", pch = 4, lty = 1, col = 4, lwd = 4)
  legend("topleft", c("BSD","SCQ", "SBS", "WBS"), pch = c(1, 2, 3, 4), lty = 1, box.col = "grey", cex = 3.5, col = c(1, 2, 3, 4), lwd = 4)
  dev.off()
  
  jpeg(filename = nameFj, width = 1200, height = 900, quality = 100)
  plot(delta, FRBSD[ ,j], type = 'b', pch = 1, lty = 1, lwd = 4, col = 1, main = mainFj, 
       ylab = "", xlab = "delta", cex.lab = 3, cex.axis = 2.8, ylim = c(0, 1), cex.main = 4)
  lines(delta, FRSCQ[ ,j], type = "b", pch = 2, lty = 1, col = 2, lwd = 4)
  lines(delta, FRSBS[ ,j], type = "b", pch = 3, lty = 1, col = 3, lwd = 4)
  lines(delta, FRWBS[ ,j], type = "b", pch = 4, lty = 1, col = 4, lwd = 4)
  legend("topleft", c("BSD","SCQ", "SBS", "WBS"), pch = c(1, 2, 3, 4), lty = 1, box.col = "grey", cex = 3.5, col = c(1, 2, 3, 4), lwd = 4)
  dev.off()
}

#############################################################################################################################################################################################

rm(list = ls())

load('Review-Scan-Cpts/MNS-BSD_Describe.RData')
load('Review-Scan-Cpts/MNS-SBS_Describe.RData')
load('Review-Scan-Cpts/MNS-WBS_Describe.RData')
load('Review-Scan-Cpts/MNS-SCQ_Describe.RData')

powerBSD = matrix(powerBSD/100, 5, 4)
powerSBS = matrix(powerSBS/100, 5, 4)
powerWBS = matrix(powerWBS/100, 5, 4)

TRBSD = matrix(apply(TrateBSD, 2, mean), 5, 4)
TRSBS = matrix(apply(TrateSBS, 2, mean), 5, 4)
TRWBS = matrix(apply(TrateWBS, 2, mean), 5, 4)
TRSCQ = matrix(apply(TrateSCQ, 2, mean), 5, 4)

FRBSD = matrix(apply(FrateBSD, 2, mean), 5, 4)
FRSBS = matrix(apply(FrateSBS, 2, mean), 5, 4)
FRWBS = matrix(apply(FrateWBS, 2, mean), 5, 4)
FRSCQ = matrix(apply(FrateSCQ, 2, mean), 5, 4)

delta = c(0.20, 0.25, 0.30, 0.35, 0.40)
for (j in 1:4)
{
  nameTj = paste0('Figures/MNS-TPR_beta', j, '.jpg')
  nameFj = paste0('Figures/MNS-FDR_beta', j, '.jpg')
  
  # nameCTj = paste0('Review-Figures/MNS-CP-TPR_beta', j, '.jpg')
  # nameCFj = paste0('Review-Figures/MNS-CP-FDR_beta', j, '.jpg')
  
  mainTj = paste0('MNS, beta = ', j, ' TPR')
  mainFj = paste0('MNS, beta = ', j, ' FDR')
  
  jpeg(filename = nameTj, width = 1200, height = 900, quality = 100)
  plot(delta, TRBSD[ ,j], type = 'b', pch = 1, lty = 1, lwd = 4, col = 1, main = mainTj, 
       ylab = "", xlab = "delta", cex.lab = 3, cex.axis = 2.8, ylim = c(0, 1), cex.main = 4)
  lines(delta, TRSCQ[ ,j], type = "b", pch = 2, lty = 1, col = 2, lwd = 4)
  lines(delta, TRSBS[ ,j], type = "b", pch = 3, lty = 1, col = 3, lwd = 4)
  lines(delta, TRWBS[ ,j], type = "b", pch = 4, lty = 1, col = 4, lwd = 4)
  legend("topleft", c("BSD","SCQ", "SBS", "WBS"), pch = c(1, 2, 3, 4), lty = 1, box.col = "grey", cex = 3.5, col = c(1, 2, 3, 4), lwd = 4)
  dev.off()
  
  jpeg(filename = nameFj, width = 1200, height = 900, quality = 100)
  plot(delta, FRBSD[ ,j], type = 'b', pch = 1, lty = 1, lwd = 4, col = 1, main = mainFj, 
       ylab = "", xlab = "delta", cex.lab = 3, cex.axis = 2.8, ylim = c(0, 1), cex.main = 4)
  lines(delta, FRSCQ[ ,j], type = "b", pch = 2, lty = 1, col = 2, lwd = 4)
  lines(delta, FRSBS[ ,j], type = "b", pch = 3, lty = 1, col = 3, lwd = 4)
  lines(delta, FRWBS[ ,j], type = "b", pch = 4, lty = 1, col = 4, lwd = 4)
  legend("topleft", c("BSD","SCQ", "SBS", "WBS"), pch = c(1, 2, 3, 4), lty = 1, box.col = "grey", cex = 3.5, col = c(1, 2, 3, 4), lwd = 4)
  dev.off()
}

#############################################################################################################################################################################################

rm(list = ls())

load('Review-Scan-Cpts/WES-BSD_Describe.RData')
load('Review-Scan-Cpts/WES-SBS_Describe.RData')
load('Review-Scan-Cpts/WES-WBS_Describe.RData')
load('Review-Scan-Cpts/WES-SCQ_Describe.RData')

powerBSD = matrix(powerBSD/100, 5, 4)
powerSBS = matrix(powerSBS/100, 5, 4)
powerWBS = matrix(powerWBS/100, 5, 4)

TRBSD = matrix(apply(TrateBSD, 2, mean), 5, 4)
TRSBS = matrix(apply(TrateSBS, 2, mean), 5, 4)
TRWBS = matrix(apply(TrateWBS, 2, mean), 5, 4)
TRSCQ = matrix(apply(TrateSCQ, 2, mean), 5, 4)

FRBSD = matrix(apply(FrateBSD, 2, mean), 5, 4)
FRSBS = matrix(apply(FrateSBS, 2, mean), 5, 4)
FRWBS = matrix(apply(FrateWBS, 2, mean), 5, 4)
FRSCQ = matrix(apply(FrateSCQ, 2, mean), 5, 4)

delta = c(0.20, 0.25, 0.30, 0.35, 0.40)
for (j in 1:4)
{
  nameTj = paste0('Figures/WES-TPR_beta', j, '.jpg')
  nameFj = paste0('Figures/WES-FDR_beta', j, '.jpg')
  
  # nameCTj = paste0('Review-Figures/WES-CP-TPR_beta', j, '.jpg')
  # nameCFj = paste0('Review-Figures/WES-CP-FDR_beta', j, '.jpg')
  
  mainTj = paste0('WES, beta = ', j, ' TPR')
  mainFj = paste0('WES, beta = ', j, ' FDR')
  
  jpeg(filename = nameTj, width = 1200, height = 900, quality = 100)
  plot(delta, TRBSD[ ,j], type = 'b', pch = 1, lty = 1, lwd = 4, col = 1, main = mainTj, 
       ylab = "", xlab = "delta", cex.lab = 3, cex.axis = 2.8, ylim = c(0, 1), cex.main = 4)
  lines(delta, TRSCQ[ ,j], type = "b", pch = 2, lty = 1, col = 2, lwd = 4)
  lines(delta, TRSBS[ ,j], type = "b", pch = 3, lty = 1, col = 3, lwd = 4)
  lines(delta, TRWBS[ ,j], type = "b", pch = 4, lty = 1, col = 4, lwd = 4)
  legend("topleft", c("BSD","SCQ", "SBS", "WBS"), pch = c(1, 2, 3, 4), lty = 1, box.col = "grey", cex = 3.5, col = c(1, 2, 3, 4), lwd = 4)
  dev.off()
  
  jpeg(filename = nameFj, width = 1200, height = 900, quality = 100)
  plot(delta, FRBSD[ ,j], type = 'b', pch = 1, lty = 1, lwd = 4, col = 1, main = mainFj, 
       ylab = "", xlab = "delta", cex.lab = 3, cex.axis = 2.8, ylim = c(0, 1), cex.main = 4)
  lines(delta, FRSCQ[ ,j], type = "b", pch = 2, lty = 1, col = 2, lwd = 4)
  lines(delta, FRSBS[ ,j], type = "b", pch = 3, lty = 1, col = 3, lwd = 4)
  lines(delta, FRWBS[ ,j], type = "b", pch = 4, lty = 1, col = 4, lwd = 4)
  legend("topleft", c("BSD","SCQ", "SBS", "WBS"), pch = c(1, 2, 3, 4), lty = 1, box.col = "grey", cex = 3.5, col = c(1, 2, 3, 4), lwd = 4)
  dev.off()
}


#############################################################################################################################################################################################

rm(list = ls())

load('Review-Scan-Cpts/WNS-BSD_Describe.RData')
load('Review-Scan-Cpts/WNS-SBS_Describe.RData')
load('Review-Scan-Cpts/WNS-WBS_Describe.RData')
load('Review-Scan-Cpts/WNS-SCQ_Describe.RData')

powerBSD = matrix(powerBSD/100, 5, 4)
powerSBS = matrix(powerSBS/100, 5, 4)
powerWBS = matrix(powerWBS/100, 5, 4)

TRBSD = matrix(apply(TrateBSD, 2, mean), 5, 4)
TRSBS = matrix(apply(TrateSBS, 2, mean), 5, 4)
TRWBS = matrix(apply(TrateWBS, 2, mean), 5, 4)
TRSCQ = matrix(apply(TrateSCQ, 2, mean), 5, 4)

FRBSD = matrix(apply(FrateBSD, 2, mean), 5, 4)
FRSBS = matrix(apply(FrateSBS, 2, mean), 5, 4)
FRWBS = matrix(apply(FrateWBS, 2, mean), 5, 4)
FRSCQ = matrix(apply(FrateSCQ, 2, mean), 5, 4)

delta = c(0.20, 0.25, 0.30, 0.35, 0.40)
for (j in 1:4)
{
  nameTj = paste0('Figures/WNS-TPR_beta', j, '.jpg')
  nameFj = paste0('Figures/WNS-FDR_beta', j, '.jpg')
  
  # nameCTj = paste0('Review-Figures/WNS-CP-TPR_beta', j, '.jpg')
  # nameCFj = paste0('Review-Figures/WNS-CP-FDR_beta', j, '.jpg')
  
  mainTj = paste0('WNS, beta = ', j, ' TPR')
  mainFj = paste0('WNS, beta = ', j, ' FDR')
  
  jpeg(filename = nameTj, width = 1200, height = 900, quality = 100)
  plot(delta, TRBSD[ ,j], type = 'b', pch = 1, lty = 1, lwd = 4, col = 1, main = mainTj, 
       ylab = "", xlab = "delta", cex.lab = 3, cex.axis = 2.8, ylim = c(0, 1), cex.main = 4)
  lines(delta, TRSCQ[ ,j], type = "b", pch = 2, lty = 1, col = 2, lwd = 4)
  lines(delta, TRSBS[ ,j], type = "b", pch = 3, lty = 1, col = 3, lwd = 4)
  lines(delta, TRWBS[ ,j], type = "b", pch = 4, lty = 1, col = 4, lwd = 4)
  legend("topleft", c("BSD","SCQ", "SBS", "WBS"), pch = c(1, 2, 3, 4), lty = 1, box.col = "grey", cex = 3.5, col = c(1, 2, 3, 4), lwd = 4)
  dev.off()
  
  jpeg(filename = nameFj, width = 1200, height = 900, quality = 100)
  plot(delta, FRBSD[ ,j], type = 'b', pch = 1, lty = 1, lwd = 4, col = 1, main = mainFj, 
       ylab = "", xlab = "delta", cex.lab = 3, cex.axis = 2.8, ylim = c(0, 1), cex.main = 4)
  lines(delta, FRSCQ[ ,j], type = "b", pch = 2, lty = 1, col = 2, lwd = 4)
  lines(delta, FRSBS[ ,j], type = "b", pch = 3, lty = 1, col = 3, lwd = 4)
  lines(delta, FRWBS[ ,j], type = "b", pch = 4, lty = 1, col = 4, lwd = 4)
  legend("topleft", c("BSD","SCQ", "SBS", "WBS"), pch = c(1, 2, 3, 4), lty = 1, box.col = "grey", cex = 3.5, col = c(1, 2, 3, 4), lwd = 4)
  dev.off()
}
