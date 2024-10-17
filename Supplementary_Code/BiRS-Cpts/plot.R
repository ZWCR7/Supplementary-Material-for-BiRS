rm(list = ls())

load('Review-BiRS-cpts/MES-BSD_Describe.RData')
load('Review-BiRS-cpts/MES-WBB_Describe.RData')

powerBSD = matrix(powerBSD/100, 5, 4)
powerWBB = matrix(powerWBB/100, 5, 4)

TRBSD = matrix(apply(TrateBSD, 2, mean), 5, 4)
TRWBB = matrix(apply(TrateWBB, 2, mean), 5, 4)

FRBSD = matrix(apply(FrateBSD, 2, mean), 5, 4)
FRWBB = matrix(apply(FrateWBB, 2, mean), 5, 4)

DTBSD = matrix(apply(DifTiBSD, 2, mean), 5, 4)
DTWBB = matrix(apply(DifTiWBB, 2, mean), 5, 4)


delta = c(0.20, 0.25, 0.30, 0.35, 0.40)
for (j in 1:4)
{
  nameTj = paste0('Figures/MES-TPR_beta', j, '.jpg')
  nameFj = paste0('Figures/MES-FDR_beta', j, '.jpg')
  
  mainTj = paste0('MES, beta = ', j, ' TPR')
  mainFj = paste0('MES, beta = ', j, ' FDR')
  
  jpeg(filename = nameTj, width = 1200, height = 900, quality = 100)
  plot(delta, TRBSD[ ,j], type = 'b', pch = 1, lty = 1, lwd = 4, col = 1, main = mainTj, 
       ylab = "", xlab = "delta", cex.lab = 3, cex.axis = 2.8, ylim = c(0, 1), cex.main = 4)
  lines(delta, TRWBB[ ,j], type = "b", pch = 2, lty = 1, col = 2, lwd = 4)
  legend("right", c("BSD", "WBB"), pch = c(1, 2), lty = 1, box.col = "grey", cex = 3.5, col = c(1, 2), lwd = 4)
  dev.off()
  
  jpeg(filename = nameFj, width = 1200, height = 900, quality = 100)
  plot(delta, FRBSD[ ,j], type = 'b', pch = 1, lty = 1, lwd = 4, col = 1, main = mainFj, 
       ylab = "", xlab = "delta", cex.lab = 3, cex.axis = 2.8, ylim = c(0, 0.5), cex.main = 4)
  lines(delta, FRWBB[ ,j], type = "b", pch = 2, lty = 1, col = 2, lwd = 4)
  legend("topright", c("BSD", "WBB"), pch = c(1, 2), lty = 1, box.col = "grey", cex = 3.5, col = c(1, 2), lwd = 4)
  dev.off()
}

#############################################################################################################################################################################################

rm(list = ls())

load('Review-BiRS-cpts/MNS-BSD_Describe.RData')
load('Review-BiRS-cpts/MNS-WBB_Describe.RData')

powerBSD = matrix(powerBSD/100, 5, 4)
powerWBB = matrix(powerWBB/100, 5, 4)

TRBSD = matrix(apply(TrateBSD, 2, mean), 5, 4)
TRWBB = matrix(apply(TrateWBB, 2, mean), 5, 4)

FRBSD = matrix(apply(FrateBSD, 2, mean), 5, 4)
FRWBB = matrix(apply(FrateWBB, 2, mean), 5, 4)

DTBSD = matrix(apply(DifTiBSD, 2, mean), 5, 4)
DTWBB = matrix(apply(DifTiWBB, 2, mean), 5, 4)


delta = c(0.20, 0.25, 0.30, 0.35, 0.40)
for (j in 1:4)
{
  nameTj = paste0('Figures/MNS-TPR_beta', j, '.jpg')
  nameFj = paste0('Figures/MNS-FDR_beta', j, '.jpg')
  
  mainTj = paste0('MNS, beta = ', j, ' TPR')
  mainFj = paste0('MNS, beta = ', j, ' FDR')
  
  jpeg(filename = nameTj, width = 1200, height = 900, quality = 100)
  plot(delta, TRBSD[ ,j], type = 'b', pch = 1, lty = 1, lwd = 4, col = 1, main = mainTj, 
       ylab = "", xlab = "delta", cex.lab = 3, cex.axis = 2.8, ylim = c(0, 1), cex.main = 4)
  lines(delta, TRWBB[ ,j], type = "b", pch = 2, lty = 1, col = 2, lwd = 4)
  legend("right", c("BSD", "WBB"), pch = c(1, 2), lty = 1, box.col = "grey", cex = 3.5, col = c(1, 2), lwd = 4)
  dev.off()
  
  jpeg(filename = nameFj, width = 1200, height = 900, quality = 100)
  plot(delta, FRBSD[ ,j], type = 'b', pch = 1, lty = 1, lwd = 4, col = 1, main = mainFj, 
       ylab = "", xlab = "delta", cex.lab = 3, cex.axis = 2.8, ylim = c(0, 0.5), cex.main = 4)
  lines(delta, FRWBB[ ,j], type = "b", pch = 2, lty = 1, col = 2, lwd = 4)
  legend("topright", c("BSD", "WBB"), pch = c(1, 2), lty = 1, box.col = "grey", cex = 3.5, col = c(1, 2), lwd = 4)
  dev.off()
}

#############################################################################################################################################################################################

rm(list = ls())

load('Review-BiRS-cpts/WES-BSD_Describe.RData')
load('Review-BiRS-cpts/WES-WBB_Describe.RData')

powerBSD = matrix(powerBSD/100, 5, 4)
powerWBB = matrix(powerWBB/100, 5, 4)

TRBSD = matrix(apply(TrateBSD, 2, mean), 5, 4)
TRWBB = matrix(apply(TrateWBB, 2, mean), 5, 4)

FRBSD = matrix(apply(FrateBSD, 2, mean), 5, 4)
FRWBB = matrix(apply(FrateWBB, 2, mean), 5, 4)

DTBSD = matrix(apply(DifTiBSD, 2, mean), 5, 4)
DTWBB = matrix(apply(DifTiWBB, 2, mean), 5, 4)


delta = c(0.20, 0.25, 0.30, 0.35, 0.40)
for (j in 1:4)
{
  nameTj = paste0('Figures/WES-TPR_beta', j, '.jpg')
  nameFj = paste0('Figures/WES-FDR_beta', j, '.jpg')
  
  mainTj = paste0('WES, beta = ', j, ' TPR')
  mainFj = paste0('WES, beta = ', j, ' FDR')
  
  jpeg(filename = nameTj, width = 1200, height = 900, quality = 100)
  plot(delta, TRBSD[ ,j], type = 'b', pch = 1, lty = 1, lwd = 4, col = 1, main = mainTj, 
       ylab = "", xlab = "delta", cex.lab = 3, cex.axis = 2.8, ylim = c(0, 1), cex.main = 4)
  lines(delta, TRWBB[ ,j], type = "b", pch = 2, lty = 1, col = 2, lwd = 4)
  legend("right", c("BSD", "WBB"), pch = c(1, 2), lty = 1, box.col = "grey", cex = 3.5, col = c(1, 2), lwd = 4)
  dev.off()
  
  jpeg(filename = nameFj, width = 1200, height = 900, quality = 100)
  plot(delta, FRBSD[ ,j], type = 'b', pch = 1, lty = 1, lwd = 4, col = 1, main = mainFj, 
       ylab = "", xlab = "delta", cex.lab = 3, cex.axis = 2.8, ylim = c(0, 0.5), cex.main = 4)
  lines(delta, FRWBB[ ,j], type = "b", pch = 2, lty = 1, col = 2, lwd = 4)
  legend("topright", c("BSD", "WBB"), pch = c(1, 2), lty = 1, box.col = "grey", cex = 3.5, col = c(1, 2), lwd = 4)
  dev.off()
}

#############################################################################################################################################################################################

rm(list = ls())

load('Review-BiRS-cpts/WNS-BSD_Describe.RData')
load('Review-BiRS-cpts/WNS-WBB_Describe.RData')

powerBSD = matrix(powerBSD/100, 5, 4)
powerWBB = matrix(powerWBB/100, 5, 4)

TRBSD = matrix(apply(TrateBSD, 2, mean), 5, 4)
TRWBB = matrix(apply(TrateWBB, 2, mean), 5, 4)

FRBSD = matrix(apply(FrateBSD, 2, mean), 5, 4)
FRWBB = matrix(apply(FrateWBB, 2, mean), 5, 4)

DTBSD = matrix(apply(DifTiBSD, 2, mean), 5, 4)
DTWBB = matrix(apply(DifTiWBB, 2, mean), 5, 4)


delta = c(0.20, 0.25, 0.30, 0.35, 0.40)
for (j in 1:4)
{
  nameTj = paste0('Figures/WNS-TPR_beta', j, '.jpg')
  nameFj = paste0('Figures/WNS-FDR_beta', j, '.jpg')
  
  mainTj = paste0('WNS, beta = ', j, ' TPR')
  mainFj = paste0('WNS, beta = ', j, ' FDR')
  
  jpeg(filename = nameTj, width = 1200, height = 900, quality = 100)
  plot(delta, TRBSD[ ,j], type = 'b', pch = 1, lty = 1, lwd = 4, col = 1, main = mainTj, 
       ylab = "", xlab = "delta", cex.lab = 3, cex.axis = 2.8, ylim = c(0, 1), cex.main = 4)
  lines(delta, TRWBB[ ,j], type = "b", pch = 2, lty = 1, col = 2, lwd = 4)
  legend("right", c("BSD", "WBB"), pch = c(1, 2), lty = 1, box.col = "grey", cex = 3.5, col = c(1, 2), lwd = 4)
  dev.off()
  
  jpeg(filename = nameFj, width = 1200, height = 900, quality = 100)
  plot(delta, FRBSD[ ,j], type = 'b', pch = 1, lty = 1, lwd = 4, col = 1, main = mainFj, 
       ylab = "", xlab = "delta", cex.lab = 3, cex.axis = 2.8, ylim = c(0, 0.5), cex.main = 4)
  lines(delta, FRWBB[ ,j], type = "b", pch = 2, lty = 1, col = 2, lwd = 4)
  legend("topright", c("BSD", "WBB"), pch = c(1, 2), lty = 1, box.col = "grey", cex = 3.5, col = c(1, 2), lwd = 4)
  dev.off()
}

#############################################################################################################################################################################################
