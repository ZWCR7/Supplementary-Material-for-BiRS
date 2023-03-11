rm(list = ls())

load('Normal_ESBA_Decay/BSD_D_Describe.RData')
load('Normal_ESBA_Decay/BSD_ND_Describe.RData')
 
TRBSD_ND = apply(TrateBSD_ND, 2, mean)
TRBSD_D = apply(TrateBSD_D, 2, mean)

FRBSD_ND = apply(FrateBSD_ND, 2, mean)
FRBSD_D = apply(FrateBSD_D, 2, mean)

delta = c(0.20, 0.25, 0.30, 0.35, 0.40)
 
nameT = 'Figures/Normal_ESBA_TPR_Decay.jpg'
mainT = 'TPR, Unequal Covariance'

nameF = 'Figures/Normal_ESBA_FDR_Decay.jpg'
mainF = 'FDR, Unequal Covariance'

jpeg(filename = nameT, width = 1200, height = 900)
plot(delta, TRBSD_ND, type = 'b', pch = 1, lty = 1, lwd = 2, col = 1, main = mainT, 
     ylab = "", xlab = "delta", cex.lab = 3, cex.axis = 2.8, ylim = c(0, 1), cex.main = 4)
lines(delta, TRBSD_D, type = "b", pch = 2, lty = 2, col = 1, lwd = 4)
legend("topleft", c("Non-Decay", "Decay"), pch = c(1, 2), lty = c(1, 2), box.col = "grey", cex = 4, col = 1, lwd = 4)
dev.off()

jpeg(filename = nameF, width = 1200, height = 900)
plot(delta, FRBSD_ND, type = 'b', pch = 1, lty = 1, lwd = 4, col = 1, main = mainF, 
     ylab = "", xlab = "delta", cex.lab = 3, cex.axis = 2.8, ylim = c(0, 0.3), cex.main = 4)
lines(delta, FRBSD_D, type = "b", pch = 2, lty = 2, col = 1, lwd = 2)
legend("topleft", c("Non-Decay", "Decay"), pch = c(1, 2), lty = c(1, 2), box.col = "grey", cex = 4, col = 1, lwd = 4)
dev.off()
############################################################################################################################

rm(list = ls())

load('Normal_NSBA_Decay/BSD_D_Describe.RData')
load('Normal_NSBA_Decay/BSD_ND_Describe.RData')

TRBSD_ND = apply(TrateBSD_ND, 2, mean)
TRBSD_D = apply(TrateBSD_D, 2, mean)

FRBSD_ND = apply(FrateBSD_ND, 2, mean)
FRBSD_D = apply(FrateBSD_D, 2, mean)

delta = c(0.20, 0.25, 0.30, 0.35, 0.40)

nameT = 'Figures/Normal_NSBA_TPR_Decay.jpg'
mainT = 'TPR, Unequal Covariance'

nameF = 'Figures/Normal_NSBA_FDR_Decay.jpg'
mainF = 'FDR, Unequal Covariance'

jpeg(filename = nameT, width = 1200, height = 900)
plot(delta, TRBSD_ND, type = 'b', pch = 1, lty = 1, lwd = 2, col = 1, main = mainT, 
     ylab = "", xlab = "delta", cex.lab = 3, cex.axis = 2.8, ylim = c(0, 1), cex.main = 4)
lines(delta, TRBSD_D, type = "b", pch = 2, lty = 2, col = 1, lwd = 4)
legend("topleft", c("Non-Decay", "Decay"), pch = c(1, 2), lty = c(1, 2), box.col = "grey", cex = 4, col = 1, lwd = 4)
dev.off()

jpeg(filename = nameF, width = 1200, height = 900)
plot(delta, FRBSD_ND, type = 'b', pch = 1, lty = 1, lwd = 4, col = 1, main = mainF, 
     ylab = "", xlab = "delta", cex.lab = 3, cex.axis = 2.8, ylim = c(0, 0.3), cex.main = 4)
lines(delta, FRBSD_D, type = "b", pch = 2, lty = 2, col = 1, lwd = 2)
legend("topleft", c("Non-Decay", "Decay"), pch = c(1, 2), lty = c(1, 2), box.col = "grey", cex = 4, col = 1, lwd = 4)
dev.off()