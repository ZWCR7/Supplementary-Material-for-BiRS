rm(list = ls())
deltav = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7)

TrateBSD = matrix(0, 500, 7)
FrateBSD = matrix(0, 500, 7)

TrateBSC = matrix(0, 500, 7)
FrateBSC = matrix(0, 500, 7)

TrateQSCAN = matrix(0, 500, 7)
FrateQSCAN = matrix(0, 500, 7)

TrateSTAARO = matrix(0, 500, 7)
FrateSTAARO = matrix(0, 500, 7)

TrateSTAARS = matrix(0, 500, 7)
FrateSTAARS = matrix(0, 500, 7)

TrateSTAARB = matrix(0, 500, 7)
FrateSTAARB = matrix(0, 500, 7)

TrateKnockoffScreen = matrix(0, 500, 7)
FrateKnockoffScreen = matrix(0, 500, 7)

Trate4S = matrix(0, 500, 7)
Frate4S = matrix(0, 500, 7)


for (i in 1:7)
{ 
  load(paste0('Genetic_RES/BiRS-DCF_delta = ', deltav[i], '.RData'))
  load(paste0('Genetic_RES/BiRS-CL_delta = ', deltav[i], '.RData'))
  load(paste0('Genetic_RES/QSCAN_delta = ', deltav[i], '.RData'))
  load(paste0('Genetic_RES/STAARO_delta = ', deltav[i], '.RData'))
  load(paste0('Genetic_RES/STAARS_delta = ', deltav[i], '.RData'))
  load(paste0('Genetic_RES/STAARB_delta = ', deltav[i], '.RData'))
  load(paste0('Genetic_RES/KnockoffScreen_delta = ', deltav[i], '.RData'))
  load(paste0('Genetic_RES/4S_delta = ', deltav[i], '.RData'))
  
  nsimu = 500; p = length(true_causal)
  EstBSD = matrix(0, nsimu, p)
  EstBSC = matrix(0, nsimu, p)
  EstQSCAN = matrix(0, nsimu, p)
  EstSTAARO = matrix(0, nsimu, p)
  EstSTAARS = matrix(0, nsimu, p)
  EstSTAARB = matrix(0, nsimu, p)
  EstKnockoffScreen = matrix(0, nsimu, p)
  Est4S = matrix(0, nsimu, p)
  
  for (j in 1:500)
  {
    if (length(BiRSDCF_Det[[j]]) != 0)
    {
      Resj = BiRSDCF_Det[[j]]
      for (k in 1:nrow(Resj))
      {
        EstBSD[j, (Resj$startind[k]:Resj$endind[k])] = 1
      }
    }
    
    #############################################################################
    if (length(BiRSCL_Det[[j]]) != 0)
    {
      Resj = BiRSCL_Det[[j]]
      for (k in 1:nrow(Resj))
      {
        EstBSC[j, (Resj$startind[k]:Resj$endind[k])] = 1
      }
    }

    #############################################################################
    if (length(QSCAN_Det[[j]]) != 0)
    {
      Resj = QSCAN_Det[[j]]
      for (k in 1:nrow(Resj))
      {
        EstQSCAN[j, (Resj$startind[k]:Resj$endind[k])] = 1
      }
    }
    
    #############################################################################
    if (length(STAARO_Det[[j]]) != 0)
    {
      Resj = STAARO_Det[[j]]
      for (k in 1:nrow(Resj))
      {
        EstSTAARO[j, (Resj$startind[k]:Resj$endind[k])] = 1
      }
    }
    
    #############################################################################
    if (length(STAARS_Det[[j]]) != 0)
    {
      Resj = STAARS_Det[[j]]
      for (k in 1:nrow(Resj))
      {
        EstSTAARS[j, (Resj$startind[k]:Resj$endind[k])] = 1
      }
    }
    
    #############################################################################
    if (length(STAARB_Det[[j]]) != 0)
    {
      Resj = STAARB_Det[[j]]
      for (k in 1:nrow(Resj))
      {
        EstSTAARB[j, (Resj$startind[k]:Resj$endind[k])] = 1
      }
    }
    
    #############################################################################
    if (length(KnockoffScreen_Det[[j]]) != 0)
    {
      Resj = KnockoffScreen_Det[[j]]
      for (k in 1:nrow(Resj))
      {
        EstKnockoffScreen[j, (Resj$startind[k]:Resj$endind[k])] = 1
      }
    }
    
    #############################################################################
    if (length(SSSS_Det[[j]]) != 0)
    {
      Resj = SSSS_Det[[j]]
      for (k in 1:nrow(Resj))
      {
        Est4S[j, (Resj$startind[k]:Resj$endind[k])] = 1
      }
    }
    
  }
  
  true.ind = rep(0, p)
  true.ind[which(true_causal != 0)] = 1
  true.ind[which(true_causal == 0)] = -1
  
  for (l in 1:25)
  {
    if (length(which(EstBSD[l, ] != 0)) != 0)
    {
      Diffl = EstBSD[l, ] - true.ind
      
      TrateBSD[l, i] = length(which(Diffl == 0))/length(which(true_causal != 0))
      FrateBSD[l, i] = length(which(Diffl == 3))/length(which(EstBSD[l, ] != 0))
    }
    
    ############################################################################
    if (length(which(EstBSC[l, ] != 0)) != 0)
    {
      Diffl = EstBSC[l, ] - true.ind
       
      TrateBSC[l, i] = length(which(Diffl == 0))/length(which(true_causal != 0))
      FrateBSC[l, i] = length(which(Diffl == 3))/length(which(EstBSC[l, ] != 0))
    }

    ############################################################################
    if (length(which(EstQSCAN[l, ] != 0)) != 0)
    {
      Diffl = EstQSCAN[l, ] - true.ind
       
      TrateQSCAN[l, i] = length(which(Diffl == 0))/length(which(true_causal != 0))
      FrateQSCAN[l, i] = length(which(Diffl == 3))/length(which(EstQSCAN[l, ] != 0))
    }
    
    ############################################################################
    if (length(which(EstSTAARO[l, ] != 0)) != 0)
    {
      Diffl = EstSTAARO[l, ] - true.ind
       
      TrateSTAARO[l, i] = length(which(Diffl == 0))/length(which(true_causal != 0))
      FrateSTAARO[l, i] = length(which(Diffl == 3))/length(which(EstSTAARO[l, ] != 0))
    }
    
    ############################################################################
    if (length(which(EstSTAARS[l, ] != 0)) != 0)
    {
      Diffl = EstSTAARS[l, ] - true.ind
       
      TrateSTAARS[l, i] = length(which(Diffl == 0))/length(which(true_causal != 0))
      FrateSTAARS[l, i] = length(which(Diffl == 3))/length(which(EstSTAARS[l, ] != 0))
    }
    
    ############################################################################
    if (length(which(EstSTAARB[l, ] != 0)) != 0)
    {
      Diffl = EstSTAARB[l, ] - true.ind
       
      TrateSTAARB[l, i] = length(which(Diffl == 0))/length(which(true_causal != 0))
      FrateSTAARB[l, i] = length(which(Diffl == 3))/length(which(EstSTAARB[l, ] != 0))
    }
    ############################################################################
    
    if (length(which(EstKnockoffScreen[l, ] != 0)) != 0)
    {
      Diffl = EstKnockoffScreen[l, ] - true.ind
       
      TrateKnockoffScreen[l, i] = length(which(Diffl == 0))/length(which(true_causal != 0))
      FrateKnockoffScreen[l, i] = length(which(Diffl == 3))/length(which(EstKnockoffScreen[l, ] != 0))
    }
    
    ############################################################################
    
    if (length(which(Est4S[l, ] != 0)) != 0)
    {
      Diffl = Est4S[l, ] - true.ind
       
      Trate4S[l, i] = length(which(Diffl == 0))/length(which(true_causal != 0))
      Frate4S[l, i] = length(which(Diffl == 3))/length(which(Est4S[l, ] != 0))
    }
    
  }
}

TRBSD = apply(TrateBSD, 2, mean)
TRBSC = apply(TrateBSC, 2, mean)
TRKSD = apply(TrateKnockoffScreen, 2, mean)
TRSCQ = apply(TrateQSCAN, 2, mean)
TRSTAARO = apply(TrateSTAARO, 2, mean)
TRSTAARB = apply(TrateSTAARB, 2, mean)
TRSTAARS = apply(TrateSTAARS, 2, mean)
TR4S = apply(Trate4S, 2, mean)

FRBSD = apply(FrateBSD, 2, mean)
FRBSC = apply(FrateBSC, 2, mean)
FRKSD = apply(FrateKnockoffScreen, 2, mean)
FRSCQ = apply(FrateQSCAN, 2, mean)
FRSTAARO = apply(FrateSTAARO, 2, mean)
FRSTAARB = apply(FrateSTAARB, 2, mean)
FRSTAARS = apply(FrateSTAARS, 2, mean)
FR4S = apply(Frate4S, 2, mean)

delta = c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70)
nameTj = paste0('Figure-Genetic/Genetic-TPR.jpg')
nameFj = paste0('Figure-Genetic/Genetic-FDR.jpg')

# nameCTj = paste0('Figures/MES-CP-TPR_beta', j, '.jpg')
# nameCFj = paste0('Figures/MES-CP-FDR_beta', j, '.jpg')

mainTj = paste0('Genetic Setting, TPR')
mainFj = paste0('Genetic Setting, FDR')

jpeg(filename = nameTj, width = 1200, height = 900, quality = 100)
plot(delta, TRBSD, type = 'b', pch = 1, lty = 1, lwd = 4, col = 1, main = mainTj, 
     ylab = "", xlab = "delta", cex.lab = 3, cex.axis = 2.8, ylim = c(0, 1), cex.main = 4)
lines(delta, TRBSC, type = "b", pch = 2, lty = 1, col = 2, lwd = 4)
lines(delta, TRSCQ, type = "b", pch = 3, lty = 1, col = 3, lwd = 4)
lines(delta, TRKSD, type = "b", pch = 4, lty = 1, col = 4, lwd = 4)
# lines(delta, TRSTAARO, type = "b", pch = 5, lty = 1, col = 5, lwd = 4)
# lines(delta, TRSTAARB, type = "b", pch = 6, lty = 1, col = 6, lwd = 4)
lines(delta, TR4S, type = "b", pch = 5, lty = 1, col = 5, lwd = 4)
lines(delta, TRSTAARS, type = "b", pch = 6, lty = 1, col = 6, lwd = 4)
#legend("topleft", c("BSD","BSC", "SCQ", "KSD", "STA-O", "STA-B", "STA-S", "4S"), pch = c(1, 2, 3, 4, 5, 6, 7, 8), lty = 1, box.col = "grey", cex = 3.5, col = c(1, 2, 3, 4, 5, 6, 7, 8), lwd = 4)
dev.off()

jpeg(filename = nameFj, width = 1200, height = 900, quality = 100)
plot(delta, FRBSD, type = 'b', pch = 1, lty = 1, lwd = 4, col = 1, main = mainFj, 
     ylab = "", xlab = "delta", cex.lab = 3, cex.axis = 2.8, ylim = c(0, 1), cex.main = 4)
lines(delta, FRBSC, type = "b", pch = 2, lty = 1, col = 2, lwd = 4)
lines(delta, FRSCQ, type = "b", pch = 3, lty = 1, col = 3, lwd = 4)
lines(delta, FRKSD, type = "b", pch = 4, lty = 1, col = 4, lwd = 4)
# lines(delta, FRSTAARO, type = "b", pch = 5, lty = 1, col = 5, lwd = 4)
# lines(delta, FRSTAARB, type = "b", pch = 6, lty = 1, col = 6, lwd = 4)
lines(delta, FR4S, type = "b", pch = 5, lty = 1, col = 5, lwd = 4)
lines(delta, FRSTAARS, type = "b", pch = 6, lty = 1, col = 6, lwd = 4)
legend("topleft", c("BiRS-DCF","BiRS-CL", "Q-SCAN", "KSD", "4S", "SCAN-STA-S"), pch = c(1, 2, 3, 4, 5, 6), lty = 1, box.col = "grey", cex = 3.5, col = c(1, 2, 3, 4, 5, 6), lwd = 4)
dev.off()

