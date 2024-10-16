rm(list = ls())

truncv = 3：6
snum = length(truncv)

Trate = matrix(0, 100, snum)
Frate = matrix(0, 100, snum)

load('Review-Sensitive/MES.RData')
nsimu = length(res); p = 8192

true.ind = rep(0, p)
true.ind[which(mu != 0)] = 1
true.ind[which(mu == 0)] = -2

for (i in 1:snum)
{
  EstR = matrix(0, nsimu, p)
  for (j in 1:nsimu)
  {
    if (!is.character(res[[j]][[i]]))
    {
      Resj = res[[j]][[i]]$BSDCF_res
      
      for (k in 1:nrow(Resj))
      {
        EstR[j, (Resj$startind[k]:Resj$endind[k])] = 1
      }
      
    }
  }
  
  for (l in 1:nsimu)
  {
    if (length(which(EstR[l, ] != 0)) != 0)
    {
      Diffl = EstR[l, ] - true.ind
      
      Trate[l, i] = length(which(Diffl == 0))/length(which(mu != 0))
      Frate[l, i] = length(which(Diffl == 3))/length(which(EstR[l, ] != 0))
    }
  }
  
}

delta = 3：6

TPR = apply(Trate, 2, mean)
FDR = apply(Frate, 2, mean)

nameTR = paste0('Figures/MES-Sensitive.jpg')
mainTR = paste0('MES, beta = 4, Sensitive')

jpeg(filename = nameTR, width = 1200, height = 900, quality = 100)
plot(delta, TPR, type = 'b', pch = 1, lty = 1, lwd = 4, col = 1, main = mainTR, 
     ylab = "", xlab = "s", cex.lab = 3, cex.axis = 2.8, ylim = c(0, 1), cex.main = 4)
lines(delta, FDR, type = "b", pch = 2, lty = 1, col = 2, lwd = 4)
legend("right", c("TPR", "FDR"), pch = c(1, 2, 3), lty = 1, box.col = "grey", cex = 3.5, col = c(1, 2, 3), lwd = 4)
dev.off()
##################################################################################################################################

rm(list = ls())

truncv = 3:6
snum = length(truncv)

Trate = matrix(0, 100, snum)
Frate = matrix(0, 100, snum)

load('Review-Sensitive/MNS.RData')
nsimu = length(res); p = 8192

true.ind = rep(0, p)
true.ind[which(mu != 0)] = 1
true.ind[which(mu == 0)] = -2

for (i in 1:snum)
{
  EstR = matrix(0, nsimu, p)
  for (j in 1:nsimu)
  {
    if (!is.character(res[[j]][[i]]))
    {
      Resj = res[[j]][[i]]$BSDCF_res
      
      for (k in 1:nrow(Resj))
      {
        EstR[j, (Resj$startind[k]:Resj$endind[k])] = 1
      }
      
    }
  }
  
  for (l in 1:nsimu)
  {
    if (length(which(EstR[l, ] != 0)) != 0)
    {
      Diffl = EstR[l, ] - true.ind
      
      Trate[l, i] = length(which(Diffl == 0))/length(which(mu != 0))
      Frate[l, i] = length(which(Diffl == 3))/length(which(EstR[l, ] != 0))
    }
  }
  
}

delta = 3:6

TPR = apply(Trate, 2, mean)
FDR = apply(Frate, 2, mean)

nameTR = paste0('Figures/MNS-Sensitive.jpg')
mainTR = paste0('MNS, beta = 4, Sensitive')

jpeg(filename = nameTR, width = 1200, height = 900, quality = 100)
plot(delta, TPR, type = 'b', pch = 1, lty = 1, lwd = 4, col = 1, main = mainTR, 
     ylab = "", xlab = "s", cex.lab = 3, cex.axis = 2.8, ylim = c(0, 1), cex.main = 4)
lines(delta, FDR, type = "b", pch = 2, lty = 1, col = 2, lwd = 4)
legend("right", c("TPR", "FDR"), pch = c(1, 2, 3), lty = 1, box.col = "grey", cex = 3.5, col = c(1, 2, 3), lwd = 4)
dev.off()
############################################################################################################################

rm(list = ls())

truncv = 3:6
snum = length(truncv)

Trate = matrix(0, 100, snum)
Frate = matrix(0, 100, snum)

load('Review-Sensitive/WES.RData')
nsimu = length(res); p = 8192

true.ind = rep(0, p)
true.ind[which(mu != 0)] = 1
true.ind[which(mu == 0)] = -2

for (i in 1:snum)
{
  EstR = matrix(0, nsimu, p)
  for (j in 1:nsimu)
  {
    if (!is.character(res[[j]][[i]]))
    {
      Resj = res[[j]][[i]]$BSDCF_res
      
      for (k in 1:nrow(Resj))
      {
        EstR[j, (Resj$startind[k]:Resj$endind[k])] = 1
      }
      
    }
  }
  
  for (l in 1:nsimu)
  {
    if (length(which(EstR[l, ] != 0)) != 0)
    {
      Diffl = EstR[l, ] - true.ind
      
      Trate[l, i] = length(which(Diffl == 0))/length(which(mu != 0))
      Frate[l, i] = length(which(Diffl == 3))/length(which(EstR[l, ] != 0))
    }
  }
  
}

delta = 3:6

TPR = apply(Trate, 2, mean)
FDR = apply(Frate, 2, mean)

nameTR = paste0('Figures/WES-Sensitive.jpg')
mainTR = paste0('WES, beta = 4, Sensitive')

jpeg(filename = nameTR, width = 1200, height = 900, quality = 100)
plot(delta, TPR, type = 'b', pch = 1, lty = 1, lwd = 4, col = 1, main = mainTR, 
     ylab = "", xlab = "s", cex.lab = 3, cex.axis = 2.8, ylim = c(0, 1), cex.main = 4)
lines(delta, FDR, type = "b", pch = 2, lty = 1, col = 2, lwd = 4)
legend("right", c("TPR", "FDR"), pch = c(1, 2, 3), lty = 1, box.col = "grey", cex = 3.5, col = c(1, 2, 3), lwd = 4)
dev.off()

################################################################################################################################

rm(list = ls())

truncv = 3:6
snum = length(truncv)

Trate = matrix(0, 100, snum)
Frate = matrix(0, 100, snum)

load('Review-Sensitive/WNS.RData')
nsimu = length(res); p = 8192

true.ind = rep(0, p)
true.ind[which(mu != 0)] = 1
true.ind[which(mu == 0)] = -2

for (i in 1:snum)
{
  EstR = matrix(0, nsimu, p)
  for (j in 1:nsimu)
  {
    if (!is.character(res[[j]][[i]]))
    {
      Resj = res[[j]][[i]]$BSDCF_res
      
      for (k in 1:nrow(Resj))
      {
        EstR[j, (Resj$startind[k]:Resj$endind[k])] = 1
      }
      
    }
  }
  
  for (l in 1:nsimu)
  {
    if (length(which(EstR[l, ] != 0)) != 0)
    {
      Diffl = EstR[l, ] - true.ind
      
      Trate[l, i] = length(which(Diffl == 0))/length(which(mu != 0))
      Frate[l, i] = length(which(Diffl == 3))/length(which(EstR[l, ] != 0))
    }
  }
  
}

delta = 3:6

TPR = apply(Trate, 2, mean)
FDR = apply(Frate, 2, mean)

nameTR = paste0('Figures/WNS-Sensitive.jpg')
mainTR = paste0('WNS, beta = 4, Sensitive')

jpeg(filename = nameTR, width = 1200, height = 900, quality = 100)
plot(delta, TPR, type = 'b', pch = 1, lty = 1, lwd = 4, col = 1, main = mainTR, 
     ylab = "", xlab = "s", cex.lab = 3, cex.axis = 2.8, ylim = c(0, 1), cex.main = 4)
lines(delta, FDR, type = "b", pch = 2, lty = 1, col = 2, lwd = 4)
legend("right", c("TPR", "FDR"), pch = c(1, 2, 3), lty = 1, box.col = "grey", cex = 3.5, col = c(1, 2, 3), lwd = 4)
dev.off()
