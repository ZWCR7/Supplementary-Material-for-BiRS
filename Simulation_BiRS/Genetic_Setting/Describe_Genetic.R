powerBSD = rep(500, 20)
TrateBSD = matrix(0, 500, 20)
FrateBSD = matrix(0, 500, 20)
DifTiBSD = matrix(0, 500, 20)

for (i in 1:20)
{ 
  load(paste0('Mulist', i, '.RData'))
  
  nsimu = length(res); p = 8192
  EstR = matrix(0, nsimu, p)
  
  for (j in 1:length(res))
  {
    if (!is.character(res[[j]]$reBSD))
    {
      Resj = res[[j]]$reBSD$BSDCF_res
      
      DifTiBSD[j, i] = as.numeric(res[[j]]$diffBSD, units = 'secs')
      
      for (k in 1:nrow(Resj))
      {
        EstR[j, (Resj$startind[k]:Resj$endind[k])] = 1
      }
    }
    else
    {
      powerBSD[i] = powerBSD[i] - 1
    }
  }
  
  true.ind = rep(0, p)
  true.ind[which(mu != 0)] = 1
  true.ind[which(mu == 0)] = -2
  
  for (l in 1:nsimu)
  {
    if (length(which(EstR[l, ] != 0)) != 0)
    {
      Diffl = EstR[l, ] - true.ind
      
      TrateBSD[l, i] = length(which(Diffl == 0))/length(which(mu != 0))
      FrateBSD[l, i] = length(which(Diffl == 3))/length(which(EstR[l, ] != 0))
    }
  }
}

save(list = c('FrateBSD', 'TrateBSD', 'powerBSD', 'DifTiBSD'), file = 'BSD_Describe.RData')

#####################################################################################################

rm(list = ls())

powerBSC = rep(500, 20)
TrateBSC = matrix(0, 500, 20)
FrateBSC = matrix(0, 500, 20)
DifTiBSC = matrix(0, 500, 20)

for (i in 1:20)
{
  load(paste0('Mulist', i, '.RData'))
  
  nsimu = length(res); p = 8192
  EstR = matrix(0, nsimu, p)
  
  for (j in 1:length(res))
  {
    
    if (!is.character(res[[j]]$reBSC))
    {
      Resj = res[[j]]$reBSC$BSCL_res
      
      DifTiBSC[j, i] = as.numeric(res[[j]]$diffBSC, units = 'secs')
      
      for (k in 1:nrow(Resj))
      {
        EstR[j, (Resj$startind[k]:Resj$endind[k])] = 1
      }
    }
    else
    {
      powerBSC[i] = powerBSC[i] - 1
    }
  }
  
  true.ind = rep(0, p)
  true.ind[which(mu != 0)] = 1
  true.ind[which(mu == 0)] = -2
  
  for (l in 1:nsimu)
  {
    if (length(which(EstR[l, ] != 0)) != 0)
    {
      Diffl = EstR[l, ] - true.ind
      
      TrateBSC[l, i] = length(which(Diffl == 0))/length(which(mu != 0))
      FrateBSC[l, i] = length(which(Diffl == 3))/length(which(EstR[l, ] != 0))
    }
  }
  
}

save(list = c('FrateBSC', 'TrateBSC', 'powerBSC', 'DifTiBSC'), file = 'BSC_Describe.RData')
#######################################################################################################

rm(list = ls())

powerSCQ = rep(500, 20)
TrateSCQ = matrix(0, 500, 20)
FrateSCQ = matrix(0, 500, 20)
DifTiSCQ = matrix(0, 500, 20)

for (i in 1:20)
{
  load(paste0('Mulist', i, '.RData'))
  
  nsimu = length(res); p = 8192
  EstR = matrix(0, nsimu, p)
  
  for (j in 1:length(res))
  {
    Resj = as.matrix(res[[j]]$reSCQ$SCAN_res)
    
    if (ncol(Resj) == 1) Resj = t(Resj)
    
    if (sum(Resj) == 1)
    {
      powerSCQ[i] = powerSCQ[i] - 1
    }
    else
    {
      DifTiSCQ[j, i] = as.numeric(res[[j]]$diffSCQ, units = 'secs')
      
      for (k in 1:nrow(Resj))
      {
        EstR[j, (Resj[k, 2]:Resj[k, 3])] = 1
      }
    }
  }
  
  true.ind = rep(0, p)
  true.ind[which(mu != 0)] = 1
  true.ind[which(mu == 0)] = -2
  
  for (l in 1:nsimu)
  {
    if (length(which(EstR[l, ] != 0)) != 0)
    {
      Diffl = EstR[l, ] - true.ind
      
      TrateSCQ[l, i] = length(which(Diffl == 0))/length(which(mu != 0))
      FrateSCQ[l, i] = length(which(Diffl == 3))/length(which(EstR[l, ] != 0))
    }
  }
  
}

save(list = c('FrateSCQ', 'TrateSCQ', 'powerSCQ', 'DifTiSCQ'), file = 'SCQ_Describe.RData')

######################################################################################################################################

rm(list = ls())

MK.q.byStat<-function (kappa,tau,M,Rej.Bound=10000){
  b<-order(tau,decreasing=T)
  c_0<-kappa[b]==0
  #calculate ratios for top Rej.Bound tau values
  ratio<-c();temp_0<-0
  for(i in 1:length(b)){
    #if(i==1){temp_0=c_0[i]}
    temp_0<-temp_0+c_0[i]
    temp_1<-i-temp_0
    temp_ratio<-(1/M+1/M*temp_1)/max(1,temp_0)
    ratio<-c(ratio,temp_ratio)
    if(i>Rej.Bound){break}
  }
  #calculate q values for top Rej.Bound values
  q<-rep(1,length(tau));index_bound<-max(which(tau[b]>0))
  for(i in 1:length(b)){
    temp.index<-i:min(length(b),Rej.Bound,index_bound)
    if(length(temp.index)==0){next}
    q[b[i]]<-min(ratio[temp.index])*c_0[i]+1-c_0[i]
    if(i>Rej.Bound){break}
  }
  return(q)
}


powerKSD = rep(500, 20)
TrateKSD = matrix(0, 500, 20)
FrateKSD = matrix(0, 500, 20)
DifTiKSD = matrix(0, 500, 20)

for (i in 1:20)
{
  load(paste0('Mulist', i, '.RData'))
  
  nsimu = length(res); p = 8192
  EstR = matrix(0, nsimu, p)
  
  for (j in 1:nsimu)
  {
    result.single = res[[j]]$reKSD$result.single
    result.window = res[[j]]$reKSD$result.window
    
    if (is.null(result.single))
    {
      result = result.window
    }
    else
    {
      temp<-result.window[,match(colnames(result.single),colnames(result.window))]
      
      result<-rbind(result.single,temp)
    }
    
    result<-result[order(result[,2]),]
    result<-result[order(result[,1]),]
    
    q<-MK.q.byStat(result[,'kappa'],result[,'tau'],M=5)
    result.summary<-cbind(result[,1:5],q,result[,-(1:5)])
    colnames(result.summary)[6]<-'Qvalue'
    
    result.summary = as.data.frame(result.summary)
    indsel = which(result.summary$Qvalue < 0.05)
    
    if (length(indsel != 0))
    {
      for (k in 1:length(indsel))
      {
        locate = result.summary$actual_start[indsel[k]]:result.summary$actual_end[indsel[k]]
        EstR[j, locate] = 1
      }
      
      DifTiKSD[j, i] = as.numeric(res[[j]]$diffKSD, units = 'secs')
    }
    else
    {
      powerKSD[i] = powerKSD[i] - 1
    }
    
  }
  
  true.ind = rep(0, 8192)
  true.ind[which(mu != 0)] = 1
  true.ind[which(mu == 0)] = -2
  
  for (l in 1:nsimu)
  {
    if (length(which(EstR[l, ] != 0)) != 0)
    {
      Diffl = EstR[l, ] - true.ind
      
      TrateKSD[l, i] = length(which(Diffl == 0))/length(which(mu != 0))
      FrateKSD[l, i] = length(which(Diffl == 3))/length(which(EstR[l, ] != 0))
    }
  }
}

save(list = c('FrateKSD', 'TrateKSD', 'powerKSD', 'DifTiKSD'), file = 'KSD_Describe.RData')


