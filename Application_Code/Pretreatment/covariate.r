# R program ukb51090.tab created 2022-05-10 by ukb2r.cpp Mar 14 2018 14:22:05

bd <- read.table("/home/zhw98/Biobank/covariate.tab", header=TRUE, sep="\t")
lvl.0009 <- c(0,1)
lbl.0009 <- c("Female","Male")
bd$f.31.0.0 <- ordered(bd$f.31.0.0, levels=lvl.0009, labels=lbl.0009)

attach(bd)
covariate = data.frame(f.54.0.0, f.21022.0.0)
save(covariate, file = 'covariates.RData')
detach(bd)
