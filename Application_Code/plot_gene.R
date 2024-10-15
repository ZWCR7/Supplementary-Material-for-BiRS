library(ggplot2)
library(ggbreak)
library(Unicode)
###########################################################################################################################

rm(list = ls())

BiRS = read.csv('C50Breast_BiRS.csv')
SCAN = read.csv('C50Breast_Scan.csv')
KSCN = read.csv('C50Breast_KSGWAS.csv')
SSSS = read.csv('C50Breast_4S.csv')
REGN = read.csv('C50Breast_Regenie.csv')
chr11_gene = read.csv('C50Breast_chr11_gene.csv')

BiRS = BiRS[, -1]; SCAN = SCAN[, -1]; KSCN = KSCN[, -1]; SSSS = SSSS[, -1]; REGN = REGN[, -1]

BiRSChr11 = BiRS[which(BiRS$chrsel == 11), ]
SCANChr11 = SCAN[which(SCAN$chrsel == 11), ]
KSCNChr11 = KSCN[which(KSCN$chr == 11), ]
SSSSChr11 = SSSS[which(SSSS$chr == 11), ]
REGNChr11 = REGN[which(REGN$CHR == 11), ]

chr11_len = nrow(chr11_gene)
startposchr11 = chr11_gene$Start
endposchr11 = chr11_gene$End
deschr11 = chr11_gene$Gene

# jpeg(filename = 'Chr11_Part3.png', width = 2000, height = 1000, quality = 100)
# plot(x = BiRSChr11$possel, y = rep(2, length(BiRSChr11$possel)), type = 'p', cex = 1.5, col = 'blue', pch = 16, ylim = c(0, 4), cex.lab = 1.5, main = 'Chromosomes 11', cex.main = 2, 
#      xlim = c(startposchr11[1711], endposchr11[1720]), xlab = 'Genomic Position', ylab = '')

# plot(x = startposchr1[1]:endposchr1[1], y = rep(1, endposchr1[1] - startposchr1[1] + 1), ylim = c(0, 4), type = 'l', lwd = 8, xlim = c(startposchr1[1], endposchr1[chr1_len]), xlab = 'Genomic Position', ylab = '', 
#      cex.lab = 1.5, main = 'Chromosomes 1', cex.main = 2)

x1start = 1850000; x1end = 2000000; 
x2start = 69200000; x2end = 69500000; 
x3start = 129200000; x3end = 129600000;

gene_sel = c('LSP1', 'TNNT3', 'LINC02953', 'LINC01488', 'CCND1', 'LTO1', 'BARX2', 'LINC01395')
xgene = NULL
textpos = rep(0, length(gene_sel))
for (i in 1:length(gene_sel))
{
  ind_seli = which(deschr11 == gene_sel[i])
  xgene = c(xgene, startposchr11[ind_seli]:endposchr11[ind_seli])
  textpos[i] = (startposchr11[ind_seli] + endposchr11[ind_seli])/2
}

notation = c(rep('BiRS-DCF', length(BiRSChr11$possel)), rep("Q-SCAN", length(SCANChr11$possel)), rep('KnockoffScreen', length(KSCNChr11$bp)), rep('4S', length(SSSSChr11$possel)) , rep("REGENIE", length(REGNChr11$BP)))
notation = factor(notation, levels = c('BiRS-DCF', 'Q-SCAN', 'KnockoffScreen', '4S', 'REGENIE'))
X = c(BiRSChr11$possel, SCANChr11$possel, KSCNChr11$bp, SSSSChr11$possel, REGNChr11$BP)
Y = c(rep(6, length(BiRSChr11$possel)), rep(5, length(SCANChr11$possel)), 
      rep(4, length(KSCNChr11$bp)), rep(3, length(SSSSChr11$possel)), rep(2, length(REGNChr11$BP))) 

CHR11Data = data.frame(X, Y, notation)
Pos = rep(1, length(xgene))
GeneData = data.frame(xgene, Pos)

region = ggplot(data = CHR11Data, aes(x = X, y = Y)) +
         geom_point(stat = 'identity', aes(color = notation), size = 4) + 
         geom_point(data = GeneData, mapping = aes(x = xgene, y = Pos), size = 4) + 
         ggtitle("B.Location of Significant Regions in Chromosome 11") +
         xlab("Genomic Location(kb)") + 
         ylab(" ") + 
         labs(color = " ")

for(i in 1:length(gene_sel))
{
  region = region + annotate("text", x = textpos[i], y = 0.6, label = gene_sel[i], 
                             size = 5, fontface = 'bold.italic')
}

region = region + theme_gray() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(axis.text.x.top = element_blank(), axis.ticks.x.top = element_blank()) + 
  theme(plot.title = element_text(face = "bold", size = 24), axis.title.x = element_text(face = "bold", size = 20), 
        legend.text = element_text(face = "bold", size = 17), axis.text.x = element_text(face = "bold", size = 17), 
        legend.position = "none") +
  scale_x_break(c(x1end, x2start)) +
  scale_x_break(c(x2end, x3start)) +
  scale_x_continuous(breaks = c(1950000, 69300000, 69400000, 129200000, 129300000, 129400000)) + 
  scale_y_continuous(limits = c(0.5, 6.5), breaks = NULL) +
  scale_color_manual(values = c("red", "blue", "green", "purple", "orange"))

png(filename = "Chr11.png", width = 1800, height = 400)
print(region)
dev.off()
################################################################################################################################################################################

rm(list = ls())

BiRS = read.csv('C50Breast_BiRS.csv')
SCAN = read.csv('C50Breast_Scan.csv')
KSCN = read.csv('C50Breast_KSGWAS.csv')
SSSS = read.csv('C50Breast_4S.csv')
REGN = read.csv('C50Breast_Regenie.csv')
chr12_gene = read.csv('C50Breast_chr12_gene.csv')

BiRS = BiRS[, -1]; SCAN = SCAN[, -1]; KSCN = KSCN[, -1]; SSSS = SSSS[, -1]; REGN = REGN[, -1]

BiRSChr12 = BiRS[which(BiRS$chrsel == 12), ]
SCANChr12 = SCAN[which(SCAN$chrsel == 12), ]
KSCNChr12 = KSCN[which(KSCN$chr == 12), ]
SSSSChr12 = SSSS[which(SSSS$chr == 12), ]
REGNChr12 = REGN[which(REGN$CHR == 12), ]

chr12_len = nrow(chr12_gene)
startposchr12 = chr12_gene$Start
endposchr12 = chr12_gene$End
deschr12 = chr12_gene$Gene

x1start = 28000000; x1end = 28800000; 
x2start = 95900000; x2end = 96200000; 
x3start = 115000000; x3end = 1168000000;

gene_sel = c('PTHLH', 'LOC729291', 'CCDC91', 'NTN4', 'TBX3', 'MED13L')
xgene = NULL
textpos = rep(0, length(gene_sel))
for (i in 1:length(gene_sel))
{
  ind_seli = which(deschr12 == gene_sel[i])
  xgene = c(xgene, startposchr12[ind_seli]:endposchr12[ind_seli])
  textpos[i] = (startposchr12[ind_seli] + endposchr12[ind_seli])/2
}

notation = c(rep('BiRS-DCF', length(BiRSChr12$possel)), rep("Q-SCAN", length(SCANChr12$possel)), rep('KnockoffScreen', length(KSCNChr12$bp)), rep('4S', length(SSSSChr12$possel)), rep("REGENIE", length(REGNChr12$BP)))
notation = factor(notation, levels = c('BiRS-DCF', 'Q-SCAN', 'KnockoffScreen', '4S', 'REGENIE'))
X = c(BiRSChr12$possel, SCANChr12$possel, KSCNChr12$bp, SSSSChr12$possel, REGNChr12$BP)
Y = c(rep(6, length(BiRSChr12$possel)), rep(5, length(SCANChr12$possel)), 
      rep(4, length(KSCNChr12$bp)), rep(3, length(SSSSChr12$possel)), rep(2, length(REGNChr12$BP))) 

CHR12Data = data.frame(X, Y, notation)
Pos = rep(1, length(xgene))
GeneData = data.frame(xgene, Pos)

region = ggplot(data = CHR12Data, aes(x = X, y = Y)) +
  geom_point(stat = 'identity', aes(color = notation), size = 4) + 
  geom_point(data = GeneData, mapping = aes(x = xgene, y = Pos), size = 4) + 
  ggtitle("C.Location of Significant Regions in Chromosome 12") +
  xlab("Genomic Location(kb)") + 
  ylab(" ") + 
  labs(color = " ")

for(i in 1:length(gene_sel))
{
  region = region + annotate("text", x = textpos[i], y = 0.6, label = gene_sel[i], 
                             size = 5, fontface = 'bold.italic')
}

region = region + theme_gray() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(axis.text.x.top = element_blank(), axis.ticks.x.top = element_blank()) + 
  theme(plot.title = element_text(face = "bold", size = 24), axis.title.x = element_text(face = "bold", size = 20), 
        legend.text = element_text(face = "bold", size = 17), axis.text.x = element_text(face = "bold", size = 17), 
        legend.position = "none") +
  scale_x_break(c(x1end, x2start)) +
  scale_x_break(c(x2end, x3start)) +
  scale_x_continuous(breaks = c(28100000, 28300000, 28500000, 28700000, 96050000, 115000000, 115500000, 116000000, 116500000)) + 
  scale_y_continuous(limits = c(0.5, 6.5), breaks = NULL) +
  scale_color_manual(values = c("red", "blue", "green", "purple", "orange"))

png(filename = "Chr12.png", width = 1800, height = 400)
print(region)
dev.off()
################################################################################################################################################################################

rm(list = ls())

BiRS = read.csv('C50Breast_BiRS.csv')
SCAN = read.csv('C50Breast_Scan.csv')
KSCN = read.csv('C50Breast_KSGWAS.csv')
SSSS = read.csv('C50Breast_4S.csv')
REGN = read.csv('C50Breast_Regenie.csv')
chr19_gene = read.csv('C50Breast_chr19_gene.csv')

BiRS = BiRS[, -1]; SCAN = SCAN[, -1]; KSCN = KSCN[, -1]; REGN = REGN[, -1]

BiRSChr19 = BiRS[which(BiRS$chrsel == 19), ]
SCANChr19 = SCAN[which(SCAN$chrsel == 19), ]
KSCNChr19 = KSCN[which(KSCN$chr == 19), ]
SSSSChr19 = KSCN[which(KSCN$chr == 19), ]
REGNChr19 = REGN[which(REGN$CHR == 19), ]

chr19_len = nrow(chr19_gene)
startposchr19 = chr19_gene$Start
endposchr19 = chr19_gene$End
deschr19 = chr19_gene$Gene

x1start = 44270000; x1end = 44400000; 

gene_sel = c('KCNN4', 'LYPD5', 'ZNF283', 'ZNF404')
xgene = NULL
textpos = rep(0, length(gene_sel))
for (i in 1:length(gene_sel))
{
  ind_seli = which(deschr19 == gene_sel[i])
  xgene = c(xgene, startposchr19[ind_seli]:endposchr19[ind_seli])
  textpos[i] = (startposchr19[ind_seli] + endposchr19[ind_seli])/2
}

notation = c(rep('BiRS-DCF', length(BiRSChr19$possel)), rep("Q-SCAN", length(SCANChr19$possel)), 
             rep('KnockoffScreen', length(KSCNChr19$bp)), rep('4S', length(SSSSChr19$possel)),  rep("REGENIE", length(REGNChr19$BP)))
notation = factor(notation, levels = c('BiRS-DCF', 'Q-SCAN', 'KnockoffScreen', '4S', 'REGENIE'))
X = c(BiRSChr19$possel, SCANChr19$possel, KSCNChr19$bp, SSSSChr19$possel, REGNChr19$BP)
Y = c(rep(6, length(BiRSChr19$possel)), rep(5, length(SCANChr19$possel)), 
      rep(4, length(KSCNChr19$bp)), rep(3, length(SSSSChr19$possel)), rep(2, length(REGNChr19$BP))) 

CHR19Data = data.frame(X, Y, notation)
Pos = rep(1, length(xgene))
GeneData = data.frame(xgene, Pos)

region = ggplot(data = CHR19Data, aes(x = X, y = Y)) +
  geom_point(stat = 'identity', aes(color = notation), size = 4) + 
  geom_point(data = GeneData, mapping = aes(x = xgene, y = Pos), size = 4) + 
  ggtitle("A.Location of Significant Regions in Chromosome 19") +
  xlab("Genomic Location(kb)") + 
  ylab(" ") + 
  labs(color = " ")

for(i in 1:length(gene_sel))
{
  region = region + annotate("text", x = textpos[i], y = 0.6, label = gene_sel[i], 
                             size = 5, fontface = 'bold.italic')
}

region = region + theme_gray() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  theme(axis.text.x.top = element_blank(), axis.ticks.x.top = element_blank()) + 
  theme(plot.title = element_text(face = "bold", size = 24), axis.title.x = element_text(face = "bold", size = 20), 
        legend.text = element_text(face = "bold", size = 17), axis.text.x = element_text(face = "bold", size = 17), 
        legend.position = "none") +
  scale_x_continuous(breaks = c(44270000, 44290000, 44310000, 44330000, 44350000, 44370000, 44390000)) + 
  scale_y_continuous(limits = c(0.5, 6.5), breaks = NULL) +
  scale_color_manual(values = c("red", "blue", "green", "purple", "orange"))

png(filename = "Chr19.png", width = 1800, height = 400)
print(region)
dev.off()

# notation = rep('BiRS-DCF', length(BiRSChr19$possel))
# X = BiRSChr19$possel
# Y = rep(5, length(BiRSChr19$possel)) 
# 
# CHR19Data = data.frame(X, Y, notation)
# Pos = rep(1, length(xgene))
# GeneData = data.frame(xgene, Pos)
# 
# region = ggplot(data = CHR19Data, aes(x = X, y = Y)) +
#   geom_point(stat = 'identity', aes(color = 'red'), size = 3) + 
#   geom_point(data = GeneData, mapping = aes(x = xgene, y = Pos), size = 3) + 
#   ggtitle("Location of Significant Regions in Chromosome 19") +
#   xlab("Genomic Location(kb)") + 
#   ylab(" ") + 
#   labs(color = " ")
# 
# for(i in 1:length(gene_sel))
# {
#   region = region + annotate("text", x = textpos[i], y = 0.8, label = gene_sel[i], 
#                              size = 5, fontface = 'bold.italic')
# }
# 
# region = region + theme_gray() +
#   theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
#   theme(axis.text.x.top = element_blank(), axis.ticks.x.top = element_blank()) + 
#   theme(plot.title = element_text(face = "bold", size = 22), axis.title.x = element_text(face = "bold", size = 18), 
#         legend.text = element_text(face = "bold", size = 15), axis.text.x = element_text(face = "bold", size = 15), 
#         legend.position = "right") +
#   scale_y_continuous(limits = c(0.5, 6), breaks = NULL) +
#   scale_color_manual(values = c("red"))
# 
# png(filename = "Chr19.png", width = 1600, height = 600)
# print(region)
# dev.off()