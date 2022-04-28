## Figure 6B
## 'The evolution of two transmissible cancers in Tasmanian devils'
## Stammnitz et al., 2022
## maxrupsta@gmail

library(readxl)

## set input path(s)
setwd('/Users/ms37/Desktop/Labwork/DFT_evolution/doc/manuscripts/The Evolutionary History of Two Transmissible Cancers in Tasmanian Devils/Tables/v6/')


# LZTR1 structural variant display #
####################################

## load LZTR1 annotation
load("/Users/mstammnitz/Desktop/DFT_evolution/doc/manuscripts/The Evolutionary History of Two Transmissible Cancers in Tasmanian Devils/Tables/v6/Supplementary_data/Sarcophilus_harrisii.mSarHar1.11.102.gtf.Rdata")
LZTR1.ensembl.out <- ensembl[grep('LZTR1', ensembl[,9]),]
LZTR1.ensembl.out.exons <- LZTR1.ensembl.out[which(LZTR1.ensembl.out[,3] == 'exon'),]
LZTR1.ensembl.out.exons <- LZTR1.ensembl.out.exons[grep('LZTR1-201', LZTR1.ensembl.out.exons[,9]),]

## load LZTR1 structural variants
DFT1.SVs <- as.matrix(read_xlsx('/Users/mstammnitz/Desktop/DFT_evolution/doc/manuscripts/The Evolutionary History of Two Transmissible Cancers in Tasmanian Devils/Tables/v6/Table-S5_v6.xlsx', sheet = 2))
colnames(DFT1.SVs) <- as.character(DFT1.SVs[2,])
DFT1.SVs <- DFT1.SVs[-c(1:2),]
DFT1.SVs.LZTR1 <- DFT1.SVs[grep('LZTR1', DFT1.SVs[,'ENSEMBL ANNOTATION2']),]

pdf("Figure6B_LZTR1_mutations.pdf",
    width = 16, height = 5)
par(mar = c(8,1,1,1))

cexs <- 2

plot(x = c(712247000, 712263000),
     y = rep(1,2), col = 'white',
     yaxt = 'n', frame = F, ylim = c(0,2), ylab = '', 
     xaxt = 'n', xlab = '', xlim = c(712246000, 712263000),
     main = '', cex.main = 2)

## add exon polygons and intron zick-zacks
for(i in 1:nrow(LZTR1.ensembl.out.exons)){
  
  ## exon "squares"
  polygon(x = as.numeric(LZTR1.ensembl.out.exons[i,c(4,5,5,4)]),
          y = c(0.5,0.5,0.75,0.75), 
          col = 'black',
          border = 'black')
  
  ## intron zick-zack
  if(i != nrow(LZTR1.ensembl.out.exons)){
    
    tmp.intron.middle <- mean(c(as.numeric(LZTR1.ensembl.out.exons[i,4]),
                                as.numeric(LZTR1.ensembl.out.exons[c(i+1),5])))
    
    lines(x = c(as.numeric(LZTR1.ensembl.out.exons[i,4]),
                tmp.intron.middle),
          y = c(0.625, 0.5),
          col = 'black')
    
    lines(x = c(tmp.intron.middle, as.numeric(LZTR1.ensembl.out.exons[c(i+1),5])),
          y = c(0.5, 0.625),
          col = 'black') 
    
  }
  
}

## display natural ATG
text(x = 712262146,
     y = 0.15, 
     labels = "ATG", 
     cex = cexs,
     pos = 4,
     offset = -0.1)
lines(x = c(712262146,712262146),
      y = c(0.3,0.5),
      col = 'black')

## display natural TGA
text(x = 712249072,
     y = 0.15, 
     labels = "TGA", 
     cex = cexs,
     pos = 2,
     offset = -0.1)
lines(x = c(712249072,712249072),
      y = c(0.3,0.5),
      col = 'black')

## display SV breakpoints
DFT1.SVs.LZTR1 <- matrix(NA, nrow = 6, ncol = 5)
colnames(DFT1.SVs.LZTR1) <- c('POS1', 'GENE1', 'POS2', 'GENE2', 'STRAND2')
DFT1.SVs.LZTR1[,'POS1'] <- c(455942825,
                             668980987,
                             706692959,
                             707843055,
                             708882343,
                             709003173)
DFT1.SVs.LZTR1[,'GENE1'] <- c('-',
                              'KCNN1',
                              'ASCC2',
                              'OSBP2',
                              '-',
                              'SLC5A1')
DFT1.SVs.LZTR1[,'POS2'] <- c(mean(c(712256122,712256224)),
                             mean(c(712256122,712256224)),
                             mean(c(712249364,712249363)),
                             mean(c(712256634,712256495)),
                             mean(c(712256634,712256495)),
                             mean(c(712249364,712249363)))
DFT1.SVs.LZTR1[,'GENE2'] <- c('LZTR1',
                              'LZTR1',
                              'LZTR1',
                              'LZTR1',
                              'LZTR1',
                              'LZTR1')
DFT1.SVs.LZTR1[,'STRAND2'] <- c('+',
                                '-',
                                '-',
                                '-',
                                '+',
                                '+')

for (i in c(1,3,5)){
  
    lines(x = rep(as.numeric(DFT1.SVs.LZTR1[i,'POS2']),2),
          y = c(0.375,1),
          col = 'cornflowerblue',
          lwd = 2,
          lty = 1)
    
    points(x = as.numeric(DFT1.SVs.LZTR1[i,'POS2']),
           y = 0.375,
           col = 'cornflowerblue', 
           bg = 'cornflowerblue', 
           cex = cexs*1.5, 
           pch = 25)

}

## display mutation #1

### local (intra-window) SVs
fit.func <- function(x){
  y <- fit$coefficients[3]*x^2+c(fit$coefficients[2]*x)+fit$coefficients[1]
  return(y)
}
bp.x <- c(as.numeric(DFT1.SVs.LZTR1[1,'POS2']), as.numeric(DFT1.SVs.LZTR1[1,'POS2']) + 5000, as.numeric(DFT1.SVs.LZTR1[1,'POS2']) + 10000)
bp.y <- c(1, 2.05, 1)
fit <- glm(bp.y ~ poly(bp.x, 2, raw = T))
x.fit.out <- seq(f = bp.x[1], t = bp.x[3], length.out = 5000)
y.fit.out <- fit.func(x.fit.out)
lines(x = x.fit.out[1:1000], 
      y = y.fit.out[1:1000], 
      col = 'cornflowerblue', 
      lwd = 2)
text(x = 712258171,
     y = 1.8,
     labels = 'Chr1:455.9',
     col = 'cornflowerblue',
     cex = 1.8)

## display mutation #2
bp.x <- c(as.numeric(DFT1.SVs.LZTR1[2,'POS2']) - 10000, as.numeric(DFT1.SVs.LZTR1[2,'POS2']) - 5000, as.numeric(DFT1.SVs.LZTR1[2,'POS2']))
bp.y <- c(1, 2.05, 1)
fit <- glm(bp.y ~ poly(bp.x, 2, raw = T))
x.fit.out <- seq(f = bp.x[1], t = bp.x[3], length.out = 5000)
y.fit.out <- fit.func(x.fit.out)
lines(x = x.fit.out[4500:5000], 
      y = y.fit.out[4500:5000], 
      col = 'cornflowerblue', 
      lwd = 2)
text(x = 712255173,
     y = 1.4,
     labels = 'Chr1:669.0',
     col = 'cornflowerblue',
     cex = 1.8,
     pos = 2)

## display mutation #3
bp.x <- c(as.numeric(DFT1.SVs.LZTR1[3,'POS2']) - 10000, as.numeric(DFT1.SVs.LZTR1[3,'POS2']) - 5000, as.numeric(DFT1.SVs.LZTR1[3,'POS2']))
bp.y <- c(1, 2.05, 1)
fit <- glm(bp.y ~ poly(bp.x, 2, raw = T))
x.fit.out <- seq(f = bp.x[1], t = bp.x[3], length.out = 5000)
y.fit.out <- fit.func(x.fit.out)
lines(x = x.fit.out[4000:5000], 
      y = y.fit.out[4000:5000], 
      col = 'cornflowerblue', 
      lwd = 2)
text(x = 712247363,
     y = 1.8,
     labels = 'Chr1:706.7',
     col = 'cornflowerblue',
     cex = 1.8)

## display mutation #4
bp.x <- c(as.numeric(DFT1.SVs.LZTR1[4,'POS2']) - 10000, as.numeric(DFT1.SVs.LZTR1[4,'POS2']) - 5000, as.numeric(DFT1.SVs.LZTR1[4,'POS2']))
bp.y <- c(1, 2.05, 1)
fit <- glm(bp.y ~ poly(bp.x, 2, raw = T))
x.fit.out <- seq(f = bp.x[1], t = bp.x[3], length.out = 5000)
y.fit.out <- fit.func(x.fit.out)
lines(x = x.fit.out[4000:5000], 
      y = y.fit.out[4000:5000], 
      col = 'cornflowerblue',
      lwd = 2)
text(x = 712254564,
     y = 1.8,
     labels = 'Chr1:707.8',
     col = 'cornflowerblue',
     cex = 1.8)

## display mutation #5
bp.x <- c(as.numeric(DFT1.SVs.LZTR1[5,'POS2']), as.numeric(DFT1.SVs.LZTR1[5,'POS2']) + 5000, as.numeric(DFT1.SVs.LZTR1[5,'POS2']) + 10000)
bp.y <- c(1, 2.05, 1)
fit <- glm(bp.y ~ poly(bp.x, 2, raw = T))
x.fit.out <- seq(f = bp.x[1], t = bp.x[3], length.out = 5000)
y.fit.out <- fit.func(x.fit.out)
lines(x = x.fit.out[1:500], 
      y = y.fit.out[1:500], 
      col = 'cornflowerblue', 
      lwd = 2)
text(x = 712257563,
     y = 1.4,
     labels = 'Chr1:708.9',
     col = 'cornflowerblue',
     cex = 1.8,
     pos = 4)

## display mutation #6
bp.x <- c(as.numeric(DFT1.SVs.LZTR1[6,'POS2']), as.numeric(DFT1.SVs.LZTR1[6,'POS2']) + 5000, as.numeric(DFT1.SVs.LZTR1[6,'POS2']) + 10000)
bp.y <- c(1, 2.05, 1)
fit <- glm(bp.y ~ poly(bp.x, 2, raw = T))
x.fit.out <- seq(f = bp.x[1], t = bp.x[3], length.out = 5000)
y.fit.out <- fit.func(x.fit.out)
lines(x = x.fit.out[1:1000], 
      y = y.fit.out[1:1000], 
      col = 'cornflowerblue', 
      lwd = 2)
text(x = 712251362,
     y = 1.8,
     labels = 'Chr1:709.0',
     col = 'cornflowerblue',
     cex = 1.8)

## axis
axis(side = 1, 
     at = seq(f = 712247000, t = 712263000, length.out = 5),
     labels = c("712.247", "712.251", "712.255", "712.259", "712.263"),
     cex.axis = 1.8)
mtext(substitute(paste(italic("LZTR1"), " locus on Chromosome 1 [MBP]")), side = 1, line = 5, cex = 2.5)
dev.off()

## clean up environment
rm(list=ls())