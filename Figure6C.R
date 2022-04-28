## Figure 6C
## 'The evolution of two transmissible cancers in Tasmanian devils'
## Stammnitz et al., 2022
## maxrupsta@gmail

library(readxl)

## set input path(s)
setwd('/Users/ms37/Desktop/Labwork/DFT_evolution/doc/manuscripts/The Evolutionary History of Two Transmissible Cancers in Tasmanian Devils/Tables/v6/')


# DFT1 and DFT2 dNdS summary #
##############################

dNdS.summary <- as.matrix(read_xlsx('/Users/mstammnitz/Desktop/DFT_evolution/doc/manuscripts/The Evolutionary History of Two Transmissible Cancers in Tasmanian Devils/Tables/v6/Table-S8_v6.xlsx', sheet = 1))
colnames(dNdS.summary) <- dNdS.summary[2,]
dNdS.summary <- dNdS.summary[-c(1:2),]
dNdS.summary[,"Genic variants considered"] <- as.numeric(dNdS.summary[,"Genic variants considered"])
rownames(dNdS.summary) <- dNdS.summary[,1]
dNdS.summary <- dNdS.summary[,-1]
dNdS.summary <- dNdS.summary[-3,-1]
class(dNdS.summary) <- 'numeric'

pdf('Figure6C_DFT1_DFT2_dNdS.pdf', 
    width = 7, height = 10)

par(mar=c(7, 17, 5, 4))
fits <- barplot(t(dNdS.summary[,1]), 
                beside = T,
                col = NA,
                xaxt = 'n',
                yaxt = 'n',
                ylab = '',
                border = 'white',
                las = 1,
                ylim = c(0,2))

## neutral line
abline(h = 1, lty = 2, lwd = 2)

## error bars
arrows(fits[1,1],
       y1 = t(dNdS.summary[1,2]), 
       y0 = t(dNdS.summary[1,3]),
       angle = 90, code = 3, length = 0.1, lwd = 2, col = 'black')
arrows(fits[1,2],
       y1 = t(dNdS.summary[2,2]), 
       y0 = t(dNdS.summary[2,3]),
       angle = 90, code = 3, length = 0.1, lwd = 2, col = 'black')

## points
points(x = fits[1,],
       y = t(dNdS.summary[,1]),
       col = c('cornflowerblue', 'red'),
       cex = 4,
       pch = 16)

## axis
axis(side = 1, at = apply(fits,2,mean)[1], tick = FALSE,
     labels = 'DFT1',
     cex.axis = 3, padj = 0.8, col.axis = 'cornflowerblue')
axis(side = 1, at = apply(fits,2,mean)[1], tick = FALSE,
     labels = '(N = 2,102)',
     cex.axis = 1.8, padj = 0.8, col.axis = 'cornflowerblue', line = 2.5)
axis(side = 1, at = apply(fits,2,mean)[2], tick = FALSE,
     labels = 'DFT2',
     cex.axis = 3, padj = 0.8, col.axis = 'red')
axis(side = 1, at = apply(fits,2,mean)[2], tick = FALSE,
     labels = '(N = 183)',
     cex.axis = 1.8, padj = 0.8, col.axis = 'red', line = 2.5)
axis(2, at = seq(f = 0, t = 2, length.out = 5),
     labels = c('0.0', '0.5', '1.0', '1.5', '2'),
     cex.axis = 2, las = 2, lwd = 3, hadj = 1.2, line = 2)
title(ylab = "dN/dS", line = 8, cex.lab = 3)

dev.off()

## clean up environment
rm(list=ls())