## Figure 2E
## 'The evolution of two transmissible cancers in Tasmanian devils'
## Stammnitz et al., 2022
## maxrupsta@gmail

library(scales)

## set input path(s)
setwd('/Users/ms37/Desktop/Labwork/DFT_evolution/doc/manuscripts/The Evolutionary History of Two Transmissible Cancers in Tasmanian Devils/Tables/v6/')


# DFT1-C subclonal VAF histograms #
###################################

## load corrected VAFs
load("/Users/mstammnitz/Desktop/DFT_evolution/doc/manuscripts/The Evolutionary History of Two Transmissible Cancers in Tasmanian Devils/Tables/v6/Supplementary_data/DFT2_subclones_corrected_VAF.Rdata")
VAF.139T4 <- DFT1.corrected.SNVs[,'139T4']
VAF.139T1 <- DFT1.corrected.SNVs[,'139T1']
VAF.140T <- DFT1.corrected.SNVs[,'140T']

## split VAFs into three separate categories

### SC1
SC1.id <- as.numeric(which(c(rowSums(DFT1.corrected.SNVs.bin) == 4 & 
                               DFT1.corrected.SNVs.bin[,'139T1'] == 1 & 
                               DFT1.corrected.SNVs.bin[,'139T4'] == 1 & 
                               DFT1.corrected.SNVs.bin[,'139T5'] == 1 & 
                               DFT1.corrected.SNVs.bin[,'139T6'] == 1) |
                             c(rowSums(DFT1.corrected.SNVs.bin) == 3 & 
                                 DFT1.corrected.SNVs.bin[,'139T4'] == 1 & 
                                 DFT1.corrected.SNVs.bin[,'139T5'] == 1 & 
                                 DFT1.corrected.SNVs.bin[,'139T6'] == 1)))

SC1.VAF.139T4 <- VAF.139T4[SC1.id]
SC1.VAF.139T1 <- VAF.139T1[SC1.id]
SC1.VAF.140T <- VAF.140T[SC1.id]

### SC2
SC2.id <- as.numeric(which(c(rowSums(DFT1.corrected.SNVs.bin) == 4 & 
                               DFT1.corrected.SNVs.bin[,'139T1'] == 1 & 
                               DFT1.corrected.SNVs.bin[,'140T'] == 1 & 
                               DFT1.corrected.SNVs.bin[,'141T'] == 1 & 
                               DFT1.corrected.SNVs.bin[,'142T'] == 1)))
SC2.VAF.139T4 <- VAF.139T4[SC2.id]
SC2.VAF.139T1 <- VAF.139T1[SC2.id]
SC2.VAF.140T <- VAF.140T[SC2.id]

### plot histograms and/or kernels
out.sc1.139T4 <- density(SC1.VAF.139T4, n = 500, kernel = 'gaussian', from = -0.1, to = 1.10)
out.sc2.139T4 <- density(SC2.VAF.139T4, n = 500, kernel = 'gaussian', from = -0.1, to = 1.10)
out.sc1.139T4$y <- scales::rescale(x = out.sc1.139T4$y, to = c(0, max(out.sc1.139T4$y,out.sc2.139T4$y)))
out.sc2.139T4$y <- scales::rescale(x = out.sc2.139T4$y, to = c(0, max(out.sc1.139T4$y,out.sc2.139T4$y)))
out.sc1.139T1 <- density(SC1.VAF.139T1, n = 500, kernel = 'gaussian', from = -0.1, to = 1.10)
out.sc2.139T1 <- density(SC2.VAF.139T1, n = 500, kernel = 'gaussian', from = -0.1, to = 1.10)
out.sc1.139T1$y <- scales::rescale(x = out.sc1.139T1$y, to = c(0, max(out.sc1.139T1$y,out.sc2.139T1$y)))
out.sc2.139T1$y <- scales::rescale(x = out.sc2.139T1$y, to = c(0, max(out.sc1.139T1$y,out.sc2.139T1$y)))
out.sc1.140T <- density(SC1.VAF.140T, n = 500, kernel = 'gaussian', from = -0.1, to = 1.10)
out.sc2.140T <- density(SC2.VAF.140T, n = 500, kernel = 'gaussian', from = -0.1, to = 1.10)
out.sc1.140T$y <- scales::rescale(x = out.sc1.140T$y, to = c(0, max(out.sc1.140T$y,out.sc2.140T$y)))
out.sc2.140T$y <- scales::rescale(x = out.sc2.140T$y, to = c(0, max(out.sc1.140T$y,out.sc2.140T$y)))

## SC colors
SC1.col <- 'grey1'
SC2.col <- 'grey80'
transp <- 0.3

### Plot
pdf('Figure2E_VAF_histograms.pdf',
    height = 20, width = 12)
par(mar = c(16.5,3,0,3), mfcol = c(3,1))

#### 139T4
plot(out.sc1.139T4, 
     xlim = c(0,1), main = '',
     bty = 'n', col = NA, 
     yaxt = 'n', xlab = '', ylab = '', lwd = 3,
     xaxt = 'n')
polygon(out.sc1.139T4, col = alpha(SC1.col, transp), border = "black")
polygon(out.sc2.139T4, col = alpha(SC2.col, transp), border = "black")
axis(1, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c(0, 20, 40, 60, 80, 100), padj = 1,
     cex.axis = 7, lwd = 3)

#### 139T1
plot(out.sc1.139T1, 
     xlim = c(0,1), main = '',
     bty = 'n', col = NA, 
     yaxt = 'n', xlab = '', ylab = '', lwd = 3,
     xaxt = 'n')
polygon(out.sc1.139T1, col = alpha(SC1.col, transp), border = "black",)
polygon(out.sc2.139T1, col = alpha(SC2.col, transp), border = "black")
axis(1, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c(0, 20, 40, 60, 80, 100), padj = 1,
     cex.axis = 7, lwd = 3)

#### 140T
plot(out.sc1.140T, 
     xlim = c(0,1), main = '',
     bty = 'n', col = NA, 
     yaxt = 'n', xlab = '', 
     ylab = '', lwd = 3,
     xaxt = 'n')
polygon(out.sc1.140T, col = alpha(SC1.col, transp), border = "black")
polygon(out.sc2.140T, col = alpha(SC2.col, transp), border = "black")
axis(1, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c(0, 20, 40, 60, 80, 100), padj = 1,
     cex.axis = 7, lwd = 3)

mtext('Variant allele fraction (%)', side = 1, line = 13.7, cex = 6)
dev.off()