## Figure 2B
## 'The evolution of two transmissible cancers in Tasmanian devils'
## Stammnitz et al., 2022
## maxrupsta@gmail

library(scales)

## set input path(s)
setwd('/Tables')


# DFT2-B subclonal VAF histograms #
###################################

## load corrected VAFs
load("DFT2_subclones_corrected_VAF.Rdata")
VAF.1334T1 <- DFT2.corrected.SNVs[,'1334T1']
VAF.1509T1 <- DFT2.corrected.SNVs[,'1509T1']
VAF.1529T2 <- DFT2.corrected.SNVs[,'1529T2']

## split VAFs into three separate categories

### SC1
cladeA <- c("202T1", "202T2", "203T1", "203T2", "203T3", "338T", "339T", 
            "637T1", "637T2", "638T1", "809T1", "812T1", "818T3", "1515T1",
            "1524T1", "1525T1", "1528T", "1538T1", "1538T2", "1538T3", "1548T1",
            "1553T1", "1553T2")
SC1.id <- as.numeric(which(c(rowSums(DFT2.corrected.SNVs.bin) == 12 & rowSums(DFT2.corrected.SNVs.bin[,colnames(DFT2.corrected.SNVs.bin) %in% cladeA]) == 0 & DFT2.corrected.SNVs.bin[,'1509T1'] == 1) |
                             c(rowSums(DFT2.corrected.SNVs.bin) == 13 & rowSums(DFT2.corrected.SNVs.bin[,colnames(DFT2.corrected.SNVs.bin) %in% cladeA]) == 0 & DFT2.corrected.SNVs.bin[,'1509T1'] == 1) |
                             c(rowSums(DFT2.corrected.SNVs.bin) == 14 & rowSums(DFT2.corrected.SNVs.bin[,colnames(DFT2.corrected.SNVs.bin) %in% cladeA]) == 0 & DFT2.corrected.SNVs.bin[,'1509T1'] == 1) |
                             c(rowSums(DFT2.corrected.SNVs.bin) == 15 & rowSums(DFT2.corrected.SNVs.bin[,colnames(DFT2.corrected.SNVs.bin) %in% cladeA]) == 0 & DFT2.corrected.SNVs.bin[,'1509T1'] == 1)))
SC1.VAF.1529T2 <- VAF.1529T2[SC1.id]
SC1.VAF.1509T1 <- VAF.1509T1[SC1.id]
SC1.VAF.1334T1 <- VAF.1334T1[SC1.id]

### SC2
SC2.id <- as.numeric(which(rowSums(DFT2.corrected.SNVs.bin) == 4 & DFT2.corrected.SNVs.bin[,'1509T1'] == 1 & DFT2.corrected.SNVs.bin[,'1509T2'] == 1 & DFT2.corrected.SNVs.bin[,'1334T1'] == 1 & DFT2.corrected.SNVs.bin[,'1545T1'] == 1))
SC2.VAF.1529T2 <- VAF.1529T2[SC2.id]
SC2.VAF.1509T1 <- VAF.1509T1[SC2.id]
SC2.VAF.1334T1 <- VAF.1334T1[SC2.id]

### plot histograms and/or kernels
out.sc1.1334T1 <- density(SC1.VAF.1334T1, n = 500, kernel = 'gaussian', from = -0.1, to = 1.10)
out.sc2.1334T1<- density(SC2.VAF.1334T1, n = 500, kernel = 'gaussian', from = -0.1, to = 1.10)
out.sc1.1334T1$y <- scales::rescale(x = out.sc1.1334T1$y, to = c(0, max(out.sc1.1334T1$y,out.sc2.1334T1$y)))
out.sc2.1334T1$y <- scales::rescale(x = out.sc2.1334T1$y, to = c(0, max(out.sc1.1334T1$y,out.sc2.1334T1$y)))
out.sc1.1509T1 <- density(SC1.VAF.1509T1, n = 500, kernel = 'gaussian', from = -0.1, to = 1.10)
out.sc2.1509T1 <- density(SC2.VAF.1509T1, n = 500, kernel = 'gaussian', from = -0.1, to = 1.10)
out.sc1.1509T1$y <- scales::rescale(x = out.sc1.1509T1$y, to = c(0, max(out.sc1.1509T1$y,out.sc2.1509T1$y)))
out.sc2.1509T1$y <- scales::rescale(x = out.sc2.1509T1$y, to = c(0, max(out.sc1.1509T1$y,out.sc2.1509T1$y)))
out.sc1.1529T2 <- density(SC1.VAF.1529T2, n = 500, kernel = 'gaussian', from = -0.1, to = 1.10)
out.sc2.1529T2 <- density(SC2.VAF.1529T2, n = 500, kernel = 'gaussian', from = -0.1, to = 1.10)
out.sc1.1529T2$y <- scales::rescale(x = out.sc1.1529T2$y, to = c(0, max(out.sc1.1529T2$y,out.sc2.1529T2$y)))
out.sc2.1529T2$y <- scales::rescale(x = out.sc2.1529T2$y, to = c(0, max(out.sc1.1529T2$y,out.sc2.1529T2$y)))

## SC colors
SC1.col <- 'grey1'
SC2.col <- 'grey80'
transp <- 0.3

### Plot
pdf('Figure2B_VAF_histograms.pdf',
    height = 20, width = 12)
par(mar = c(16.5,3,0,3), mfcol = c(3,1))

#### 1529T2
plot(out.sc1.1529T2, 
     xlim = c(0,1), main = '',
     bty = 'n', col = NA, 
     yaxt = 'n', xlab = '', ylab = '', lwd = 3,
     xaxt = 'n')
polygon(out.sc1.1529T2, col = alpha(SC1.col, transp), border="black")
polygon(out.sc2.1529T2, col = alpha(SC2.col, transp), border="black")
axis(1, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c(0, 20, 40, 60, 80, 100), padj = 1,
     cex.axis = 7, lwd = 3)

#### 1509T1
plot(out.sc1.1509T1, 
     xlim = c(0,1), main = '',
     bty = 'n', col = NA, 
     yaxt = 'n', xlab = '', ylab = '', lwd = 3,
     xaxt = 'n')
axis(1, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c(0, 20, 40, 60, 80, 100), padj = 1,
     cex.axis = 7, lwd = 3)
polygon(out.sc1.1509T1, col = alpha(SC1.col, transp), border = "black",)
polygon(out.sc2.1509T1, col = alpha(SC2.col, transp), border = "black")

#### 1334T1
plot(out.sc1.1334T1, 
     xlim = c(0,1), main = '',
     bty = 'n', col = NA, 
     yaxt = 'n', xlab = '', 
     ylab = '', lwd = 3,
     xaxt = 'n')
axis(1, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c(0, 20, 40, 60, 80, 100), padj = 1,
     cex.axis = 7, lwd = 3)
polygon(out.sc1.1334T1, col = alpha(SC1.col, transp), border = "black")
polygon(out.sc2.1334T1, col = alpha(SC2.col, transp), border = "black")

mtext('Variant allele fraction (%)', side = 1, line = 13.7, cex = 6)
dev.off()
