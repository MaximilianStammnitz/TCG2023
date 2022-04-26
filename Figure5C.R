## Figure 5C
## 'The evolution of two transmissible cancers in Tasmanian devils'
## Stammnitz et al., 2022
## maxrupsta@gmail

library(readxl)
library(stringr)
library(GenomicRanges)
library(circlize)
library(scales)

## set input path(s)
setwd('/Users/ms37/Desktop/Labwork/DFT_evolution/doc/manuscripts/The Evolutionary History of Two Transmissible Cancers in Tasmanian Devils/Tables/v6/')


# Chromothripsis in 384T1 - circos plot #
#########################################

## load and subset 384T1 SVs
DFT1.SVs <- as.matrix(read_xlsx('/Users/mstammnitz/Desktop/DFT_evolution/doc/manuscripts/The Evolutionary History of Two Transmissible Cancers in Tasmanian Devils/Tables/v6/Table-S5_v6.xlsx', sheet = 2))
colnames(DFT1.SVs) <- as.character(DFT1.SVs[2,])
DFT1.SVs <- DFT1.SVs[-c(1:2),]
DFT1.SVs.384T1.MSG <- DFT1.SVs[grep('MSG', DFT1.SVs[,'CALLER']),]
DFT1.SVs.384T1.MSG <- DFT1.SVs.384T1.MSG[which(apply(DFT1.SVs.384T1.MSG[,grep('384T1', colnames(DFT1.SVs.384T1.MSG)),drop=F], 1, function(x){out <- as.numeric(str_split_fixed(x, '/', 2)[,1]); all(out >= 1)}) == T),]
DFT1.SVs.384T1.SvABA <- DFT1.SVs[-grep('MSG', DFT1.SVs[,'CALLER']),]
DFT1.SVs.384T1.SvABA <- DFT1.SVs.384T1.SvABA[which(apply(DFT1.SVs.384T1.SvABA[,grep('384T1', colnames(DFT1.SVs.384T1.SvABA)),drop=F], 1, function(x){out <- as.numeric(x); all(out >= 3)}) == T),]

### 384T1 SVs (all non-unique)
DFT1.SVs.384T1.MSG.non.unique <- DFT1.SVs.384T1.MSG[which(apply(DFT1.SVs.384T1.MSG[,c(41:88,90:105)], 1, function(x){out <- as.numeric(str_split_fixed(x, '/', 2)[,1]); any(out >= 1)}) == T),]
DFT1.SVs.384T1.SvABA.non.unique <- DFT1.SVs.384T1.SvABA[which(apply(DFT1.SVs.384T1.SvABA[,c(41:88,90:105)], 1, function(x){out <- as.numeric(x); any(out >= 3)}) == T),]
DFT1.SVs.384T1.non.unique <- rbind(DFT1.SVs.384T1.MSG.non.unique, DFT1.SVs.384T1.SvABA.non.unique)
DFT1.SVs.384T1.non.unique <- DFT1.SVs.384T1.non.unique[order(as.character(DFT1.SVs.384T1.non.unique[,"CHROM1"]), as.numeric(DFT1.SVs.384T1.non.unique[,"ORIGINAL POS1"])),]
DFT1.SVs.384T1.non.unique.left.GR <- GRanges(seqnames = as.character(DFT1.SVs.384T1.non.unique[,'CHROM1']), 
                                             ranges = IRanges::IRanges(start = as.numeric(as.character(DFT1.SVs.384T1.non.unique[,'ORIGINAL POS1'])),
                                                                       end = as.numeric(as.character(DFT1.SVs.384T1.non.unique[,'ORIGINAL POS1']))))
DFT1.SVs.384T1.non.unique.right.GR <- GRanges(seqnames = as.character(DFT1.SVs.384T1.non.unique[,'CHROM2']), 
                                              ranges = IRanges::IRanges(start = as.numeric(as.character(DFT1.SVs.384T1.non.unique[,'ORIGINAL POS2'])),
                                                                        end = as.numeric(as.character(DFT1.SVs.384T1.non.unique[,'ORIGINAL POS2']))))
OL.left <- findOverlaps(DFT1.SVs.384T1.non.unique.left.GR, GRanges(seqnames = 1, 
                                                                   ranges = IRanges::IRanges(start = 410000000,
                                                                                             end = 470000000)))
OL.right <- findOverlaps(DFT1.SVs.384T1.non.unique.right.GR, GRanges(seqnames = 1, 
                                                                     ranges = IRanges::IRanges(start = 410000000,
                                                                                               end = 470000000)))
OL.left <- as.matrix(OL.left); OL.right <- as.matrix(OL.right)
DFT1.SVs.384T1.non.unique.local <- DFT1.SVs.384T1.non.unique[OL.left[OL.left[,1] %in% OL.right[,1],1],]
DFT1.SVs.384T1.non.unique.inter.right <- DFT1.SVs.384T1.non.unique[OL.left[!OL.left[,1] %in% OL.right[,1],1],]
DFT1.SVs.384T1.non.unique.inter.left <- DFT1.SVs.384T1.non.unique[OL.right[!OL.right[,1] %in% OL.left[,1],1],]

### 384T1 SVs (all unique)
DFT1.SVs.384T1.MSG.unique <- DFT1.SVs.384T1.MSG[which(apply(DFT1.SVs.384T1.MSG[,c(41:88,90:105)], 1, function(x){out <- as.numeric(str_split_fixed(x, '/', 2)[,1]); all(out < 1)}) == T),]
DFT1.SVs.384T1.SvABA.unique <- DFT1.SVs.384T1.SvABA[which(apply(DFT1.SVs.384T1.SvABA[,c(41:88,90:105)], 1, function(x){out <- as.numeric(x); all(out < 3)}) == T),]
DFT1.SVs.384T1.unique <- rbind(DFT1.SVs.384T1.MSG.unique, DFT1.SVs.384T1.SvABA.unique)
DFT1.SVs.384T1.unique <- DFT1.SVs.384T1.unique[order(as.character(DFT1.SVs.384T1.unique[,"CHROM1"]), as.numeric(DFT1.SVs.384T1.unique[,"ORIGINAL POS1"])),]
DFT1.SVs.384T1.unique.left.GR <- GRanges(seqnames = as.character(DFT1.SVs.384T1.unique[,'CHROM1']), 
                                         ranges = IRanges::IRanges(start = as.numeric(as.character(DFT1.SVs.384T1.unique[,'ORIGINAL POS1'])),
                                                                   end = as.numeric(as.character(DFT1.SVs.384T1.unique[,'ORIGINAL POS1']))))
DFT1.SVs.384T1.unique.right.GR <- GRanges(seqnames = as.character(DFT1.SVs.384T1.unique[,'CHROM2']), 
                                          ranges = IRanges::IRanges(start = as.numeric(as.character(DFT1.SVs.384T1.unique[,'ORIGINAL POS2'])),
                                                                    end = as.numeric(as.character(DFT1.SVs.384T1.unique[,'ORIGINAL POS2']))))
OL.left <- findOverlaps(DFT1.SVs.384T1.unique.left.GR, GRanges(seqnames = 1, 
                                                               ranges = IRanges::IRanges(start = 410000000,
                                                                                         end = 470000000)))
OL.right <- findOverlaps(DFT1.SVs.384T1.unique.right.GR, GRanges(seqnames = 1, 
                                                                 ranges = IRanges::IRanges(start = 410000000,
                                                                                           end = 470000000)))
OL.left <- as.matrix(OL.left); OL.right <- as.matrix(OL.right)
DFT1.SVs.384T1.unique.local <- DFT1.SVs.384T1.unique[OL.left[OL.left[,1] %in% OL.right[,1],1],]
DFT1.SVs.384T1.unique.inter.right <- DFT1.SVs.384T1.unique[OL.left[!OL.left[,1] %in% OL.right[,1],1],]
DFT1.SVs.384T1.unique.inter.left <- DFT1.SVs.384T1.unique[OL.right[!OL.right[,1] %in% OL.left[,1],1],]

## chromosomes
chromosome.ranges <- matrix(0, ncol = 2, nrow = 7)
chromosome.ranges[,1] <- 1
chromosome.ranges[,2] <- c(716413629, 662751787, 611347268, 464895054, 288121652, 254895979, 83081154)
rownames(chromosome.ranges) <- c("1", "2", "3", "4", "5", "6", "X")
chromosome.ranges.GR <- GRanges(seqnames = rownames(chromosome.ranges), 
                                ranges = IRanges(start = as.numeric(chromosome.ranges[,1]),
                                                 end = as.numeric(chromosome.ranges[,2])))

## plot DFT1 circos
pdf('Figure5C_384T1_chromothripsis_circos.pdf', 
    width = 20, height = 19)

circos.par(track.height = 0.2, 
           cell.padding = c(0, 0, 0, 0), 
           start.degree = 90, gap.degree = 6,
           track.margin = c(0.1,0.1))
circos.initialize(factors = rownames(chromosome.ranges), 
                  xlim = chromosome.ranges)

## add outer ring
circos.trackPlotRegion(ylim = c(0, 1), 
                       panel.fun = function(x, y) {get.cell.meta.data("xlim")}, 
                       track.height = 0.2, bg.col = 'black', bg.border = NA,
                       track.index = 1, bg.lwd = 2)

## add circos labels
for (i in 1:length(rownames(chromosome.ranges))){
  circos.axis(h = 'top', 
              sector.index = rownames(chromosome.ranges)[i], 
              major.at = chromosome.ranges[i,2]/2,
              labels = rownames(chromosome.ranges)[i],
              direction = "outside", 
              labels.cex = 10,
              lwd = 2,
              labels.niceFacing = F,
              major.tick = F,
              major.tick.length = 0.2)
}

## add SVs (all in 384T1, minus locus-specific)
DFT1.SVs.384T1.plot <- as.data.frame(DFT1.SVs.384T1.non.unique[,c(4,5,18,19)])
DFT1.SVs.384T1.plot[,2] <- as.numeric(as.character(DFT1.SVs.384T1.plot[,2]))
DFT1.SVs.384T1.plot[,4] <- as.numeric(as.character(DFT1.SVs.384T1.plot[,4]))
circos.genomicLink(region1 = DFT1.SVs.384T1.plot[,c(1,2,2)],
                   region2 = DFT1.SVs.384T1.plot[,c(3,4,4)],
                   col = alpha("black", 0.05), lwd = 5,
                   rou = 0.69)

## add SVs (384T1 locus-specific)
DFT1.SVs.384T1.chr1.plot <- as.data.frame(DFT1.SVs.384T1.unique.local[,c(4,5,18,19)])
DFT1.SVs.384T1.chr1.plot[,2] <- as.numeric(as.character(DFT1.SVs.384T1.chr1.plot[,2]))
DFT1.SVs.384T1.chr1.plot[,4] <- as.numeric(as.character(DFT1.SVs.384T1.chr1.plot[,4]))
circos.genomicLink(region1 = DFT1.SVs.384T1.chr1.plot[,c(1,2,2)],
                   region2 = DFT1.SVs.384T1.chr1.plot[,c(3,4,4)],
                   col = "cornflowerblue", lwd = 7,
                   rou = 0.69)

dev.off()

## clean up environment
rm(list=ls())


# Chromothripsis in 384T1 - Campbellgram #
##########################################

## load and subset 384T1 SVs
DFT1.SVs <- as.matrix(read_xlsx('/Users/mstammnitz/Desktop/DFT_evolution/doc/manuscripts/The Evolutionary History of Two Transmissible Cancers in Tasmanian Devils/Tables/v6/Table-S5_v6.xlsx', sheet = 2))
colnames(DFT1.SVs) <- as.character(DFT1.SVs[2,])
DFT1.SVs <- DFT1.SVs[-c(1:2),]
DFT1.SVs.384T1.MSG <- DFT1.SVs[grep('MSG', DFT1.SVs[,'CALLER']),]
DFT1.SVs.384T1.MSG <- DFT1.SVs.384T1.MSG[which(apply(DFT1.SVs.384T1.MSG[,grep('384T1', colnames(DFT1.SVs.384T1.MSG)),drop=F], 1, function(x){out <- as.numeric(str_split_fixed(x, '/', 2)[,1]); all(out >= 1)}) == T),]
DFT1.SVs.384T1.SvABA <- DFT1.SVs[-grep('MSG', DFT1.SVs[,'CALLER']),]
DFT1.SVs.384T1.SvABA <- DFT1.SVs.384T1.SvABA[which(apply(DFT1.SVs.384T1.SvABA[,grep('384T1', colnames(DFT1.SVs.384T1.SvABA)),drop=F], 1, function(x){out <- as.numeric(x); all(out >= 3)}) == T),]

### 384T1 SVs (all non-unique)
DFT1.SVs.384T1.MSG.non.unique <- DFT1.SVs.384T1.MSG[which(apply(DFT1.SVs.384T1.MSG[,c(41:88,90:105)], 1, function(x){out <- as.numeric(str_split_fixed(x, '/', 2)[,1]); any(out >= 1)}) == T),]
DFT1.SVs.384T1.SvABA.non.unique <- DFT1.SVs.384T1.SvABA[which(apply(DFT1.SVs.384T1.SvABA[,c(41:88,90:105)], 1, function(x){out <- as.numeric(x); any(out >= 3)}) == T),]
DFT1.SVs.384T1.non.unique <- rbind(DFT1.SVs.384T1.MSG.non.unique, DFT1.SVs.384T1.SvABA.non.unique)
DFT1.SVs.384T1.non.unique <- DFT1.SVs.384T1.non.unique[order(as.character(DFT1.SVs.384T1.non.unique[,"CHROM1"]), as.numeric(DFT1.SVs.384T1.non.unique[,"ORIGINAL POS1"])),]
DFT1.SVs.384T1.non.unique.left.GR <- GRanges(seqnames = as.character(DFT1.SVs.384T1.non.unique[,'CHROM1']), 
                                             ranges = IRanges::IRanges(start = as.numeric(as.character(DFT1.SVs.384T1.non.unique[,'ORIGINAL POS1'])),
                                                                       end = as.numeric(as.character(DFT1.SVs.384T1.non.unique[,'ORIGINAL POS1']))))
DFT1.SVs.384T1.non.unique.right.GR <- GRanges(seqnames = as.character(DFT1.SVs.384T1.non.unique[,'CHROM2']), 
                                              ranges = IRanges::IRanges(start = as.numeric(as.character(DFT1.SVs.384T1.non.unique[,'ORIGINAL POS2'])),
                                                                        end = as.numeric(as.character(DFT1.SVs.384T1.non.unique[,'ORIGINAL POS2']))))
OL.left <- findOverlaps(DFT1.SVs.384T1.non.unique.left.GR, GRanges(seqnames = 1, 
                                                                   ranges = IRanges::IRanges(start = 410000000,
                                                                                             end = 470000000)))
OL.right <- findOverlaps(DFT1.SVs.384T1.non.unique.right.GR, GRanges(seqnames = 1, 
                                                                     ranges = IRanges::IRanges(start = 410000000,
                                                                                               end = 470000000)))
OL.left <- as.matrix(OL.left); OL.right <- as.matrix(OL.right)
DFT1.SVs.384T1.non.unique.local <- DFT1.SVs.384T1.non.unique[OL.left[OL.left[,1] %in% OL.right[,1],1],]
DFT1.SVs.384T1.non.unique.inter.right <- DFT1.SVs.384T1.non.unique[OL.left[!OL.left[,1] %in% OL.right[,1],1],]
DFT1.SVs.384T1.non.unique.inter.left <- DFT1.SVs.384T1.non.unique[OL.right[!OL.right[,1] %in% OL.left[,1],1],]

### 384T1 SVs (all unique)
DFT1.SVs.384T1.MSG.unique <- DFT1.SVs.384T1.MSG[which(apply(DFT1.SVs.384T1.MSG[,c(41:88,90:105)], 1, function(x){out <- as.numeric(str_split_fixed(x, '/', 2)[,1]); all(out < 1)}) == T),]
DFT1.SVs.384T1.SvABA.unique <- DFT1.SVs.384T1.SvABA[which(apply(DFT1.SVs.384T1.SvABA[,c(41:88,90:105)], 1, function(x){out <- as.numeric(x); all(out < 3)}) == T),]
DFT1.SVs.384T1.unique <- rbind(DFT1.SVs.384T1.MSG.unique, DFT1.SVs.384T1.SvABA.unique)
DFT1.SVs.384T1.unique <- DFT1.SVs.384T1.unique[order(as.character(DFT1.SVs.384T1.unique[,"CHROM1"]), as.numeric(DFT1.SVs.384T1.unique[,"ORIGINAL POS1"])),]
DFT1.SVs.384T1.unique.left.GR <- GRanges(seqnames = as.character(DFT1.SVs.384T1.unique[,'CHROM1']), 
                                         ranges = IRanges::IRanges(start = as.numeric(as.character(DFT1.SVs.384T1.unique[,'ORIGINAL POS1'])),
                                                                   end = as.numeric(as.character(DFT1.SVs.384T1.unique[,'ORIGINAL POS1']))))
DFT1.SVs.384T1.unique.right.GR <- GRanges(seqnames = as.character(DFT1.SVs.384T1.unique[,'CHROM2']), 
                                          ranges = IRanges::IRanges(start = as.numeric(as.character(DFT1.SVs.384T1.unique[,'ORIGINAL POS2'])),
                                                                    end = as.numeric(as.character(DFT1.SVs.384T1.unique[,'ORIGINAL POS2']))))
OL.left <- findOverlaps(DFT1.SVs.384T1.unique.left.GR, GRanges(seqnames = 1, 
                                                               ranges = IRanges::IRanges(start = 410000000,
                                                                                         end = 470000000)))
OL.right <- findOverlaps(DFT1.SVs.384T1.unique.right.GR, GRanges(seqnames = 1, 
                                                                 ranges = IRanges::IRanges(start = 410000000,
                                                                                           end = 470000000)))
OL.left <- as.matrix(OL.left); OL.right <- as.matrix(OL.right)
DFT1.SVs.384T1.unique.local <- DFT1.SVs.384T1.unique[OL.left[OL.left[,1] %in% OL.right[,1],1],]
DFT1.SVs.384T1.unique.inter.right <- DFT1.SVs.384T1.unique[OL.left[!OL.left[,1] %in% OL.right[,1],1],]
DFT1.SVs.384T1.unique.inter.left <- DFT1.SVs.384T1.unique[OL.right[!OL.right[,1] %in% OL.left[,1],1],]

## load 384T1 normalised LogR data
copynumber.384T1 <- read.csv("/Users/mstammnitz/Desktop/DFT_evolution/doc/manuscripts/The Evolutionary History of Two Transmissible Cancers in Tasmanian Devils/Tables/v6/Supplementary_data/Copynumber_384T1_chr1.csv")
copynumber.384T1.window <- copynumber.384T1[which(copynumber.384T1[,"excluded_from_segmentation"] == F),]

## display CN estimates for SNPs & SNVs
ylims <- c(0.5, 3.5)
lwds <- 2
transp <- 0.6

## plot
pdf("Figure5C_384T1_chromothripsis_Campbellgram.pdf", 
    width = 23, height = 14)
par(mar = c(15, 12, 2, 2))
plot(x = c(410000000,470000000),
     y = c(1,1),
     pch = 16, 
     col = alpha('black', transp),
     cex = 1,
     ylab = '', 
     xlab = '',
     yaxt = 'n', 
     xaxt = 'n',
     bty = 'n',
     xlim = c(410000000,470000000),
     ylim = ylims)

### local non-unique SVs
fit.func <- function(x){
  y <- fit$coefficients[3]*x^2+c(fit$coefficients[2]*x)+fit$coefficients[1]
  return(y)
}
for (i in 1:nrow(DFT1.SVs.384T1.non.unique.local)){
  
  ## arcs
  bp <- DFT1.SVs.384T1.non.unique.local[i,]
  bp <- bp[grep('ORIGINAL', colnames(DFT1.SVs.384T1.non.unique.local))]
  bp.x <- c(as.numeric(bp[1]),c(as.numeric(bp[1])+as.numeric(bp[2]))/2,as.numeric(bp[2]))
  bp.y <- c(3, 3.5, 3)
  fit <- glm(bp.y ~ poly(bp.x, 2, raw = T))
  x.fit.out <- seq(f = bp.x[1], t = bp.x[3], length.out = 5000)
  y.fit.out <- fit.func(x.fit.out)
  lines(x = x.fit.out,
        y = y.fit.out,
        col = alpha('black', transp),
        pch = 16, 
        lwd = lwds)
  
  ## borders
  lines(x = rep(bp[1],2), y = c(0.5,3), col = alpha('black', transp), lwd = lwds)
  lines(x = rep(bp[2],2), y = c(0.5,3), col = alpha('black', transp), lwd = lwds)
  
}

### distant non-unique SVs
for (i in 1:nrow(DFT1.SVs.384T1.non.unique.inter.right)){

  ## arcs
  bp <- DFT1.SVs.384T1.non.unique.inter.right[i,]
  bp <- bp[c(grep('ORIGINAL POS1', colnames(DFT1.SVs.384T1.non.unique.inter.right)),
             grep('CHROM2', colnames(DFT1.SVs.384T1.non.unique.inter.right)))]
  end <- as.numeric(bp[1]) + 10*c(470000000 - as.numeric(bp[1]))

  if(bp[2] != 'X' & bp[2] != '6' & bp[2] != '5'){

    bp.x <- c(as.numeric(bp[1]), 480000000, 480000000 + c(480000000 - as.numeric(bp[1])))

  }else{

    bp.x <- c(400000000 - c(as.numeric(bp[1]) - 400000000), 400000000, as.numeric(bp[1]))

  }

  fit <- glm(c(3, 3.5, 3) ~ poly(bp.x, 2, raw=TRUE))
  x.fit.out <- seq(f = bp.x[1], t = bp.x[3], length.out = 5000)
  y.fit.out <- fit.func(seq(f = bp.x[1], t = bp.x[3], length.out = 5000))
  lines(x = x.fit.out,
        y = y.fit.out,
        col = alpha('black', transp),
        pch = 16,
        lwd = lwds)

  ## borders
  lines(x = rep(bp[1],2), y = c(0.5,3), col = alpha('black', transp), lwd = lwds)

}

### local unique SVs
fit.func <- function(x){
  y <- fit$coefficients[3]*x^2+c(fit$coefficients[2]*x)+fit$coefficients[1]
  return(y)
}
for (i in 1:nrow(DFT1.SVs.384T1.unique.local)){
  
  ## arcs
  bp <- DFT1.SVs.384T1.unique.local[i,]
  bp <- bp[grep('ORIGINAL', colnames(DFT1.SVs.384T1.unique.local))]
  bp.x <- c(as.numeric(bp[1]),c(as.numeric(bp[1])+as.numeric(bp[2]))/2,as.numeric(bp[2]))
  bp.y <- c(3, 3.5, 3)
  fit <- glm(bp.y ~ poly(bp.x, 2, raw = T))
  x.fit.out <- seq(f = bp.x[1], t = bp.x[3], length.out = 5000)
  y.fit.out <- fit.func(x.fit.out)
  lines(x = x.fit.out,
        y = y.fit.out,
        col = 'cornflowerblue', 
        pch = 16, 
        lwd = lwds)
  
  ## borders
  lines(x = rep(bp[1],2), y = c(0.5,3), col = 'cornflowerblue', lwd = lwds)
  lines(x = rep(bp[2],2), y = c(0.5,3), col = 'cornflowerblue', lwd = lwds)
}

### actual data points
points(x = as.numeric(copynumber.384T1.window[,'START']),
       y = as.numeric(copynumber.384T1.window[,'total_cn']),
       pch = 16, 
       col = alpha('black', transp),
       cex = 0.6)

### axes
axis(1, at = seq(from = 410000000, to = 470000000, length = 5), 
     labels = paste0(seq(from = 410000000, to = 470000000, length = 5)/1000000),
     cex.axis = 5, 
     padj = 1, 
     lwd = 5)
mtext(side = 1, 
      text = 'Chromosome 1 [MBP]', 
      line = 12, 
      cex = 7)

axis(2, 
     at = c(1, 2, 3), 
     labels = c(1, 2, 3), 
     cex.axis = 5, 
     las = 2, 
     lwd
     = 5, 
     line = -1)
mtext(side = 2, 
      text = 'Copy number', 
      line = 7, 
      cex = 7)

dev.off()

## clean up environment
rm(list=ls())