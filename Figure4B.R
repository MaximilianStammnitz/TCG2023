## Figure 4B
## 'The evolution of two transmissible cancers in Tasmanian devils'
## Stammnitz et al., 2022
## maxrupsta@gmail

library(readxl)
library(stringr)
library(GenomicRanges)
library(circlize)

## set input path(s)
setwd('/Users/ms37/Desktop/Labwork/DFT_evolution/doc/manuscripts/The Evolutionary History of Two Transmissible Cancers in Tasmanian Devils/Tables/v6/')


# Chr1:517 MB LINE-1 source element - circos plot # 
###################################################

## import LINE-1 list
DFT2.LINE1 <- as.matrix(read_xlsx('/Users/mstammnitz/Desktop/DFT_evolution/doc/manuscripts/The Evolutionary History of Two Transmissible Cancers in Tasmanian Devils/Tables/v6/Table-S4_v6.xlsx', sheet = 3))
colnames(DFT2.LINE1) <- DFT2.LINE1[2,]
DFT2.LINE1 <- DFT2.LINE1[-c(1:2),]

## create binary presence/absence table
DFT2.LINE1.MSG <- DFT2.LINE1[grep('/', DFT2.LINE1[,"202T1"]),]
DFT2.LINE1.MSG <- cbind(DFT2.LINE1.MSG[,1:27],
                        apply(DFT2.LINE1.MSG[,28:68], 2, function(x){out <- as.numeric(str_split_fixed(x, '/', 2)[,1]); out[which(out >= 1)] <- 1; return(out)}))
DFT2.LINE1.SvABA <- DFT2.LINE1[-grep('/', DFT2.LINE1[,"202T1"]),]
DFT2.LINE1.SvABA <- cbind(DFT2.LINE1.SvABA[,1:27],
                          apply(DFT2.LINE1.SvABA[,28:68], 2, function(x){out <- as.numeric(x); out[which(out < 3)] <- 0; out[which(out >= 3)] <- 1; return(out)}))
DFT2.LINE1 <- rbind(DFT2.LINE1.MSG, DFT2.LINE1.SvABA)
DFT2.LINE1 <- DFT2.LINE1[order(as.character(DFT2.LINE1[,"CHROM1"]), as.numeric(DFT2.LINE1[,"ORIGINAL POS1"])),]
for(i in 27:68){
  DFT2.LINE1[which(is.na(DFT2.LINE1[,i]) == T),i] <- 0
}
DFT2.LINE1.transductions <- DFT2.LINE1[which(DFT2.LINE1[,"LINE-1 INSERTION TYPE"] == "3' transduction"),]

## mSarHar1.11 chromosomes
chromosome.ranges <- matrix(0, ncol = 2, nrow = 7)
chromosome.ranges[,1] <- 1
chromosome.ranges[,2] <- c(716413629, 662751787, 611347268, 464895054, 288121652, 254895979, 83081154)
rownames(chromosome.ranges) <- c("1", "2", "3", "4", "5", "6", "X")
chromosome.ranges.GR <- GRanges(seqnames = rownames(chromosome.ranges), 
                                ranges = IRanges(start = as.numeric(chromosome.ranges[,1]),
                                                 end = as.numeric(chromosome.ranges[,2])))

## subset events involving the highly active L1 element (N = 29)
DFT2.active.L1.source <- DFT2.LINE1.transductions[which(DFT2.LINE1.transductions[,"CHROM1"] == 1 & abs(as.numeric(DFT2.LINE1.transductions[,"ORIGINAL POS1"]) - 516590773) < 10000 |
                                                        DFT2.LINE1.transductions[,"CHROM2"] == 1 & abs(as.numeric(DFT2.LINE1.transductions[,"ORIGINAL POS2"]) - 516590773) < 10000),]


## plot DFT2 circos
pdf('Figure4B_circos_DFT2_source_element.pdf', 
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

## foreground (combination-unique)
DFT2.active.L1.source.plot <- as.data.frame(DFT2.active.L1.source[,c(5,6,12,13)])
DFT2.active.L1.source.plot[,2] <- as.numeric(as.character(DFT2.active.L1.source.plot[,2]))
DFT2.active.L1.source.plot[,4] <- as.numeric(as.character(DFT2.active.L1.source.plot[,4]))
circos.genomicLink(region1 = DFT2.active.L1.source.plot[,c(1,2,2)],
                   region2 = DFT2.active.L1.source.plot[,c(3,4,4)],
                   col = "red", lwd = 10,
                   rou = 0.69)

dev.off()

## clean up environment
rm(list=ls())