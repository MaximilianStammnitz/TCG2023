## Figure 5A
## 'The evolution of two transmissible cancers in Tasmanian devils'
## Stammnitz et al., 2022
## maxrupsta@gmail

library(readxl)
library(stringr)
library(GenomicRanges)
library(circlize)

## set input path(s)
setwd('/Tables')


# DFT1 ancestral events - circos plot #
#######################################

## import structural variant table
DFT1.SVs <- as.matrix(read_xlsx('/Table-S5.xlsx', sheet = 2))
colnames(DFT1.SVs) <- as.character(DFT1.SVs[2,])
DFT1.SVs <- DFT1.SVs[-c(1:2),]

## import copy-number variants
DFT1.CNVs <- as.matrix(read_xlsx('Table-S6.xlsx', sheet = 2))
colnames(DFT1.CNVs) <- as.character(DFT1.CNVs[2,])
DFT1.CNVs <- DFT1.CNVs[-c(1:2),]

## subset ancestral events
DFT1.SVs.ancestral <- DFT1.SVs
DFT1.SVs.ancestral.MSG <- DFT1.SVs.ancestral[grep('/', DFT1.SVs.ancestral[,'377T1']),]
DFT1.SVs.ancestral.MSG <- DFT1.SVs.ancestral.MSG[which(apply(DFT1.SVs.ancestral.MSG[,41:105], 1, function(x){out <- as.numeric(str_split_fixed(x, '/', 2)[,1]); all(out >= 1)}) == T),]
DFT1.SVs.ancestral.SvABA <- DFT1.SVs.ancestral[-grep('/', DFT1.SVs.ancestral[,'377T1']),]
DFT1.SVs.ancestral.SvABA <- DFT1.SVs.ancestral.SvABA[which(apply(DFT1.SVs.ancestral.SvABA[,41:105], 1, function(x){out <- as.numeric(x); all(out >= 3)}) == T),]
DFT1.SVs.ancestral <- rbind(DFT1.SVs.ancestral.MSG, DFT1.SVs.ancestral.SvABA)
DFT1.SVs.ancestral <- DFT1.SVs.ancestral[order(as.character(DFT1.SVs.ancestral[,"CHROM1"]), as.numeric(DFT1.SVs.ancestral[,"ORIGINAL POS1"])),]

DFT1.CNVs.ancestral <- DFT1.CNVs[which(DFT1.CNVs[,'N'] == '78'),]

## reconstruct total copy-number states of MRCAs
chromosome.ranges <- matrix(0, ncol = 2, nrow = 7)
chromosome.ranges[,1] <- 1
chromosome.ranges[,2] <- c(716413629, 662751787, 611347268, 464895054, 288121652, 254895979, 83081154)
rownames(chromosome.ranges) <- c("1", "2", "3", "4", "5", "6", "X")
chromosome.ranges.GR <- GRanges(seqnames = rownames(chromosome.ranges), 
                                ranges = IRanges(start = as.numeric(chromosome.ranges[,1]),
                                                 end = as.numeric(chromosome.ranges[,2])))
DFT1.CNVs.ancestral.total <- cbind(DFT1.CNVs.ancestral[,1:12], 'MRCA' = rep(NA, nrow(DFT1.CNVs.ancestral)))
DFT1.CNVs.ancestral.total.GR <- GRanges(seqnames = as.character(DFT1.CNVs.ancestral.total[,"CHROM"]),
                                        ranges = IRanges(start = as.numeric(DFT1.CNVs.ancestral.total[,"START"]), 
                                                         end = as.numeric(DFT1.CNVs.ancestral.total[,"END"])))
DFT1.CNVs.ancestral.total.disjoin.GR <- disjoin(sort(c(chromosome.ranges.GR,DFT1.CNVs.ancestral.total.GR)))
elementMetadata(DFT1.CNVs.ancestral.total.disjoin.GR)['MRCA'] <- 2
OL <- as.matrix(findOverlaps(DFT1.CNVs.ancestral.total.GR, DFT1.CNVs.ancestral.total.disjoin.GR))
for(j in 1:nrow(OL)){
  
  tmp.event <- as.character(DFT1.CNVs.ancestral.total[as.numeric(OL[j,1]),'EVENT TYPE'])
  
  if(tmp.event == 'GAIN'){
    
    elementMetadata(DFT1.CNVs.ancestral.total.disjoin.GR)[as.numeric(OL[j,2]),'MRCA'] <- elementMetadata(DFT1.CNVs.ancestral.total.disjoin.GR)[as.numeric(OL[j,2]),'MRCA'] + 1
    
  }else if(tmp.event == 'LOSS'){
    
    elementMetadata(DFT1.CNVs.ancestral.total.disjoin.GR)[as.numeric(OL[j,2]),'MRCA'] <- elementMetadata(DFT1.CNVs.ancestral.total.disjoin.GR)[as.numeric(OL[j,2]),'MRCA'] - 1
    
  }
  
}
DFT1.CNVs.ancestral.total <- as.data.frame(DFT1.CNVs.ancestral.total.disjoin.GR)
DFT1.CNVs.ancestral.total <- DFT1.CNVs.ancestral.total[,-5]
colnames(DFT1.CNVs.ancestral.total)[1:4] <- c('CHROM', 'START', 'END', 'WIDTH')

## plot DFT1 circos
pdf('Figure5A_DFT1_ancestral_circos.pdf', 
    width = 30, height = 28)
  
circos.par("track.height" = 0.2, cell.padding = c(0, 0, 0, 0), start.degree = 90, gap.degree = 10,
           track.margin = c(0.06,0.06))
circos.initialize(factors = rownames(chromosome.ranges), 
                  xlim = chromosome.ranges)

## add outer ring
circos.trackPlotRegion(ylim = c(-0.2, 3.2), 
                       panel.fun = function(x, y) {get.cell.meta.data("xlim")}, 
                       track.height = 0.5, bg.col = 'grey90', bg.border = NA,
                       track.index = 1, bg.lwd = 2)

## add circos labels
for (i in 1:length(rownames(chromosome.ranges))){
  circos.axis(h = 'top', sector.index = rownames(chromosome.ranges)[i], 
              major.at = chromosome.ranges[i,2]/2,
              labels = rownames(chromosome.ranges)[i],
              direction = "outside", labels.cex = 10,
              major.tick = F,
              major.tick.length = 0.25,
              lwd = 2,
              labels.niceFacing = F)
}

## add Y-axes
for(i in c(1:6,'X')){

  circos.yaxis(side = "left",
               at = c(0,1,2,3),
               labels = c(0,1,2,3),
               tick = T,
               sector.index = i,
               labels.cex = 8,
               labels.niceFacing = TRUE,
               lwd = 4,
               col = 'black',
               labels.col = 'black')

}

## add ancestral total copy-number lines (horizontal)
for (i in 1:nrow(DFT1.CNVs.ancestral.total)){
  circos.trackLines(factors = rep(DFT1.CNVs.ancestral.total$CHROM[i],2),
                    x = c(DFT1.CNVs.ancestral.total$START[i],DFT1.CNVs.ancestral.total$END[i]),
                    y = c(DFT1.CNVs.ancestral.total$MRCA[i],DFT1.CNVs.ancestral.total$MRCA[i]),
                    col = 'cornflowerblue', lwd = 6, type = 'l')
}

## add ancestral total copy-number lines (vertical)
for (i in 1:length(unique(DFT1.CNVs.ancestral.total$CHROM))){

  tmp.chr <- unique(DFT1.CNVs.ancestral.total$CHROM)[i]
  tmp.chr.id <- which(DFT1.CNVs.ancestral.total$CHROM == tmp.chr)

  for(j in tmp.chr.id[-which.max(tmp.chr.id)]){

    circos.trackLines(factors = rep(tmp.chr,2),
                      x = c(DFT1.CNVs.ancestral.total$END[j],DFT1.CNVs.ancestral.total$END[j]),
                      y = c(DFT1.CNVs.ancestral.total$MRCA[j],DFT1.CNVs.ancestral.total$MRCA[j+1]),
                      col = 'cornflowerblue', lwd = 6, type = 'l')

  }

}

## add SVs
DFT1.SVs.ancestral.plot <- as.data.frame(DFT1.SVs.ancestral[,c(4,5,18,19)])
DFT1.SVs.ancestral.plot[,2] <- as.numeric(as.character(DFT1.SVs.ancestral.plot[,2]))
DFT1.SVs.ancestral.plot[,4] <- as.numeric(as.character(DFT1.SVs.ancestral.plot[,4]))
circos.genomicLink(region1 = DFT1.SVs.ancestral.plot[,c(1,2,2)],
                   region2 = DFT1.SVs.ancestral.plot[,c(3,4,4)],
                   col = "cornflowerblue", lwd = 6,
                   rou = 0.4)

dev.off()

## clean up environment
rm(list=ls())


# DFT2 ancestral events - circos plot #
#######################################

## import structural variant table
DFT2.SVs <- as.matrix(read_xlsx('Table-S5.xlsx', sheet = 3))
colnames(DFT2.SVs) <- as.character(DFT2.SVs[2,])
DFT2.SVs <- DFT2.SVs[-c(1:2),]

## import copy-number variants
DFT2.CNVs <- as.matrix(read_xlsx('Table-S6.xlsx', sheet = 3))
colnames(DFT2.CNVs) <- as.character(DFT2.CNVs[2,])
DFT2.CNVs <- DFT2.CNVs[-c(1:2),]

## subset ancestral events
DFT2.SVs.ancestral <- DFT2.SVs
DFT2.SVs.ancestral.MSG <- DFT2.SVs.ancestral[grep('/', DFT2.SVs.ancestral[,'202T1']),]
DFT2.SVs.ancestral.MSG <- DFT2.SVs.ancestral.MSG[which(apply(DFT2.SVs.ancestral.MSG[,41:81], 1, function(x){out <- as.numeric(str_split_fixed(x, '/', 2)[,1]); all(out >= 1)}) == T),]
DFT2.SVs.ancestral.SvABA <- DFT2.SVs.ancestral[-grep('/', DFT2.SVs.ancestral[,'202T1']),]
DFT2.SVs.ancestral.SvABA <- DFT2.SVs.ancestral.SvABA[which(apply(DFT2.SVs.ancestral.SvABA[,41:81], 1, function(x){out <- as.numeric(x); all(out >= 3)}) == T),]
DFT2.SVs.ancestral <- rbind(DFT2.SVs.ancestral.MSG, DFT2.SVs.ancestral.SvABA)
DFT2.SVs.ancestral <- DFT2.SVs.ancestral[order(as.character(DFT2.SVs.ancestral[,"CHROM1"]), as.numeric(DFT2.SVs.ancestral[,"ORIGINAL POS1"])),]

DFT2.CNVs.ancestral <- DFT2.CNVs[which(DFT2.CNVs[,'N'] == '41'),]

## reconstruct total copy-number states of MRCAs
chromosome.ranges <- matrix(0, ncol = 2, nrow = 7)
chromosome.ranges[,1] <- 1
chromosome.ranges[,2] <- c(716413629, 662751787, 611347268, 464895054, 288121652, 254895979, 83081154)
rownames(chromosome.ranges) <- c("1", "2", "3", "4", "5", "6", "X")
chromosome.ranges.GR <- GRanges(seqnames = rownames(chromosome.ranges), 
                                ranges = IRanges(start = as.numeric(chromosome.ranges[,1]),
                                                 end = as.numeric(chromosome.ranges[,2])))

DFT2.CNVs.ancestral.total <- cbind(DFT2.CNVs.ancestral[,1:10], 'MRCA' = rep(NA, nrow(DFT2.CNVs.ancestral)))
DFT2.CNVs.ancestral.total.GR <- GRanges(seqnames = as.character(DFT2.CNVs.ancestral.total[,"CHROM"]),
                                        ranges = IRanges(start = as.numeric(DFT2.CNVs.ancestral.total[,"START"]), 
                                                         end = as.numeric(DFT2.CNVs.ancestral.total[,"END"])))
DFT2.CNVs.ancestral.total.disjoin.GR <- disjoin(sort(c(chromosome.ranges.GR,DFT2.CNVs.ancestral.total.GR)))
elementMetadata(DFT2.CNVs.ancestral.total.disjoin.GR)['MRCA'] <- c(rep(2,91),1)
OL <- as.matrix(findOverlaps(DFT2.CNVs.ancestral.total.GR, DFT2.CNVs.ancestral.total.disjoin.GR))
for(j in 1:nrow(OL)){
  
  tmp.event <- as.character(DFT2.CNVs.ancestral.total[as.numeric(OL[j,1]),'EVENT TYPE'])
  
  if(tmp.event == 'GAIN'){
    
    elementMetadata(DFT2.CNVs.ancestral.total.disjoin.GR)[as.numeric(OL[j,2]),'MRCA'] <- elementMetadata(DFT2.CNVs.ancestral.total.disjoin.GR)[as.numeric(OL[j,2]),'MRCA'] + 1
    
  }else if(tmp.event == 'LOSS'){
    
    elementMetadata(DFT2.CNVs.ancestral.total.disjoin.GR)[as.numeric(OL[j,2]),'MRCA'] <- elementMetadata(DFT2.CNVs.ancestral.total.disjoin.GR)[as.numeric(OL[j,2]),'MRCA'] - 1
    
  }
  
}
DFT2.CNVs.ancestral.total <- as.data.frame(DFT2.CNVs.ancestral.total.disjoin.GR)
DFT2.CNVs.ancestral.total <- DFT2.CNVs.ancestral.total[,-5]
colnames(DFT2.CNVs.ancestral.total)[1:4] <- c('CHROM', 'START', 'END', 'WIDTH')

## plot DFT2 circos
pdf('Figure5A_DFT2_ancestral_circos.pdf', 
    width = 30, height = 28)

circos.par("track.height" = 0.2, cell.padding = c(0, 0, 0, 0), start.degree = 90, gap.degree = 10,
           track.margin = c(0.06,0.06))
circos.initialize(factors = rownames(chromosome.ranges), 
                  xlim = chromosome.ranges)

## add outer ring
circos.trackPlotRegion(ylim = c(-0.2, 3.2), 
                       panel.fun = function(x, y) {get.cell.meta.data("xlim")}, 
                       track.height = 0.5, bg.col = 'grey90', bg.border = NA,
                       track.index = 1, bg.lwd = 2)

## add circos labels
for (i in 1:length(rownames(chromosome.ranges))){
  circos.axis(h = 'top', sector.index = rownames(chromosome.ranges)[i], 
              major.at = chromosome.ranges[i,2]/2,
              labels = rownames(chromosome.ranges)[i],
              direction = "outside", labels.cex = 10,
              major.tick = F,
              major.tick.length = 0.25,
              lwd = 2,
              labels.niceFacing = F)
}

## add Y-axes
for(i in c(1:6,'X')){
  
  circos.yaxis(side = "left",
               at = c(0,1,2,3),
               labels = c(0,1,2,3),
               tick = T,
               sector.index = i,
               labels.cex = 8,
               labels.niceFacing = TRUE,
               lwd = 4,
               col = 'black',
               labels.col = 'black')
  
}

## add ancestral total copy-number lines (horizontal)
for (i in 1:nrow(DFT2.CNVs.ancestral.total)){
  circos.trackLines(factors = rep(DFT2.CNVs.ancestral.total$CHROM[i],2),
                    x = c(DFT2.CNVs.ancestral.total$START[i],DFT2.CNVs.ancestral.total$END[i]),
                    y = c(DFT2.CNVs.ancestral.total$MRCA[i],DFT2.CNVs.ancestral.total$MRCA[i]),
                    col = 'red', lwd = 6, type = 'l')
}

## add ancestral total copy-number lines (vertical)
for (i in 1:length(unique(DFT2.CNVs.ancestral.total$CHROM))){
  
  tmp.chr <- unique(DFT2.CNVs.ancestral.total$CHROM)[i]
  tmp.chr.id <- which(DFT2.CNVs.ancestral.total$CHROM == tmp.chr)
  
  for(j in tmp.chr.id[-which.max(tmp.chr.id)]){
    
    circos.trackLines(factors = rep(tmp.chr,2),
                      x = c(DFT2.CNVs.ancestral.total$END[j],DFT2.CNVs.ancestral.total$END[j]),
                      y = c(DFT2.CNVs.ancestral.total$MRCA[j],DFT2.CNVs.ancestral.total$MRCA[j+1]),
                      col = 'red', lwd = 6, type = 'l')
    
  }
  
}

## add SVs
DFT2.SVs.ancestral.plot <- as.data.frame(DFT2.SVs.ancestral[,c(4,5,18,19)])
DFT2.SVs.ancestral.plot[,2] <- as.numeric(as.character(DFT2.SVs.ancestral.plot[,2]))
DFT2.SVs.ancestral.plot[,4] <- as.numeric(as.character(DFT2.SVs.ancestral.plot[,4]))
circos.genomicLink(region1 = DFT2.SVs.ancestral.plot[,c(1,2,2)],
                   region2 = DFT2.SVs.ancestral.plot[,c(3,4,4)],
                   col = "red", lwd = 6,
                   rou = 0.4)

dev.off()

## clean up environment
rm(list=ls())


# Anal-sac carcinoma (340T) circos plot #
#########################################

## import structural variant table
T340.SVs <- as.matrix(read_xlsx('Table-S5.xlsx', sheet = 4))
colnames(T340.SVs) <- as.character(T340.SVs[2,])
T340.SVs <- T340.SVs[-c(1:2),]

## import copy-number variants
T340.CNVs <- as.matrix(read_xlsx('Table-S6.xlsx', sheet = 4))
colnames(T340.CNVs) <- as.character(T340.CNVs[2,])
T340.CNVs <- T340.CNVs[-c(1:2),]

## subset ancestral events
T340.SVs.ancestral <- T340.SVs

T340.CNVs.ancestral <- T340.CNVs

## reconstruct total copy-number states of MRCAs
chromosome.ranges <- matrix(0, ncol = 2, nrow = 7)
chromosome.ranges[,1] <- 1
chromosome.ranges[,2] <- c(716413629, 662751787, 611347268, 464895054, 288121652, 254895979, 83081154)
rownames(chromosome.ranges) <- c("1", "2", "3", "4", "5", "6", "X")
chromosome.ranges.GR <- GRanges(seqnames = rownames(chromosome.ranges), 
                                ranges = IRanges(start = as.numeric(chromosome.ranges[,1]),
                                                 end = as.numeric(chromosome.ranges[,2])))
T340.CNVs.ancestral.total <- T340.CNVs.ancestral[,c(2,3,4,5,9)]
T340.CNVs.ancestral.total[,5] <- as.numeric(T340.CNVs.ancestral.total[,5])
colnames(T340.CNVs.ancestral.total) <- c('CHROM', 'START', 'END', 'WIDTH', 'MRCA')
T340.CNVs.ancestral.total <- T340.CNVs.ancestral.total[-405,] ## remove Y-loss
T340.CNVs.ancestral.total <- as.data.frame(T340.CNVs.ancestral.total)
T340.CNVs.ancestral.total$START <- as.numeric(as.character(T340.CNVs.ancestral.total$START))
T340.CNVs.ancestral.total$END <- as.numeric(as.character(T340.CNVs.ancestral.total$END))
T340.CNVs.ancestral.total$WIDTH <- as.numeric(as.character(T340.CNVs.ancestral.total$WIDTH))
T340.CNVs.ancestral.total$MRCA <- as.numeric(as.character(T340.CNVs.ancestral.total$MRCA))
T340.CNVs.ancestral.total$MRCA[T340.CNVs.ancestral.total$MRCA > 50] <- 50 ## cap at CN >= 50

## plot 340T circos
pdf('Figure5A_340T_circos.pdf', 
    width = 30, height = 28)

circos.par("track.height" = 0.2, cell.padding = c(0, 0, 0, 0), start.degree = 90, gap.degree = 10,
           track.margin = c(0.06,0.06))
circos.initialize(factors = rownames(chromosome.ranges), 
                  xlim = chromosome.ranges)

## add outer ring
circos.trackPlotRegion(ylim = c(-2, 32), 
                       panel.fun = function(x, y) {get.cell.meta.data("xlim")}, 
                       track.height = 0.5, bg.col = 'grey90', bg.border = NA,
                       track.index = 1, bg.lwd = 2)

## add circos labels
for (i in 1:length(rownames(chromosome.ranges))){
  circos.axis(h = 'top', sector.index = rownames(chromosome.ranges)[i], 
              major.at = chromosome.ranges[i,2]/2,
              labels = rownames(chromosome.ranges)[i],
              direction = "outside", labels.cex = 10,
              major.tick.length = 0.25*10,
              major.tick = F,
              lwd = 2,
              labels.niceFacing = F)
}

## add Y-axes
for(i in c(1:6,'X')){
  
  circos.yaxis(side = "left",
               at = c(0,10,20,30),
               labels = c(0,10,20,30),
               tick = T,
               sector.index = i,
               labels.cex = 8,
               labels.niceFacing = TRUE,
               lwd = 4,
               col = 'black',
               labels.col = 'black')
  
}

## add ancestral total copy-number lines (horizontal)
for (i in 1:nrow(T340.CNVs.ancestral.total)){
  circos.trackLines(factors = rep(T340.CNVs.ancestral.total$CHROM[i],2),
                    x = c(T340.CNVs.ancestral.total$START[i],T340.CNVs.ancestral.total$END[i]),
                    y = c(T340.CNVs.ancestral.total$MRCA[i],T340.CNVs.ancestral.total$MRCA[i]),
                    col = 'darkorange', lwd = 6, type = 'l')
}

## add ancestral total copy-number lines (vertical)
for (i in 1:length(unique(T340.CNVs.ancestral.total$CHROM))){
  
  tmp.chr <- unique(T340.CNVs.ancestral.total$CHROM)[i]
  tmp.chr.id <- which(T340.CNVs.ancestral.total$CHROM == tmp.chr)
  
  for(j in tmp.chr.id[-which.max(tmp.chr.id)]){
    
    circos.trackLines(factors = rep(tmp.chr,2),
                      x = c(T340.CNVs.ancestral.total$END[j],T340.CNVs.ancestral.total$END[j]),
                      y = c(T340.CNVs.ancestral.total$MRCA[j],T340.CNVs.ancestral.total$MRCA[j+1]),
                      col = 'darkorange', lwd = 6, type = 'l')
    
  }
  
}

## add SVs
T340.SVs.ancestral.plot <- as.data.frame(T340.SVs.ancestral[,c(4,5,11,12)])
T340.SVs.ancestral.plot[,2] <- as.numeric(as.character(T340.SVs.ancestral.plot[,2]))
T340.SVs.ancestral.plot[,4] <- as.numeric(as.character(T340.SVs.ancestral.plot[,4]))
circos.genomicLink(region1 = T340.SVs.ancestral.plot[,c(1,2,2)],
                   region2 = T340.SVs.ancestral.plot[,c(3,4,4)],
                   col = "darkorange", lwd = 6,
                   rou = 0.4)

dev.off()

## clean up environment
rm(list=ls())
