## Figure 5E
## 'The evolution of two transmissible cancers in Tasmanian devils'
## Stammnitz et al., 2022
## maxrupsta@gmail

library(readxl)
library(treeio)
library(ggplot2)
library(ggtree)
library(stringr)
library(GenomicRanges)
library(circlize)

## set input path(s)
setwd('/Users/ms37/Desktop/Labwork/DFT_evolution/doc/manuscripts/The Evolutionary History of Two Transmissible Cancers in Tasmanian Devils/Tables/v6/')


# Examples of (late) DFT2 chromoplexy - tree examples #
#######################################################

## number of substitutions
DFT2.truncal <- 1335 ### number of genome-wide substitutions in the DFT2 trunk
DFT2.truncal.median <- DFT2.truncal - 818 ### median expected number of somatic substitutions in the DFT2 trunk

## import and process BEAST MCMC draws on DFT1 mutation rate and date of origin (generated with Tracer)
masked.mSarHar11.1 <- c('A' = 951740967, 'C' = 540077344, 'G' = 539833602, 'T' = 952098282)
DFT2.beast.trees.MRCAs.rates <- cbind(read.table('/Users/mstammnitz/Desktop/DFT_evolution/doc/manuscripts/The Evolutionary History of Two Transmissible Cancers in Tasmanian Devils/Tables/v6/Supplementary_data/BEAST_DFT2_MRCAs.txt', header = T),
                                      'rates' = read.table('/Users/mstammnitz/Desktop/DFT_evolution/doc/manuscripts/The Evolutionary History of Two Transmissible Cancers in Tasmanian Devils/Tables/v6/Supplementary_data/BEAST_DFT2_rates.txt', header = T)[,2])
DFT2.beast.trees.MRCAs.rates[,'rates'] <- DFT2.beast.trees.MRCAs.rates[,'rates']*
  sum(masked.mSarHar11.1) ### convert to genome-wide mutations per year
DFT2.beast.trees.MRCAs.rates <- cbind(DFT2.beast.trees.MRCAs.rates,
                                      'age at 0 singletons' = rep(NA, nrow(DFT2.beast.trees.MRCAs.rates)),
                                      'age at 818 singletons' = rep(NA, nrow(DFT2.beast.trees.MRCAs.rates)))
DFT2.beast.trees.MRCAs.rates[,'age at 0 singletons'] <- DFT2.beast.trees.MRCAs.rates[,'age.root.'] - 
  c(DFT2.truncal/DFT2.beast.trees.MRCAs.rates[,'rates'])
DFT2.beast.trees.MRCAs.rates[,'age at 818 singletons'] <- DFT2.beast.trees.MRCAs.rates[,'age.root.'] - 
  c(DFT2.truncal.median/DFT2.beast.trees.MRCAs.rates[,'rates'])

## import and process DFT2 maximum clade consensus tree
DFT2.beast.tree.hpd <- read.beast('/Users/mstammnitz/Desktop/DFT_evolution/doc/manuscripts/The Evolutionary History of Two Transmissible Cancers in Tasmanian Devils/Tables/v6/Supplementary_data/BEAST_DFT2.mcc')
DFT2.beast.tree.hpd@data[which.max(unlist(DFT2.beast.tree.hpd@data[,'height'])),"height_0.95_HPD"][[1]][[1]][2] <- 2018.66575342466 - 
  as.numeric(quantile(DFT2.beast.trees.MRCAs.rates[,'age at 0 singletons'], probs = c(0.05)))

## mark DFT2 clades
clade_a <- MRCA(DFT2.beast.tree.hpd, c('339T', '638T1'))
clade_b <- MRCA(DFT2.beast.tree.hpd, c('1529T2', '1545T2'))

## plot
pdf("Figure5E_DFT2_chromoplexy_tree_examples.pdf",
    width = 4, height = 8)

## hack in the tip edge colors
tip.edges <- grep('^202T1$', DFT2.beast.tree.hpd@phylo$tip.label)

## hack in the internal edge colors
internal.edges1 <- MRCA(DFT2.beast.tree.hpd, c('1524T1', '1548T1'))
internal.edges2 <- MRCA(DFT2.beast.tree.hpd, c('338T', '1515T1'))
internal.edges3 <- MRCA(DFT2.beast.tree.hpd, c('1545T2', '1529T4'))

## plot actual tree
ggtree(DFT2.beast.tree.hpd, 
       root.position = 2018.66575342466 - max(DFT2.beast.tree.hpd@data[,'height']), 
       layout = 'rectangular', 
       lwd = 1, 
       ladderize = T) +
  scale_x_continuous(breaks = seq(f=2008, t=2020, by=3),
                     labels = seq(f=2008, t=2020, by=3),
                     limits = c(2008,2020)) +
  geom_point2(aes(subset=(node %in% c(internal.edges1,internal.edges2,internal.edges3,grep('^202T1$', DFT2.beast.tree.hpd@phylo$tip.label)))), 
              shape = 16, size = 7, 
              fill = 'red',
              col = 'red') +
  geom_text2(aes(subset = node %in% internal.edges1), 
             label = 'I', nudge_x = 0.8, size = 12, col = 'red') +
  geom_text2(aes(subset = node %in% internal.edges2), 
             label = 'II', nudge_x = 1, size = 12, col = 'red') +
  geom_text2(aes(subset = node %in% grep('^202T1$', DFT2.beast.tree.hpd@phylo$tip.label)), 
             label = 'III',nudge_x = 1.15, size = 12, col = 'red') +
  geom_text2(aes(subset = node %in% internal.edges3), 
             label = 'IV', nudge_x = -1.1, nudge_y = 1.5, size = 12, col = 'red') +
  geom_rootedge(rootedge = 2018.66575342466 - max(DFT2.beast.tree.hpd@data[,'height']) - mean(DFT2.beast.trees.MRCAs.rates[,'age at 818 singletons']),
                lwd = 1, col = 'black') +
  geom_point2(aes(subset = (node %in% 42),
                  x = x - c(2018.66575342466 - max(DFT2.beast.tree.hpd@data[,'height']) - mean(DFT2.beast.trees.MRCAs.rates[,'age at 818 singletons']))),
              shape = 16, size = 3, col = 'black') +
  theme_classic(base_size = 20) +
  theme(axis.text=element_blank(), 
        legend.position = "none",
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), 
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        line = element_blank())
dev.off()

## clean up environment
rm(list=ls())


# Examples of (late) DFT2 chromoplexy - circos plots #
######################################################

## load SVs
DFT2.SVs <- as.matrix(read_xlsx('/Users/mstammnitz/Desktop/DFT_evolution/doc/manuscripts/The Evolutionary History of Two Transmissible Cancers in Tasmanian Devils/Tables/v6/Table-S5_v6.xlsx', sheet = 3))
colnames(DFT2.SVs) <- as.character(DFT2.SVs[2,])
DFT2.SVs <- DFT2.SVs[-c(1:2),]
DFT2.SVs.summary <- DFT2.SVs
DFT2.SVs.MSG <- DFT2.SVs.summary[grep('MSG', DFT2.SVs.summary[,"CALLER"]),]
DFT2.SVs.MSG <- cbind(DFT2.SVs.MSG[,1:40], apply(DFT2.SVs.MSG[,41:81], 2, function(x){out <- as.numeric(str_split_fixed(x, '/', 2)[,1]); out[which(out >= 1)] <- 1; return(out)}))
DFT2.SVs.SvABA <- DFT2.SVs.summary[-grep('MSG', DFT2.SVs.summary[,"CALLER"]),]
DFT2.SVs.SvABA <- cbind(DFT2.SVs.SvABA[,1:40], apply(DFT2.SVs.SvABA[,41:81], 2, function(x){out <- as.numeric(x); out[which(out < 3)] <- 0; out[which(out >= 3)] <- 1; return(out)}))
DFT2.SVs.summary <- rbind(DFT2.SVs.MSG, DFT2.SVs.SvABA)
DFT2.SVs.summary <- DFT2.SVs.summary[order(as.character(DFT2.SVs.summary[,"CHROM1"]), as.numeric(DFT2.SVs.summary[,"ORIGINAL POS1"])),]

## chromosomes
chromosome.ranges <- matrix(0, ncol = 2, nrow = 7)
chromosome.ranges[,1] <- 1
chromosome.ranges[,2] <- c(716413629, 662751787, 611347268, 464895054, 288121652, 254895979, 83081154)
rownames(chromosome.ranges) <- c("1", "2", "3", "4", "5", "6", "X")
chromosome.ranges.GR <- GRanges(seqnames = rownames(chromosome.ranges), 
                                ranges = IRanges(start = as.numeric(chromosome.ranges[,1]),
                                                 end = as.numeric(chromosome.ranges[,2])))

## subset 1524T1 & 1548T1 ones
DFT2.SVs.1524T1.1548T1 <- DFT2.SVs.summary[which(apply(DFT2.SVs.summary[,c(51,52)], 1, function(x){all(x == 1)}) == T),]
DFT2.SVs.1524T1.1548T1.unique <- DFT2.SVs.1524T1.1548T1[which(apply(DFT2.SVs.1524T1.1548T1[,c(41:50,53:81)], 1, function(x){all(x == 0)}) == T),]
DFT2.SVs.1524T1.1548T1.unique <- DFT2.SVs.1524T1.1548T1.unique[c(5,6,7,10,11,12,13,15,16,20,21,22,23,24,25,28,35),]

## plot DFT2 circos
pdf('Figure5E_circos_I_1524T1_1548T1.pdf', 
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

## unique unclustered SVs
DFT2.SVs.1524T1.1548T1.unique.plot <- as.data.frame(DFT2.SVs.1524T1.1548T1.unique[,c(4,5,18,19)])
DFT2.SVs.1524T1.1548T1.unique.plot[,2] <- as.numeric(as.character(DFT2.SVs.1524T1.1548T1.unique.plot[,2]))
DFT2.SVs.1524T1.1548T1.unique.plot[,4] <- as.numeric(as.character(DFT2.SVs.1524T1.1548T1.unique.plot[,4]))
circos.genomicLink(region1 = DFT2.SVs.1524T1.1548T1.unique.plot[,c(1,2,2)],
                   region2 = DFT2.SVs.1524T1.1548T1.unique.plot[,c(3,4,4)],
                   col = "red", lwd = 20,
                   rou = 0.69)

dev.off()

## subset shared 338T & 1515T1 ones
DFT2.SVs.338T.1515T1 <- DFT2.SVs.summary[which(apply(DFT2.SVs.summary[,c(43,44)], 1, function(x){all(x == 1)}) == T),]
DFT2.SVs.338T.1515T1unique <- DFT2.SVs.338T.1515T1[which(apply(DFT2.SVs.338T.1515T1[,c(41:42,45:81)], 1, function(x){all(x == 0)}) == T),]
DFT2.SVs.338T.1515T1unique <- DFT2.SVs.338T.1515T1unique[c(1,2,3,8,9),] ## chromoplexy cluster

## plot DFT2 circos
pdf('Figure5E_circos_II_338T_1515T1.pdf', 
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

## unique unclustered SVs
DFT2.SVs.338T.1515T1unique.plot <- as.data.frame(DFT2.SVs.338T.1515T1unique[,c(4,5,18,19)])
DFT2.SVs.338T.1515T1unique.plot[,2] <- as.numeric(as.character(DFT2.SVs.338T.1515T1unique.plot[,2]))
DFT2.SVs.338T.1515T1unique.plot[,4] <- as.numeric(as.character(DFT2.SVs.338T.1515T1unique.plot[,4]))
circos.genomicLink(region1 = DFT2.SVs.338T.1515T1unique.plot[,c(1,2,2)],
                   region2 = DFT2.SVs.338T.1515T1unique.plot[,c(3,4,4)],
                   col = "red", lwd = 20,
                   rou = 0.69)

dev.off()

## subset unique 202T1 ones
DFT2.SVs.202T1 <- DFT2.SVs.summary[which(apply(DFT2.SVs.summary[,41,drop=F], 1, function(x){all(x == 1)}) == T),]
DFT2.SVs.202T1unique <- DFT2.SVs.202T1[which(apply(DFT2.SVs.202T1[,c(42:81)], 1, function(x){all(x == 0)}) == T),]
DFT2.SVs.202T1unique <- DFT2.SVs.202T1unique[c(1:7,9:21),] ## chromoplexy cluster

## plot DFT2 circos
pdf('Figure5E_circos_III_202T1.pdf', width = 20, height = 19)

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

## unique unclustered SVs
DFT2.SVs.202T1unique.plot <- as.data.frame(DFT2.SVs.202T1unique[,c(4,5,18,19)])
DFT2.SVs.202T1unique.plot[,2] <- as.numeric(as.character(DFT2.SVs.202T1unique.plot[,2]))
DFT2.SVs.202T1unique.plot[,4] <- as.numeric(as.character(DFT2.SVs.202T1unique.plot[,4]))
circos.genomicLink(region1 = DFT2.SVs.202T1unique.plot[,c(1,2,2)],
                   region2 = DFT2.SVs.202T1unique.plot[,c(3,4,4)],
                   col = "red", lwd = 20,
                   rou = 0.69)

dev.off()

## subset shared clade B ones
DFT2.SVs.cladeB <- DFT2.SVs.summary[which(apply(DFT2.SVs.summary[,c(64:81)], 1, function(x){all(x == 1)}) == T),]
DFT2.SVs.cladeBunique <- DFT2.SVs.cladeB[which(apply(DFT2.SVs.cladeB[,c(41:63)], 1, function(x){all(x == 0)}) == T),]
DFT2.SVs.cladeBunique <- DFT2.SVs.cladeBunique[c(2,3,4,5,11),] ## chromoplexy cluster

## plot DFT2 circos
pdf('Figure5E_circos_IV_cladeB.pdf', width = 20, height = 19)

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

## unique clustered SVs
DFT2.SVs.cladeBunique.plot <- as.data.frame(DFT2.SVs.cladeBunique[,c(4,5,18,19)])
DFT2.SVs.cladeBunique.plot[,2] <- as.numeric(as.character(DFT2.SVs.cladeBunique.plot[,2]))
DFT2.SVs.cladeBunique.plot[,4] <- as.numeric(as.character(DFT2.SVs.cladeBunique.plot[,4]))
circos.genomicLink(region1 = DFT2.SVs.cladeBunique.plot[,c(1,2,2)],
                   region2 = DFT2.SVs.cladeBunique.plot[,c(3,4,4)],
                   col = "red", lwd = 20,
                   rou = 0.69)

dev.off()

## clean up environment
rm(list=ls())
