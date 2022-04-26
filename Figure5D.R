## Figure 5D
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


# Examples of (late) DFT1 chromoplexy - tree trajectory #
#########################################################

## number of substitutions
DFT1.truncal <- 1311 ### number of genome-wide substitutions in the DFT1 trunk
DFT1.truncal.median <- DFT1.truncal - 818 ### median expected number of somatic substitutions in the DFT1 trunk

## import and process BEAST MCMC draws on DFT1 mutation rate and date of origin (generated with Tracer)
masked.mSarHar11.1 <- c('A' = 951740967, 'C' = 540077344, 'G' = 539833602, 'T' = 952098282)
DFT1.beast.trees.MRCAs.rates <- cbind(read.table('/Users/mstammnitz/Desktop/DFT_evolution/doc/manuscripts/The Evolutionary History of Two Transmissible Cancers in Tasmanian Devils/Tables/v6/Supplementary_data/BEAST_DFT1_MRCAs.txt', header = T),
                                      'rates' = read.table('/Users/mstammnitz/Desktop/DFT_evolution/doc/manuscripts/The Evolutionary History of Two Transmissible Cancers in Tasmanian Devils/Tables/v6/Supplementary_data/BEAST_DFT1_rates.txt', header = T)[,2])
DFT1.beast.trees.MRCAs.rates[,'rates'] <- DFT1.beast.trees.MRCAs.rates[,'rates']*
  sum(masked.mSarHar11.1) ### convert to genome-wide mutations per year
DFT1.beast.trees.MRCAs.rates <- cbind(DFT1.beast.trees.MRCAs.rates,
                                      'age at 0 singletons' = rep(NA, nrow(DFT1.beast.trees.MRCAs.rates)),
                                      'age at 818 singletons' = rep(NA, nrow(DFT1.beast.trees.MRCAs.rates)))
DFT1.beast.trees.MRCAs.rates[,'age at 0 singletons'] <- DFT1.beast.trees.MRCAs.rates[,'age.root.'] - 
  c(DFT1.truncal/DFT1.beast.trees.MRCAs.rates[,'rates'])
DFT1.beast.trees.MRCAs.rates[,'age at 818 singletons'] <- DFT1.beast.trees.MRCAs.rates[,'age.root.'] - 
  c(DFT1.truncal.median/DFT1.beast.trees.MRCAs.rates[,'rates'])

## import and process DFT1 maximum clade consensus tree
DFT1.beast.tree.hpd <- read.beast('/Users/mstammnitz/Desktop/DFT_evolution/doc/manuscripts/The Evolutionary History of Two Transmissible Cancers in Tasmanian Devils/Tables/v6/Supplementary_data/BEAST_DFT1.mcc')
DFT1.beast.tree.hpd@data[which.max(unlist(DFT1.beast.tree.hpd@data[,'height'])),'height_0.95_HPD'][[1]][[1]][2] <- 2018.13424657534 - 
  as.numeric(quantile(DFT1.beast.trees.MRCAs.rates[,'age at 0 singletons'], probs = c(0.05)))

## mark DFT1 clades
clade_a1 <- MRCA(DFT1.beast.tree.hpd, c('78T', '88T'))
clade_a2 <- MRCA(DFT1.beast.tree.hpd, c('54T1', '372T1'))
clade_b <- MRCA(DFT1.beast.tree.hpd, c('1191T1', '1059T2'))
clade_c <- MRCA(DFT1.beast.tree.hpd, c('140T', '46T1'))
clade_cnot1 <- MRCA(DFT1.beast.tree.hpd, c("140T", "356T1"))
clade_d <- MRCA(DFT1.beast.tree.hpd, '102T2')
clade_e <- MRCA(DFT1.beast.tree.hpd, '377T1')

## plot
pdf("Figure5D_DFT1_chromoplexy_tree_trajectory.pdf",
    width = 4, height = 8)

## hack in the tip edge colors
edge.colors <- rep('black', DFT1.beast.tree.hpd@phylo$Nnode+length(DFT1.beast.tree.hpd@phylo$tip.label))
special.edge.color <- 'cornflowerblue'

## hack in the internal edge colors
internal.edges1 <- MRCA(DFT1.beast.tree.hpd, c('421T1', '451T1'))
edge.colors[internal.edges1] <- special.edge.color
internal.edges2 <- MRCA(DFT1.beast.tree.hpd, c('54T1', '1071T1'))
edge.colors[internal.edges2] <- special.edge.color
internal.edges3 <- MRCA(DFT1.beast.tree.hpd, c('54T1', '398T1'))
edge.colors[internal.edges3] <- special.edge.color
internal.edges4 <- MRCA(DFT1.beast.tree.hpd, c('1011T1', '158T'))
edge.colors[internal.edges4] <- special.edge.color
internal.edges5 <- MRCA(DFT1.beast.tree.hpd, c('54T1', '158T'))
edge.colors[internal.edges5] <- special.edge.color
internal.edges6 <- MRCA(DFT1.beast.tree.hpd, c('54T1', '378T1'))
edge.colors[internal.edges6] <- special.edge.color
internal.edges7 <- MRCA(DFT1.beast.tree.hpd, c('377T3', '378T1'))
edge.colors[internal.edges7] <- special.edge.color
internal.edges8 <- MRCA(DFT1.beast.tree.hpd, c('54T1', '372T1'))
edge.colors[internal.edges8] <- special.edge.color
internal.edges9 <- MRCA(DFT1.beast.tree.hpd, c('377T3', '18T'))
edge.colors[internal.edges9] <- special.edge.color
internal.edges10 <- MRCA(DFT1.beast.tree.hpd, c('377T3', '421T1'))
edge.colors[internal.edges10] <- special.edge.color
internal.edges11 <- MRCA(DFT1.beast.tree.hpd, c('377T3', '84T1'))
edge.colors[internal.edges11] <- special.edge.color
internal.edges12 <- MRCA(DFT1.beast.tree.hpd, c('377T3', '140T'))
edge.colors[internal.edges12] <- special.edge.color

## summarise edge colors in co-supplied dataframe
d <- data.frame(node=1:(DFT1.beast.tree.hpd@phylo$Nnode+length(DFT1.beast.tree.hpd@phylo$tip.label)), 
                color = edge.colors)
p <- ggtree(DFT1.beast.tree.hpd, 
       root.position = 2018.13424657534 - max(DFT1.beast.tree.hpd@data[,'height']), 
       layout = 'rectangular', 
       lwd = 1, 
       ladderize = T) %<+% d + aes(color=I(color)) + 
  scale_x_continuous(breaks = seq(f=1980, t=2020, by=20),
                     labels = seq(f=1980, t=2020, by=20),
                     limits = c(1975,2025)) +
  geom_rootedge(rootedge = 2018.13424657534 - max(DFT1.beast.tree.hpd@data[,'height']) - mean(DFT1.beast.trees.MRCAs.rates[,'age at 818 singletons']),
                lwd = 1, col = 'black') +
  geom_point2(aes(subset = (node %in% 79),
                  x = x - c(2018.13424657534 - max(DFT1.beast.tree.hpd@data[,'height']) - mean(DFT1.beast.trees.MRCAs.rates[,'age at 818 singletons']))),
              shape = 16, size = 3, col = 'black') +
  geom_point2(aes(subset=(node %in% c(internal.edges1,internal.edges2,internal.edges3,internal.edges4))), 
              shape = 16, size = 7, 
              fill = 'cornflowerblue',
              col = 'cornflowerblue') +
  theme_classic(base_size = 20) +
  theme(axis.text=element_blank(), 
        legend.position = "none",
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(), 
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        line = element_blank())
flip(p, clade_c, clade_c) %>% rotate(MRCA(DFT1.beast.tree.hpd, c("140T", "88T"))) %>% flip(clade_b, clade_d)

dev.off()

## clean up environment
rm(list=ls())


# Examples of (late) DFT1 chromoplexy - circos plots #
######################################################

## import structural variant table
DFT1.SVs <- as.matrix(read_xlsx('/Users/mstammnitz/Desktop/DFT_evolution/doc/manuscripts/The Evolutionary History of Two Transmissible Cancers in Tasmanian Devils/Tables/v6/Table-S5_v6.xlsx', sheet = 2))
colnames(DFT1.SVs) <- as.character(DFT1.SVs[2,])
DFT1.SVs <- DFT1.SVs[-c(1:2),]
DFT1.SVs.summary <- DFT1.SVs
DFT1.SVs.MSG <- DFT1.SVs.summary[grep('MSG', DFT1.SVs.summary[,"CALLER"]),]
DFT1.SVs.MSG <- cbind(DFT1.SVs.MSG[,1:40], apply(DFT1.SVs.MSG[,41:105], 2, function(x){out <- as.numeric(str_split_fixed(x, '/', 2)[,1]); out[which(out >= 1)] <- 1; return(out)}))
DFT1.SVs.SvABA <- DFT1.SVs.summary[-grep('MSG', DFT1.SVs.summary[,"CALLER"]),]
DFT1.SVs.SvABA <- cbind(DFT1.SVs.SvABA[,1:40],apply(DFT1.SVs.SvABA[,41:105], 2, function(x){out <- as.numeric(x); out[which(out < 3)] <- 0; out[which(out >= 3)] <- 1; return(out)}))
DFT1.SVs.summary <- rbind(DFT1.SVs.MSG, DFT1.SVs.SvABA)
DFT1.SVs.summary <- DFT1.SVs.summary[order(as.character(DFT1.SVs.summary[,"CHROM1"]), as.numeric(DFT1.SVs.summary[,"ORIGINAL POS1"])),]
DFT1.SVs.summary <- DFT1.SVs.summary[,-grep('^458T1$', colnames(DFT1.SVs.summary))]
DFT1.SVs.summary <- DFT1.SVs.summary[which(apply(DFT1.SVs.summary[,41:104], 1, function(x){all(x == 0)}) == F),]

## chromosomes
chromosome.ranges <- matrix(0, ncol = 2, nrow = 7)
chromosome.ranges[,1] <- 1
chromosome.ranges[,2] <- c(716413629, 662751787, 611347268, 464895054, 288121652, 254895979, 83081154)
rownames(chromosome.ranges) <- c("1", "2", "3", "4", "5", "6", "X")
chromosome.ranges.GR <- GRanges(seqnames = rownames(chromosome.ranges), 
                                ranges = IRanges(start = as.numeric(chromosome.ranges[,1]),
                                                 end = as.numeric(chromosome.ranges[,2])))

## subset T1 ones
DFT1.SVs.T1 <- DFT1.SVs.summary[which(apply(DFT1.SVs.summary[,c(53:71)], 1, function(x){all(x == 1)}) == T),]
DFT1.SVs.T1.unique <- DFT1.SVs.T1[which(apply(DFT1.SVs.T1[,c(41:52,72:104)], 1, function(x){all(x == 0)}) == T),]
DFT1.SVs.T1.unique <- DFT1.SVs.T1.unique[c(1, 2, 3, 7, 10, 11, 12),]

## plot
pdf('Figure5D_T1_circos.pdf', 
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

## add T1 - unique clustered SVs
DFT1.SVs.T1.unique.plot <- as.data.frame(DFT1.SVs.T1.unique[,c(4,5,18,19)])
DFT1.SVs.T1.unique.plot[,2] <- as.numeric(as.character(DFT1.SVs.T1.unique.plot[,2]))
DFT1.SVs.T1.unique.plot[,4] <- as.numeric(as.character(DFT1.SVs.T1.unique.plot[,4]))
circos.genomicLink(region1 = DFT1.SVs.T1.unique.plot[,c(1,2,2)],
                   region2 = DFT1.SVs.T1.unique.plot[,c(3,4,4)],
                   col = "cornflowerblue", lwd = 20,
                   rou = 0.69)

dev.off()

## subset T2 ones
DFT1.SVs.T2 <- DFT1.SVs.summary[which(apply(DFT1.SVs.summary[,c(63:71)], 1, function(x){all(x == 1)}) == T),]
DFT1.SVs.T2.unique <- DFT1.SVs.T2[which(apply(DFT1.SVs.T2[,c(41:62,72:104)], 1, function(x){all(x == 0)}) == T),]
DFT1.SVs.T2.unique <- DFT1.SVs.T2.unique[c(1:5),]

## plot
pdf('Figure5D_T2_circos.pdf',  width = 20, height = 19)

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

## add T2 - unique clustered SVs
DFT1.SVs.T2.unique.plot <- as.data.frame(DFT1.SVs.T2.unique[,c(4,5,18,19)])
DFT1.SVs.T2.unique.plot[,2] <- as.numeric(as.character(DFT1.SVs.T2.unique.plot[,2]))
DFT1.SVs.T2.unique.plot[,4] <- as.numeric(as.character(DFT1.SVs.T2.unique.plot[,4]))
circos.genomicLink(region1 = DFT1.SVs.T2.unique.plot[,c(1,2,2)],
                   region2 = DFT1.SVs.T2.unique.plot[,c(3,4,4)],
                   col = "cornflowerblue", lwd = 20,
                   rou = 0.69)

dev.off()

## subset T3 ones
DFT1.SVs.T3 <- DFT1.SVs.summary[which(apply(DFT1.SVs.summary[,c(66:71)], 1, function(x){all(x == 1)}) == T),]
DFT1.SVs.T3.unique <- DFT1.SVs.T3[which(apply(DFT1.SVs.T3[,c(41:65,72:104)], 1, function(x){all(x == 0)}) == T),]

## plot DFT1 circos
pdf('Figure5D_T3_circos.pdf', 
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

## add T3 - unique clustered SVs
DFT1.SVs.T3.unique.plot <- as.data.frame(DFT1.SVs.T3.unique[,c(4,5,18,19)])
DFT1.SVs.T3.unique.plot[,2] <- as.numeric(as.character(DFT1.SVs.T3.unique.plot[,2]))
DFT1.SVs.T3.unique.plot[,4] <- as.numeric(as.character(DFT1.SVs.T3.unique.plot[,4]))
circos.genomicLink(region1 = DFT1.SVs.T3.unique.plot[,c(1,2,2)],
                   region2 = DFT1.SVs.T3.unique.plot[,c(3,4,4)],
                   col = "cornflowerblue", lwd = 20,
                   rou = 0.69)

dev.off()

## subset T4 ones
DFT1.SVs.T4 <- DFT1.SVs.summary[which(apply(DFT1.SVs.summary[,c(68,69)], 1, function(x){all(x == 1)}) == T),]
DFT1.SVs.T4.unique <- DFT1.SVs.T4[which(apply(DFT1.SVs.T4[,c(41:67,70:104)], 1, function(x){all(x == 0)}) == T),]
DFT1.SVs.T4.unique <- DFT1.SVs.T4.unique[c(1,2,3,4,5,7,8,12),]

## plot DFT1 circos
pdf('Figure5D_T4_circos.pdf', 
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

## add T4 - unique clustered SVs
DFT1.SVs.T4.unique.plot <- as.data.frame(DFT1.SVs.T4.unique[,c(4,5,18,19)])
DFT1.SVs.T4.unique.plot[,2] <- as.numeric(as.character(DFT1.SVs.T4.unique.plot[,2]))
DFT1.SVs.T4.unique.plot[,4] <- as.numeric(as.character(DFT1.SVs.T4.unique.plot[,4]))
circos.genomicLink(region1 = DFT1.SVs.T4.unique.plot[,c(1,2,2)],
                   region2 = DFT1.SVs.T4.unique.plot[,c(3,4,4)],
                   col = "cornflowerblue", lwd = 20,
                   rou = 0.69)

dev.off()

## clean up environment
rm(list=ls())