## Figure 5F
## 'The evolution of two transmissible cancers in Tasmanian devils'
## Stammnitz et al., 2022
## maxrupsta@gmail

library(readxl)
library(treeio)
library(stringr)
library(ggtree)
library(ggplot2)

## set input path(s)
setwd('/Tables')


# Whole genome doubling in DFT1 tree #
######################################

## import pre vs. post-WGD mutation numbers
WGD.dates <- as.matrix(read_xlsx('Table-S6.xlsx', sheet = 5))
colnames(WGD.dates) <- WGD.dates[2,]
WGD.dates <- WGD.dates[-c(1:2),]
DFT1.WGD.dates <- WGD.dates[WGD.dates[,"LINEAGE"] == 'DFT1', ]

## number of substitutions
DFT1.truncal <- 1311 ### number of genome-wide substitutions in the DFT1 trunk
DFT1.truncal.median <- DFT1.truncal - 818 ### median expected number of somatic substitutions in the DFT1 trunk

## import and process BEAST MCMC draws on DFT1 mutation rate and date of origin (generated with Tracer)
masked.mSarHar11.1 <- c('A' = 951740967, 'C' = 540077344, 'G' = 539833602, 'T' = 952098282)
DFT1.beast.trees.MRCAs.rates <- cbind(read.table('BEAST_DFT1_MRCAs.txt', header = T),
                                      'rates' = read.table('BEAST_DFT1_rates.txt', header = T)[,2])
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
DFT1.beast.tree.hpd <- read.beast('BEAST_DFT1.mcc')
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
pdf("Figure5F_DFT1_tree_WGD.pdf", height = 7.5, width = 6)
p <- ggtree(DFT1.beast.tree.hpd, 
       root.position = 2018.13424657534 - max(DFT1.beast.tree.hpd@data[,'height']), 
       layout = 'rectangular', 
       lwd = 1, 
       ladderize = T) +
  scale_x_continuous(breaks = seq(f=1980, t=2020, by=10),
                     labels = seq(f=1980, t=2020, by=10),
                     limits = c(1978,2022)) +
  geom_point2(aes(subset = (node %in% grep(paste0('^',as.character(DFT1.WGD.dates[1,2]),'$'), DFT1.beast.tree.hpd@phylo$tip.label)),
                  x = x - c(as.numeric(DFT1.WGD.dates[1,"SAMPLING DATE"]) - as.numeric(DFT1.WGD.dates[1,"TETRAPLOIDISATION DATE [MEAN]"]))),
              shape = 16, size = 8,
              fill = 'cornflowerblue',
              col = 'cornflowerblue') +
  geom_point2(aes(subset = (node %in% grep(paste0('^',as.character(DFT1.WGD.dates[2,2]),'$'), DFT1.beast.tree.hpd@phylo$tip.label)),
                  x = x - c(as.numeric(DFT1.WGD.dates[2,"SAMPLING DATE"]) - as.numeric(DFT1.WGD.dates[2,"TETRAPLOIDISATION DATE [MEAN]"]))),
              shape = 16, size = 8,
              fill = 'cornflowerblue',
              col = 'cornflowerblue') +
  geom_point2(aes(subset = (node %in% grep(paste0('^',as.character(DFT1.WGD.dates[3,2]),'$'), DFT1.beast.tree.hpd@phylo$tip.label)),
                  x = x - c(as.numeric(DFT1.WGD.dates[3,"SAMPLING DATE"]) - as.numeric(DFT1.WGD.dates[3,"TETRAPLOIDISATION DATE [MEAN]"]))),
              shape = 16, size = 8,
              fill = 'cornflowerblue',
              col = 'cornflowerblue') +
  geom_point2(aes(subset = (node %in% grep(paste0('^',as.character(DFT1.WGD.dates[4,2]),'$'), DFT1.beast.tree.hpd@phylo$tip.label)),
                  x = x - c(as.numeric(DFT1.WGD.dates[4,"SAMPLING DATE"]) - as.numeric(DFT1.WGD.dates[4,"TETRAPLOIDISATION DATE [MEAN]"]))),
              shape = 16, size = 8,
              fill = 'cornflowerblue',
              col = 'cornflowerblue') +
  geom_point2(aes(subset = (node %in% grep(paste0('^',as.character(DFT1.WGD.dates[5,2]),'$'), DFT1.beast.tree.hpd@phylo$tip.label)),
                  x = x - c(as.numeric(DFT1.WGD.dates[5,"SAMPLING DATE"]) - as.numeric(DFT1.WGD.dates[5,"TETRAPLOIDISATION DATE [MEAN]"]))),
              shape = 16, size = 8,
              fill = 'cornflowerblue',
              col = 'cornflowerblue') +
  geom_point2(aes(subset = (node %in% grep(paste0('^',as.character(DFT1.WGD.dates[6,2]),'$'), DFT1.beast.tree.hpd@phylo$tip.label)),
                  x = x - c(as.numeric(DFT1.WGD.dates[6,"SAMPLING DATE"]) - as.numeric(DFT1.WGD.dates[6,"TETRAPLOIDISATION DATE [MEAN]"]))),
              shape = 16, size = 8,
              fill = 'cornflowerblue',
              col = 'cornflowerblue') +
  geom_point2(aes(subset = (node %in% grep(paste0('^',as.character(DFT1.WGD.dates[7,2]),'$'), DFT1.beast.tree.hpd@phylo$tip.label)),
                  x = x - c(as.numeric(DFT1.WGD.dates[7,"SAMPLING DATE"]) - as.numeric(DFT1.WGD.dates[7,"TETRAPLOIDISATION DATE [MEAN]"]))),
              shape = 16, size = 8,
              fill = 'cornflowerblue',
              col = 'cornflowerblue') +
  geom_point2(aes(subset = (node %in% grep(paste0('^',as.character(DFT1.WGD.dates[8,2]),'$'), DFT1.beast.tree.hpd@phylo$tip.label)),
                  x = x - c(as.numeric(DFT1.WGD.dates[8,"SAMPLING DATE"]) - as.numeric(DFT1.WGD.dates[8,"TETRAPLOIDISATION DATE [MEAN]"]))),
              shape = 16, size = 8,
              fill = 'cornflowerblue',
              col = 'cornflowerblue') +
  geom_point2(aes(subset = (node %in% grep(paste0('^',as.character(DFT1.WGD.dates[9,2]),'$'), DFT1.beast.tree.hpd@phylo$tip.label)),
                  x = x - c(as.numeric(DFT1.WGD.dates[9,"SAMPLING DATE"]) - as.numeric(DFT1.WGD.dates[9,"TETRAPLOIDISATION DATE [MEAN]"]))),
              shape = 16, size = 8,
              fill = 'cornflowerblue',
              col = 'cornflowerblue') +
  geom_point2(aes(subset = (node %in% grep(paste0('^',as.character(DFT1.WGD.dates[10,2]),'$'), DFT1.beast.tree.hpd@phylo$tip.label)),
                  x = x - c(as.numeric(DFT1.WGD.dates[10,"SAMPLING DATE"]) - as.numeric(DFT1.WGD.dates[10,"TETRAPLOIDISATION DATE [MEAN]"]))),
              shape = 16, size = 8,
              fill = 'cornflowerblue',
              col = 'cornflowerblue') +
  geom_point2(aes(subset = (node %in% grep(paste0('^',as.character(DFT1.WGD.dates[11,2]),'$'), DFT1.beast.tree.hpd@phylo$tip.label)),
                  x = x - c(as.numeric(DFT1.WGD.dates[11,"SAMPLING DATE"]) - as.numeric(DFT1.WGD.dates[11,"TETRAPLOIDISATION DATE [MEAN]"]))),
              shape = 16, size = 8,
              fill = 'cornflowerblue',
              col = 'cornflowerblue') +
  geom_point2(aes(subset = (node %in% grep(paste0('^',as.character(DFT1.WGD.dates[12,2]),'$'), DFT1.beast.tree.hpd@phylo$tip.label)),
                  x = x - c(as.numeric(DFT1.WGD.dates[12,"SAMPLING DATE"]) - as.numeric(DFT1.WGD.dates[12,"TETRAPLOIDISATION DATE [MEAN]"]))),
              shape = 16, size = 8,
              fill = 'cornflowerblue',
              col = 'cornflowerblue') +
  geom_point2(aes(subset = (node %in% grep(paste0('^',as.character(DFT1.WGD.dates[13,2]),'$'), DFT1.beast.tree.hpd@phylo$tip.label)),
                  x = x - c(as.numeric(DFT1.WGD.dates[13,"SAMPLING DATE"]) - as.numeric(DFT1.WGD.dates[13,"TETRAPLOIDISATION DATE [MEAN]"]))),
              shape = 16, size = 8,
              fill = 'cornflowerblue',
              col = 'cornflowerblue') +
  geom_point2(aes(subset = (node %in% grep(paste0('^',as.character(DFT1.WGD.dates[14,2]),'$'), DFT1.beast.tree.hpd@phylo$tip.label)),
                  x = x - c(as.numeric(DFT1.WGD.dates[14,"SAMPLING DATE"]) - as.numeric(DFT1.WGD.dates[14,"TETRAPLOIDISATION DATE [MEAN]"]))),
              shape = 16, size = 8,
              fill = 'cornflowerblue',
              col = 'cornflowerblue') +
  geom_point2(aes(subset = (node %in% MRCA(DFT1.beast.tree.hpd, c('2694Ta', '2694Tb')))),
              shape = 16, size = 8,
              fill = 'cornflowerblue',
              col = 'cornflowerblue') +
  geom_rootedge(rootedge = 2018.13424657534 - max(DFT1.beast.tree.hpd@data[,'height']) - mean(DFT1.beast.trees.MRCAs.rates[,'age at 818 singletons']),
                lwd = 1, col = 'black') +
  geom_point2(aes(subset = (node %in% 79),
                  x = x - c(2018.13424657534 - max(DFT1.beast.tree.hpd@data[,'height']) - mean(DFT1.beast.trees.MRCAs.rates[,'age at 818 singletons']))),
              shape = 16, size = 3, col = 'black') +
  theme_classic(base_size = 27) +
  labs(x = "Year", y = "") + 
  theme(axis.text = element_text(size = 27),
        axis.line.x = element_line(size = 1, color = 'black'),
        axis.text.y = element_blank(), 
        axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = 40, vjust = -0.6),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        text = element_text(family="Helvetica"),
        plot.margin = unit(c(t = 0, r = 20, b = 20, l = 0), unit ='pt'))
flip(p, clade_c, clade_c) %>% rotate(MRCA(DFT1.beast.tree.hpd, c("140T", "88T"))) %>% flip(clade_b, clade_d)
dev.off()

## clean up environment
rm(list=ls())


# Whole genome doubling in DFT2 tree #
######################################

## import pre vs. post-WGD mutation numbers
WGD.dates <- as.matrix(read_xlsx('Table-S6_v6.xlsx', sheet = 5))
colnames(WGD.dates) <- WGD.dates[2,]
WGD.dates <- WGD.dates[-c(1:2),]
DFT2.WGD.dates <- WGD.dates[WGD.dates[,"LINEAGE"] == 'DFT2', ]

## number of substitutions
DFT2.truncal <- 1335 ### number of genome-wide substitutions in the DFT2 trunk
DFT2.truncal.median <- DFT2.truncal - 818 ### median expected number of somatic substitutions in the DFT2 trunk

## import and process BEAST MCMC draws on DFT1 mutation rate and date of origin (generated with Tracer)
masked.mSarHar11.1 <- c('A' = 951740967, 'C' = 540077344, 'G' = 539833602, 'T' = 952098282)
DFT2.beast.trees.MRCAs.rates <- cbind(read.table('BEAST_DFT2_MRCAs.txt', header = T),
                                      'rates' = read.table('BEAST_DFT2_rates.txt', header = T)[,2])
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
DFT2.beast.tree.hpd <- read.beast('BEAST_DFT2.mcc')
DFT2.beast.tree.hpd@data[which.max(unlist(DFT2.beast.tree.hpd@data[,'height'])),"height_0.95_HPD"][[1]][[1]][2] <- 2018.66575342466 - 
  as.numeric(quantile(DFT2.beast.trees.MRCAs.rates[,'age at 0 singletons'], probs = c(0.05)))

## mark DFT2 clades
clade_a <- MRCA(DFT2.beast.tree.hpd, c('339T', '638T1'))
clade_b <- MRCA(DFT2.beast.tree.hpd, c('1529T2', '1545T2'))

## plot
pdf("Figure5F_DFT2_tree_WGD.pdf", height = 7.5, width = 6)
ggtree(DFT2.beast.tree.hpd, 
       root.position = 2018.66575342466 - max(DFT2.beast.tree.hpd@data[,'height']), 
       layout = 'rectangular', 
       lwd = 1, 
       ladderize = T) +
  scale_x_continuous(breaks = seq(f=2011, t=2019, by=2),
                     labels = seq(f=2011, t=2019, by=2),
                     limits = c(2011,2019)) +
  geom_point2(aes(subset = (node %in% grep(as.character(DFT2.WGD.dates[1,2]), DFT2.beast.tree.hpd@phylo$tip.label)),
                  x = x - c(as.numeric(DFT2.WGD.dates[1,"SAMPLING DATE"]) - as.numeric(DFT2.WGD.dates[1,"TETRAPLOIDISATION DATE [MEAN]"]))),
              shape = 16, size = 8,
              fill = 'red',
              col = 'red') +
  geom_point2(aes(subset = (node %in% grep(as.character(DFT2.WGD.dates[2,2]), DFT2.beast.tree.hpd@phylo$tip.label)),
                  x = x - c(as.numeric(DFT2.WGD.dates[2,"SAMPLING DATE"]) - as.numeric(DFT2.WGD.dates[2,"TETRAPLOIDISATION DATE [MEAN]"]))),
              shape = 16, size = 8,
              fill = 'red',
              col = 'red') +
  geom_point2(aes(subset = (node %in% grep(as.character(DFT2.WGD.dates[3,2]), DFT2.beast.tree.hpd@phylo$tip.label)),
                  x = x - c(as.numeric(DFT2.WGD.dates[3,"SAMPLING DATE"]) - as.numeric(DFT2.WGD.dates[3,"TETRAPLOIDISATION DATE [MEAN]"]))),
              shape = 16, size = 8,
              fill = 'red',
              col = 'red') +
  geom_rootedge(rootedge = 2018.66575342466 - max(DFT2.beast.tree.hpd@data[,'height']) - mean(DFT2.beast.trees.MRCAs.rates[,'age at 818 singletons']),
                lwd = 1, col = 'black') +
  geom_point2(aes(subset = (node %in% 42),
                  x = x - c(2018.66575342466 - max(DFT2.beast.tree.hpd@data[,'height']) - mean(DFT2.beast.trees.MRCAs.rates[,'age at 818 singletons']))),
              shape = 16, size = 3, col = 'black') +
  theme_classic(base_size = 27) +
  labs(x = "Year", y = "") + 
  theme(axis.text = element_text(size = 27),
        axis.line.x = element_line(size = 1, color = 'black'),
        axis.text.y = element_blank(), 
        axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = 40, vjust = -0.6),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        text = element_text(family="Helvetica"),
        plot.margin = unit(c(t = 0, r = 20, b = 20, l = 0), unit ='pt'))
dev.off()

## clean up environment
rm(list=ls())
