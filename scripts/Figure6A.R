## Figure 6A
## 'The evolution of two transmissible cancers in Tasmanian devils'
## Stammnitz et al., 2022
## maxrupsta@gmail

library(treeio)
library(ggtree)
library(ggplot2)
library(caper)

## set input path(s)
setwd('/Tables')


# Drivers on DFT1 trees #
#########################

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

### LZTR1
## (1) 1st hit in ancestral node: 79

pdf("Figure6A_DFT1_I_LZTR1.pdf",
    width = 4, height = 8)
tree2 <- groupClade(DFT1.beast.tree.hpd, 79)
p <- ggtree(tree2, aes(color = group),
            root.position = 2018.13424657534 - max(DFT1.beast.tree.hpd@data[,'height']), 
            layout = 'rectangular', 
            lwd = 1, 
            ladderize = T) + 
  theme(legend.position='none') +
  scale_color_manual(values=c("cornflowerblue", 'black')) + 
  scale_x_continuous(breaks = seq(f=1980, t=2020, by=10),
                     labels = seq(f=1980, t=2020, by=10),
                     limits = c(1980,2024)) +
  geom_rootedge(rootedge = 2018.13424657534 - max(DFT1.beast.tree.hpd@data[,'height']) - mean(DFT1.beast.trees.MRCAs.rates[,'age at 818 singletons']),
                lwd = 1, col = 'cornflowerblue') +
  geom_point2(aes(subset = (node %in% 79),
                  x = x - c(2018.13424657534 - max(DFT1.beast.tree.hpd@data[,'height']) - mean(DFT1.beast.trees.MRCAs.rates[,'age at 818 singletons']))),
              shape = 16, size = 3, col = 'cornflowerblue') +
  geom_point2(aes(subset=(node %in% c(79))), shape = 25, size = 10, 
              fill = alpha('cornflowerblue', 0.8),
              col = alpha('cornflowerblue', 0.8))
ggtree::flip(p, clade_c, clade_c) %>% ggtree::rotate(MRCA(DFT1.beast.tree.hpd, c("140T", "88T"))) %>% ggtree::flip(clade_b, clade_d)
dev.off()

### MGA:
## (1) ancestral CN2 -> CN1: all but 377T1, node: 80
## (2) 1st hit in 377T1, node: 46
## (3) 2nd hit in 52T2 & 88T, node: 148
## (4) 2nd hit in 458T1, node: 57
## (5) 2nd hit in 356T1, node: 41
## (6) 2nd hit in 78T1, node: 69
## (7) 2nd hit in 63T3, node: 68

pdf("Figure6A_DFT1_II_MGA.pdf",
    width = 4, height = 8)
tree2 <- groupClade(DFT1.beast.tree.hpd, 80)
p <- ggtree(tree2, 
            root.position = 2018.13424657534 - max(DFT1.beast.tree.hpd@data[,'height']), 
            layout = 'rectangular', 
            lwd = 1, 
            ladderize = T, aes(color = I(group))) +
  theme(legend.position='none') +
  scale_color_manual(values=c("cornflowerblue", "cornflowerblue")) + 
  scale_x_continuous(breaks = seq(f=1980, t=2020, by=10),
                     labels = seq(f=1980, t=2020, by=10),
                     limits = c(1980,2024)) +
  geom_rootedge(rootedge = 2018.13424657534 - max(DFT1.beast.tree.hpd@data[,'height']) - mean(DFT1.beast.trees.MRCAs.rates[,'age at 818 singletons']),
                lwd = 1, col = 'black') +
  geom_point2(aes(subset = (node %in% 79),
                  x = x - c(2018.13424657534 - max(DFT1.beast.tree.hpd@data[,'height']) - mean(DFT1.beast.trees.MRCAs.rates[,'age at 818 singletons']))),
              shape = 16, size = 3, col = 'black') +
  geom_point2(aes(subset=(node %in% c(41, 46, 57, 68, 69, 80, 148))), shape = 25, size = 10, 
              fill = alpha('cornflowerblue', 0.8),
              col = alpha('cornflowerblue', 0.8))
ggtree::flip(p, clade_c, clade_c) %>% ggtree::rotate(MRCA(DFT1.beast.tree.hpd, c("140T", "88T"))) %>% ggtree::flip(clade_b, clade_d)
dev.off()

### PDGFRB
## (1) gain in 356T1, node: 41
## (2) post-tetraploid gain in 1191T1, node: 6
## (3) gain(s) in 451T1, node: 56
## (4) gain(s) in 236T3, node: 30
## (5) gain(s) in 52T2, 88T, node: 148
## (6) gain(s) in 52T2, node: 65
## (7) gain(s) in 88T, node: 77

pdf("Figure6A_DFT1_III_PDGFRB.pdf",
    width = 4, height = 8)

## hack in the tip edge colors
tip.edges <- grep('356T1|1191T1|451T1|236T3|52T2|88T', DFT1.beast.tree.hpd@phylo$tip.label)
edge.colors <- rep('black', DFT1.beast.tree.hpd@phylo$Nnode+length(DFT1.beast.tree.hpd@phylo$tip.label))
edge.colors[tip.edges] <- "cornflowerblue"

## hack in the internal edge colors
internal.edges <- MRCA(DFT1.beast.tree.hpd, c('52T2', '88T'))
edge.colors[internal.edges] <- "cornflowerblue"

## summarise edge colors in co-supplied dataframe
d <- data.frame(node=1:(DFT1.beast.tree.hpd@phylo$Nnode+length(DFT1.beast.tree.hpd@phylo$tip.label)), 
                color = edge.colors)

## plot actual tree
p <- ggtree(DFT1.beast.tree.hpd, 
       root.position = 2018.13424657534 - max(DFT1.beast.tree.hpd@data[,'height']), 
       layout = 'rectangular', 
       lwd = 1, 
       ladderize = T) %<+% d + aes(color=I(color)) + 
  scale_x_continuous(breaks = seq(f=1980, t=2020, by=10),
                     labels = seq(f=1980, t=2020, by=10),
                     limits = c(1980,2024)) +
  geom_rootedge(rootedge = 2018.13424657534 - max(DFT1.beast.tree.hpd@data[,'height']) - mean(DFT1.beast.trees.MRCAs.rates[,'age at 818 singletons']),
                lwd = 1, col = 'black') +
  geom_point2(aes(subset = (node %in% 79),
                  x = x - c(2018.13424657534 - max(DFT1.beast.tree.hpd@data[,'height']) - mean(DFT1.beast.trees.MRCAs.rates[,'age at 818 singletons']))),
              shape = 16, size = 3, col = 'black') +
  geom_point2(aes(subset=(node %in% c(6, 30, 41, 56, 65, 77, 148))), shape = 24, size = 10, 
              fill = alpha('cornflowerblue', 0.8),
              col = alpha('cornflowerblue', 0.8))
ggtree::flip(p, clade_c, clade_c) %>% ggtree::rotate(MRCA(DFT1.beast.tree.hpd, c("140T", "88T"))) %>% ggtree::flip(clade_b, clade_d)
dev.off()

### Chr4:301.6Mb
## (1) gain in 998T1, node: 78
## (2) two-step gain in 155T, 155Ta, 155Tb, 54T1, node: 141
## (3) extra gain in 54T1, node: 66
## (4) extra gains in 155Ta, 155Tb, node: 143
## (5) gains in 1058T1, 46T1, node: 123
pdf("Figure6A_DFT1_IV_Chr4_302Mb.pdf",
    width = 4, height = 8)

## hack in the tip edge colors
tip.edges <- grep('^998T1$|^155T$|^155Ta$|^155Tb$|^54T1$|^1058T1$|^46T1$', DFT1.beast.tree.hpd@phylo$tip.label)
edge.colors <- rep('black', DFT1.beast.tree.hpd@phylo$Nnode+length(DFT1.beast.tree.hpd@phylo$tip.label))
edge.colors[tip.edges] <- "cornflowerblue"

## hack in the internal edge colors
internal.edges.1 <- MRCA(DFT1.beast.tree.hpd, c('54T1', '155T'))
edge.colors[internal.edges.1] <- "cornflowerblue"
internal.edges.2 <- MRCA(DFT1.beast.tree.hpd, c('155T', '155Ta'))
edge.colors[internal.edges.2] <- "cornflowerblue"
internal.edges.3 <- MRCA(DFT1.beast.tree.hpd, c('155Ta', '155Tb'))
edge.colors[internal.edges.3] <- "cornflowerblue"
internal.edges.4 <- MRCA(DFT1.beast.tree.hpd, c('155T', '155Tb'))
edge.colors[internal.edges.4] <- "cornflowerblue"
internal.edges.5 <- MRCA(DFT1.beast.tree.hpd, c('1058T1', '46T1'))
edge.colors[internal.edges.5] <- "cornflowerblue"

## summarise edge colors in co-supplied dataframe
d <- data.frame(node=1:(DFT1.beast.tree.hpd@phylo$Nnode+length(DFT1.beast.tree.hpd@phylo$tip.label)), 
                color = edge.colors)

## plot actual tree
p <- ggtree(DFT1.beast.tree.hpd, 
       root.position = 2018.13424657534 - max(DFT1.beast.tree.hpd@data[,'height']), 
       layout = 'rectangular', 
       lwd = 1, 
       ladderize = T) %<+% d + aes(color=I(color)) + 
  scale_x_continuous(breaks = seq(f=1980, t=2020, by=10),
                     labels = seq(f=1980, t=2020, by=10),
                     limits = c(1980,2024)) +
  geom_rootedge(rootedge = 2018.13424657534 - max(DFT1.beast.tree.hpd@data[,'height']) - mean(DFT1.beast.trees.MRCAs.rates[,'age at 818 singletons']),
                lwd = 1, col = 'black') +
  geom_point2(aes(subset = (node %in% 79),
                  x = x - c(2018.13424657534 - max(DFT1.beast.tree.hpd@data[,'height']) - mean(DFT1.beast.trees.MRCAs.rates[,'age at 818 singletons']))),
              shape = 16, size = 3, col = 'black') +
  geom_point2(aes(subset=(node %in% c(66, 78, 123, 141, 143))), shape = 24, size = 10, 
              fill = alpha('cornflowerblue', 0.8), 
              col = alpha('cornflowerblue', 0.8))
ggtree::flip(p, clade_c, clade_c) %>% ggtree::rotate(MRCA(DFT1.beast.tree.hpd, c("140T", "88T"))) %>% ggtree::flip(clade_b, clade_d)
dev.off()

### HMGA2
## (1) gain(s) in 18T, node: 26
## (2) two-step gain in 1011T1, 158T, node: 144
## also add gains in 86T, 88T, 134T6 as well?
pdf("Figure6A_DFT1_V_HMGA2.pdf",
    width = 4, height = 8)

## hack in the tip edge colors
tip.edges <- grep('^18T$|^1011T1$|^158T$', DFT1.beast.tree.hpd@phylo$tip.label)
edge.colors <- rep('black', DFT1.beast.tree.hpd@phylo$Nnode+length(DFT1.beast.tree.hpd@phylo$tip.label))
edge.colors[tip.edges] <- "cornflowerblue"

## hack in the internal edge colors
internal.edges <- MRCA(DFT1.beast.tree.hpd, c('1011T1', '158T'))
edge.colors[internal.edges] <- "cornflowerblue"

## summarise edge colors in co-supplied dataframe
d <- data.frame(node=1:(DFT1.beast.tree.hpd@phylo$Nnode+length(DFT1.beast.tree.hpd@phylo$tip.label)), 
                color = edge.colors)

## plot actual tree
p <- ggtree(DFT1.beast.tree.hpd, 
            root.position = 2018.13424657534 - max(DFT1.beast.tree.hpd@data[,'height']), 
            layout = 'rectangular', 
            lwd = 1, 
            ladderize = T) %<+% d + aes(color=I(color)) + 
  scale_x_continuous(breaks = seq(f=1980, t=2020, by=10),
                     labels = seq(f=1980, t=2020, by=10),
                     limits = c(1980,2024)) +
  geom_rootedge(rootedge = 2018.13424657534 - max(DFT1.beast.tree.hpd@data[,'height']) - mean(DFT1.beast.trees.MRCAs.rates[,'age at 818 singletons']),
                lwd = 1, col = 'black') +
  geom_point2(aes(subset = (node %in% 79),
                  x = x - c(2018.13424657534 - max(DFT1.beast.tree.hpd@data[,'height']) - mean(DFT1.beast.trees.MRCAs.rates[,'age at 818 singletons']))),
              shape = 16, size = 3, col = 'black') +
  geom_point2(aes(subset=(node %in% c(26, 144))), shape = 24, size = 10, 
              fill = alpha('cornflowerblue', 0.8), 
              col = alpha('cornflowerblue', 0.8))
ggtree::flip(p, clade_c, clade_c) %>% ggtree::rotate(MRCA(DFT1.beast.tree.hpd, c("140T", "88T"))) %>% ggtree::flip(clade_b, clade_d)
dev.off()

## clean up environment
rm(list=ls())


# Drivers on DFT2 trees #
#########################

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

## DFT2: highlight driver gene events

### PDGFRA
## (1) gain in all, node: 42
## (2) gain in all clade As, node: 43
## (3) gains in all 638T1, node: 38
## (4) gains in all 1515T1, node: 6
## (5) gains in all 1524T1, node: 7
## (6) gains in all 812T1, node: 40
pdf("Figure6A_DFT2_I_PDGFRA.pdf",
    width = 4, height = 8)
ggtree(DFT2.beast.tree.hpd, 
       root.position = 2018.66575342466 - max(DFT2.beast.tree.hpd@data[,'height']), 
       layout = 'rectangular', 
       lwd = 1, 
       col = 'red',
       ladderize = T) + 
  scale_x_continuous(breaks = seq(f=2008, t=2020, by=3),
                     labels = seq(f=2008, t=2020, by=3),
                     limits = c(2008,2020)) +
  geom_rootedge(rootedge = 2018.66575342466 - max(DFT2.beast.tree.hpd@data[,'height']) - mean(DFT2.beast.trees.MRCAs.rates[,'age at 818 singletons']),
                lwd = 1, col = 'red') +
  geom_point2(aes(subset = (node %in% 42),
                  x = x - c(2018.66575342466 - max(DFT2.beast.tree.hpd@data[,'height']) - mean(DFT2.beast.trees.MRCAs.rates[,'age at 818 singletons']))),
              shape = 16, size = 3, col = 'red') +
  geom_point2(aes(subset=(node %in% c(6,7,38,40,42,43))), shape = 24, size = 10, 
              fill = alpha('red', 0.8),
              col = alpha('red', 0.8))
dev.off()

### chrY
## (1) loss in all clade Bs, node: 65
## (2) loss in 1538T1,1538T2,1538T3,637T1,637T2, node: 61
## (3) loss in 1524T1, 1548T1: 55
## (4) loss in 1525T1, node: 8
## (5) loss in 1528T, node: 9
pdf("Figure6A_DFT2_II_chrY.pdf",
    width = 4, height = 8)

## hack in the tip edge colors
tip.edges <- grep('^1524T1$|^1548T1$|^1528T$|^1525T1$|^1538T1$|^1538T2$|^1538T3$|^637T1$|^637T2$', DFT2.beast.tree.hpd@phylo$tip.label)
edge.colors <- rep('black', DFT2.beast.tree.hpd@phylo$Nnode+length(DFT2.beast.tree.hpd@phylo$tip.label))
edge.colors[tip.edges] <- "red"

## hack in the internal edge colors
internal.edges.1 <- clade.members(65, DFT2.beast.tree.hpd@phylo, include.nodes = T)
internal.edges.1 <- sort(as.numeric(unlist(internal.edges.1)))
edge.colors[internal.edges.1] <- "red"
internal.edges.2 <- clade.members(61, DFT2.beast.tree.hpd@phylo, include.nodes = T)
internal.edges.2 <- sort(as.numeric(unlist(internal.edges.2)))
edge.colors[internal.edges.2] <- "red"
internal.edges.3 <- clade.members(55, DFT2.beast.tree.hpd@phylo, include.nodes = T)
internal.edges.3 <- sort(as.numeric(unlist(internal.edges.3)))
edge.colors[internal.edges.3] <- "red"

## summarise edge colors in co-supplied dataframe
d <- data.frame(node=1:(DFT2.beast.tree.hpd@phylo$Nnode+length(DFT2.beast.tree.hpd@phylo$tip.label)), 
                color = edge.colors)

## plot actual tree
ggtree(DFT2.beast.tree.hpd, 
       root.position = 2018.66575342466 - max(DFT2.beast.tree.hpd@data[,'height']), 
       layout = 'rectangular', 
       lwd = 1, 
       ladderize = T) %<+% d + aes(color=I(color)) +
  scale_x_continuous(breaks = seq(f=2008, t=2020, by=3),
                     labels = seq(f=2008, t=2020, by=3),
                     limits = c(2008,2020)) +
  geom_rootedge(rootedge = 2018.66575342466 - max(DFT2.beast.tree.hpd@data[,'height']) - mean(DFT2.beast.trees.MRCAs.rates[,'age at 818 singletons']),
                lwd = 1, col = 'black') +
  geom_point2(aes(subset = (node %in% 42),
                  x = x - c(2018.66575342466 - max(DFT2.beast.tree.hpd@data[,'height']) - mean(DFT2.beast.trees.MRCAs.rates[,'age at 818 singletons']))),
              shape = 16, size = 3, col = 'black') +
  geom_point2(aes(subset=(node %in% c(8,9,55,61,65))), shape = 25, size = 10, 
              fill = alpha('red', 0.8),
              col = alpha('red', 0.8))
dev.off()

## clean up environment
rm(list=ls())
