## Figure 1E
## 'The evolution of two transmissible cancers in Tasmanian devils'
## Stammnitz et al., 2022
## maxrupsta@gmail

## libraries
library(treeio)
library(ggtree)
library(ggsn)

## set input path(s)
setwd('/Users/ms37/Desktop/Labwork/DFT_evolution/doc/manuscripts/The Evolutionary History of Two Transmissible Cancers in Tasmanian Devils/Tables/v6/')


# Dated DFT2 tree #
###################

## number of substitutions
DFT2.truncal <- 1335 ### number of genome-wide substitutions in the DFT2 trunk
DFT2.truncal.median <- DFT2.truncal - 818 ### median expected number of somatic substitutions in the DFT2 trunk

## import and process BEAST MCMC draws on DFT1 mutation rate and date of origin (generated with Tracer)
masked.mSarHar11.1 <- c('A' = 951740967, 'C' = 540077344, 'G' = 539833602, 'T' = 952098282)
DFT2.beast.trees.MRCAs.rates <- cbind(read.table('Supplementary_data/BEAST_DFT2_MRCAs.txt', header = T),
                                      'rates' = read.table('Supplementary_data/BEAST_DFT2_rates.txt', header = T)[,2])
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
DFT2.beast.tree.hpd <- read.beast('Supplementary_data/BEAST_DFT2.mcc')
DFT2.beast.tree.hpd@data[which.max(unlist(DFT2.beast.tree.hpd@data[,'height'])),"height_0.95_HPD"][[1]][[1]][2] <- 2018.66575342466 - 
  as.numeric(quantile(DFT2.beast.trees.MRCAs.rates[,'age at 0 singletons'], probs = c(0.05)))

## mark DFT2 clades
clade_a <- MRCA(DFT2.beast.tree.hpd, c('339T', '638T1'))
clade_b <- MRCA(DFT2.beast.tree.hpd, c('1529T2', '1545T2'))

## plot
pdf('Figure1E_DFT2_tree.pdf', width = 8, height = 12)
ggtree(DFT2.beast.tree.hpd, 
       root.position = 2018.66575342466 - max(DFT2.beast.tree.hpd@data[,'height']), 
       layout = 'rectangular', 
       lwd = 1, 
       col = 'black',
       ladderize = T) + 
  scale_x_continuous(breaks = seq(f=2008, t=2020, by=3),
                     labels = seq(f=2008, t=2020, by=3),
                     limits = c(2008,2020)) +
  scale_y_continuous(expand = c(0.05,0)) +
  geom_tippoint(pch = 16, size = 4.5,
                col = 'red') +
  geom_range(range = 'height_0.95_HPD', center='height',
             color = 'red', 
             alpha = 0.8, size = 4.5) +
  geom_rootedge(rootedge = 2018.66575342466 - max(DFT2.beast.tree.hpd@data[,'height']) - mean(DFT2.beast.trees.MRCAs.rates[,'age at 818 singletons']),
                lwd = 1) +
  geom_point2(aes(subset = (node %in% 42),
                  x = x - c(2018.66575342466 - max(DFT2.beast.tree.hpd@data[,'height']) - mean(DFT2.beast.trees.MRCAs.rates[,'age at 818 singletons']))),
              shape = 16, size = 4.5, col = 'black') +
  theme_classic(base_size = 20, base_family = "Helvetica") +
  labs(x = "Year", y = "") + 
  theme(axis.text = element_text(size = 25),
        axis.line.x = element_line(size = 1, color = 'black'),
        axis.text.y = element_blank(), 
        axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = 40, vjust = -1),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        text = element_text(family="Helvetica")) +
  geom_cladelabel(clade_a, label = "A", align = T, colour = "red", 
                  offset = -40, offset.text = -19.85, fontsize = 10, family = 'Helvetica', barsize = 2) +
  geom_cladelabel(clade_b, label = "B", align = T, colour = "red", 
                  offset = -40, offset.text = -19.85, fontsize = 10, family = 'Helvetica', barsize = 2)
dev.off()

##
quantile(DFT2.beast.trees.MRCAs.rates[,'age at 0 singletons'], probs = c(0.05)) ## earliest inferred date of origin
mean(DFT2.beast.trees.MRCAs.rates[,'age at 818 singletons']) ## median inferred date of origin

## clean up environment
rm(list=ls())