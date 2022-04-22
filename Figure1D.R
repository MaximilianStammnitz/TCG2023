## Figure 1D
## 'The evolution of two transmissible cancers in Tasmanian devils'
## Stammnitz et al., 2022
## maxrupsta@gmail

## libraries
library(treeio)
library(ggtree)
library(ggsn)

## set input path(s)
setwd('/Tables')


# Dated DFT1 tree #
###################

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
clade_d <- MRCA(DFT1.beast.tree.hpd, '102T2')
clade_e <- MRCA(DFT1.beast.tree.hpd, '377T1')

## plot
pdf('Figure1D_DFT1_tree.pdf', width = 8, height = 12)
p <- ggtree(DFT1.beast.tree.hpd, 
            root.position = 2018.13424657534 - max(DFT1.beast.tree.hpd@data[,'height']), 
            layout = 'rectangular', 
            lwd = 1, 
            col = 'black',
            right = F,
            ladderize = T) + 
  scale_x_continuous(breaks = seq(f=1980, t=2020, by=10),
                     labels = seq(f=1980, t=2020, by=10),
                     limits = c(1980,2022)) +
  scale_y_continuous(expand = c(0.05,0)) +
  geom_tippoint(pch = 16, size = 3.5,
                col = 'cornflowerblue') +
  geom_range(range = 'height_0.95_HPD', center='height',
             color = 'cornflowerblue', 
             alpha = 0.8, size = 3.5) +
  geom_rootedge(rootedge = 2018.13424657534 - 
                  max(DFT1.beast.tree.hpd@data[,'height']) - 
                  mean(DFT1.beast.trees.MRCAs.rates[,'age at 818 singletons']),
                lwd = 1) +
  geom_point2(aes(subset = (node %in% 79),
                  x = x - c(2018.13424657534 - 
                              max(DFT1.beast.tree.hpd@data[,'height']) - 
                              mean(DFT1.beast.trees.MRCAs.rates[,'age at 818 singletons']))), 
              shape = 16, size = 3.5, col = 'black') +
  theme_classic(base_size = 20, base_family = 'Helvetica') +
  labs(x = 'Year', y = '') + 
  theme(axis.text = element_text(size = 25), 
        axis.text.y = element_blank(), 
        axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = 40, vjust = -1),
        axis.line.x = element_line(size = 1, color = 'black'),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = 'none',
        text = element_text(family='Helvetica')) +
  geom_cladelabel(clade_a1, label = 'A1', align = T, colour = 'cornflowerblue', 
                  offset = -39.5, offset.text = -19.3, fontsize = 10, family = 'Helvetica', barsize = 2) +
  geom_cladelabel(clade_a2, label = 'A2', align = T, colour = 'cornflowerblue', 
                  offset = -39.5, offset.text = -19.3, fontsize = 10, family = 'Helvetica', barsize = 2) +
  geom_cladelabel(clade_b, label = 'B', align = T, colour = 'cornflowerblue', 
                  offset = -39.5, offset.text = -19.3, fontsize = 10, family = 'Helvetica', barsize = 2) +
  geom_cladelabel(clade_c, label = 'C', align = T, colour = 'cornflowerblue', 
                  offset = -39.5, offset.text = -19.3, fontsize = 10, family = 'Helvetica', barsize = 2) +
  geom_cladelabel(clade_d, label = 'D', align = T, colour = 'cornflowerblue', 
                  offset = -48.2, offset.text = -19.3, fontsize = 10, family = 'Helvetica', barsize = 2) +
  geom_cladelabel(clade_e, label = 'E', align = T, colour = 'cornflowerblue', 
                  offset = -54, offset.text = -19.3, fontsize = 10, family = 'Helvetica', barsize = 2)
flip(p, clade_c, clade_c) %>% rotate(MRCA(DFT1.beast.tree.hpd, c('140T', '88T'))) %>% flip(clade_b, clade_d)
dev.off()

##
quantile(DFT1.beast.trees.MRCAs.rates[,'age at 0 singletons'], probs = c(0.05)) ## earliest inferred date of origin
mean(DFT1.beast.trees.MRCAs.rates[,'age at 818 singletons']) ## median inferred date of origin

## clean up environment
rm(list=ls())
