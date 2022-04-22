## Figure 1C
## 'The evolution of two transmissible cancers in Tasmanian devils'
## Stammnitz et al., 2022
## maxrupsta@gmail

## libraries
library(treeio)
library(ggtree)
library(scales)

## set input path(s)
setwd('/Users/ms37/Desktop/Labwork/DFT_evolution/doc/manuscripts/The Evolutionary History of Two Transmissible Cancers in Tasmanian Devils/Tables/v6/')


# Devil ancestry tree #
#######################

## input highest scoring RAxML-NG tree
Ancestry.tree <- read.tree('Supplementary_data/RAxML_Ancestry_tree.bestTree') ### provide on Github and/or Apollo
Ancestry.tree <- treeio::root(Ancestry.tree, node = 231)
class(Ancestry.tree$edge) <- 'integer'

## mark DFT1/DFT2 clades
clade_DFT1 <- MRCA(Ancestry.tree, c('155Tb', '366T1'))
clade_DFT2 <- MRCA(Ancestry.tree, c('1334T1', '203T3'))

## plot
pdf('Figure1C_devil_ancestry_tree.pdf', width = 14, height = 8)
ggtree(Ancestry.tree, 
       layout = 'daylight', 
       lwd = 1, 
       col = 'black',
       MAX_COUNT = 1) + 
  geom_hilight(node = clade_DFT1, fill = 'cornflowerblue', extend = 0.02) +
  geom_hilight(node = clade_DFT2, fill = 'red', extend = 0.02) +
  geom_tippoint(aes(subset=(node %in% grep('340T', Ancestry.tree$tip.label))), 
                shape = 16, size = 6, 
                fill = alpha('darkorange',0.8),
                col = alpha('darkorange',0.8))
dev.off()

## clean up environment
rm(list=ls())