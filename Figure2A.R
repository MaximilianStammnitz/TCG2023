## Figure 2A
## 'The evolution of two transmissible cancers in Tasmanian devils'
## Stammnitz et al., 2022
## maxrupsta@gmail

library(treeio)
library(phytools)
library(ggtree)
library(scales)

## set input path(s)
setwd('/Users/ms37/Desktop/Labwork/DFT_evolution/doc/manuscripts/The Evolutionary History of Two Transmissible Cancers in Tasmanian Devils/Tables/v6/')


# DFT2-B subclonal tree #
#########################

## Import and plot DFT2-B RAxML tree
DFT2.raxml.tree <- read.tree("/Users/mstammnitz/Desktop/DFT_evolution/doc/manuscripts/The Evolutionary History of Two Transmissible Cancers in Tasmanian Devils/Tables/v6/Supplementary_data/RAxML_DFT2_subclones.raxml.support")
outgroup_node <- which(DFT2.raxml.tree$tip.label == "1545T2")
outgroup_edge <- which(DFT2.raxml.tree$edge[,2]==outgroup_node)
outgroup_edge_length <- DFT2.raxml.tree$edge.length[outgroup_edge]
rerooting_position <- outgroup_edge_length / 2
DFT2.raxml.tree <- phytools::reroot(DFT2.raxml.tree, outgroup_node, position = rerooting_position)

## Define splits
clade_b_split <- MRCA(DFT2.raxml.tree, c("1529T7", "1334T1"))

## SC colors
SC1.col <- 'grey1'
SC2.col <- 'grey80'
transp <- 0.3

### Plot
pdf('Figure2A_DFT2_subclone_tree.pdf', 
    width = 7, height = 6)
p <- ggtree(DFT2.raxml.tree, 
       layout = 'rectangular', 
       ladderize = T,
       right = F,
       lwd = 1) +
  geom_taxalink(taxa1 = "1509T1_SC1", 
                taxa2 = "1509T1_SC2", 
                color='black', curvature = -1.3, lwd = 1, lty = 3) +
  geom_balance(node=clade_b_split, fill=c(SC1.col, SC2.col), col = 'white', alpha = transp, extend = 0.03) +
  geom_tippoint(pch = 16, colour = "red", size = 6) + 
  geom_rootpoint(size = 3) +
  geom_rootedge(rootedge = 0.01, lwd = 1)
flip(p, 1, 2) %>% flip(16, 17)
dev.off()