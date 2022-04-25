## Figure 2D
## 'The evolution of two transmissible cancers in Tasmanian devils'
## Stammnitz et al., 2022
## maxrupsta@gmail

library(treeio)
library(phytools)
library(ggtree)
library(scales)

## set input path(s)
setwd('/Users/ms37/Desktop/Labwork/DFT_evolution/doc/manuscripts/The Evolutionary History of Two Transmissible Cancers in Tasmanian Devils/Tables/v6/')


# DFT1-C subclonal tree #
#########################

## Import and plot DFT2-B RAxML tree
DFT1.raxml.tree <- read.tree("/Users/mstammnitz/Desktop/DFT_evolution/doc/manuscripts/The Evolutionary History of Two Transmissible Cancers in Tasmanian Devils/Tables/v6/Supplementary_data/RAxML_DFT1_subclones.raxml.support")
outgroup_node <- which(DFT1.raxml.tree$tip.label == "1439T7")
outgroup_edge <- which(DFT1.raxml.tree$edge[,2]==outgroup_node)
outgroup_edge_length <- DFT1.raxml.tree$edge.length[outgroup_edge]
rerooting_position <- outgroup_edge_length / 2
DFT1.raxml.tree <- phytools::reroot(DFT1.raxml.tree, outgroup_node, position = rerooting_position)

## Define splits
clade_c_split <- MRCA(DFT1.raxml.tree, c("139T1_SC1", "139T1_SC2"))

## SC colors
SC1.col <- 'grey1'
SC2.col <- 'grey80'
transp <- 0.3

### Plot
pdf('Figure2D_DFT1_subclone_tree.pdf', 
    width = 7, height = 6)
p <- ggtree(DFT1.raxml.tree, 
            layout = 'rectangular', 
            ladderize = T,
            right = T,
            lwd = 1) +
  geom_taxalink(taxa1 = "139T1_SC1",
                taxa2 = "139T1_SC2", 
                color='black', curvature = -1.3,
                lwd = 1, lty = 3) +
  geom_balance(node = clade_c_split, fill = c(SC2.col, SC1.col), col = 'white', 
               alpha = transp, extend = 0.03) +
  geom_tippoint(pch = 16, colour = "cornflowerblue", size = 6) + 
  geom_rootpoint(size = 3) +
  geom_rootedge(rootedge = 0.01, lwd = 1)
rotate(p, 14) %>% rotate(15) %>% rotate(16)
dev.off()