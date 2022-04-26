## Figure 3H
## 'The evolution of two transmissible cancers in Tasmanian devils'
## Stammnitz et al., 2022
## maxrupsta@gmail

library(readxl)
library(treeio)
library(ggtree)
library(ggsn)
library(stringr)
library(ggplot2)
library(ggmap)
library(rstudioapi)
library(lubridate)
library(Biostrings)
library(GenomicRanges)

## set input path(s)
setwd('/Users/ms37/Desktop/Labwork/DFT_evolution/doc/manuscripts/The Evolutionary History of Two Transmissible Cancers in Tasmanian Devils/Tables/v6/')


# DFT1 tree #
#############

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

pdf("Figure3H_DFT1_E_tree.pdf",
    width = 13, height = 11)

## hack in the edge colors and width
edge.sizes <- rep(1, DFT1.beast.tree.hpd@phylo$Nnode+length(DFT1.beast.tree.hpd@phylo$tip.label))

## highlight branch
edge.sizes[clade_e] <- 4
d <- data.frame(node=1:(DFT1.beast.tree.hpd@phylo$Nnode+length(DFT1.beast.tree.hpd@phylo$tip.label)), 
                size = edge.sizes)

p <- ggtree(DFT1.beast.tree.hpd, 
            root.position = 2018.13424657534 - max(DFT1.beast.tree.hpd@data[,'height']), 
            layout = 'rectangular', 
            col = 'black', 
            ladderize = T)  %<+% d + aes(lwd=I(size)) +
  scale_x_continuous(breaks = seq(f=1990, t=2020, by=10),
                     labels = seq(f=1990, t=2020, by=10),
                     limits = c(1986,2020)) +
  geom_rootedge(rootedge = 2018.13424657534 - max(DFT1.beast.tree.hpd@data[,'height']) - mean(DFT1.beast.trees.MRCAs.rates[,'age at 818 singletons']),
                lwd = 1, col = 'black') +
  geom_point2(aes(subset = (node %in% 79),
                  x = x - c(2018.13424657534 - max(DFT1.beast.tree.hpd@data[,'height']) - mean(DFT1.beast.trees.MRCAs.rates[,'age at 818 singletons']))),
              shape = 16, size = 3, col = 'black')  +
  geom_point2(aes(subset=(node %in% 46)), 
              shape = 24, size = 16, 
              fill = 'cornflowerblue',
              col = 'cornflowerblue') +
  geom_cladelabel(clade_e, label = "E", align = T, colour = "cornflowerblue", 
                  offset = -54, offset.text = -18.5, fontsize = 23, family = 'Helvetica') +
  labs(x = "Year", y = "") + 
  theme_classic(base_size = 45) +
  theme(axis.text = element_text(size = 45), 
        axis.text.y = element_blank(), 
        axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = 55, vjust = -1),
        axis.line.x = element_line(size = 1, color = 'black'),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        text = element_text(family="Helvetica"),
        plot.margin = unit(c(t = 0, r = 20, b = 20, l = 0), unit ='pt'))
flip(p, clade_c, clade_c) %>% rotate(MRCA(DFT1.beast.tree.hpd, c("140T", "88T"))) %>% flip(clade_b, clade_d)
dev.off()


# DFT1-E sampling map #
#######################

## clade E separation
samples <- as.matrix(read_xlsx('/Users/mstammnitz/Desktop/DFT_evolution/doc/manuscripts/The Evolutionary History of Two Transmissible Cancers in Tasmanian Devils/Tables/v6/Table-S2_v6.xlsx', sheet = 1))
dft1.samples <- samples[samples[,9] == 'DFT1',]
dft1.samples <- dft1.samples[which(is.na(dft1.samples[,9]) == F),]
dft1.samples <- dft1.samples[,c(8,5,6,7)]
dft1.samples[,1] <- str_split_fixed(dft1.samples[,1], ' ', 2)[,1]
colnames(dft1.samples) <- c('Tumour', 'Native Location', 'Latitude', 'Longitude')
summary <- dft1.samples[grep('377T1', dft1.samples[,"Tumour"]),,drop=F]
rownames(summary) <- summary[,1]
summary.sf <- summary[,3:4,drop=F]
class(summary.sf) <- 'numeric'
summary.sf <- as.data.frame(summary.sf)

## downlaod Tasmania map via ggmap package, requires a Google API key:
### video link ###
## register_google(key = '') fill in
tasmania.map.googlemaps_terrain <- get_googlemap(center = c(146.8087, -42.0409), 
                                                 zoom = 7,
                                                 size = c(640, 640),
                                                 scale = 2,
                                                 maptype = 'terrain',
                                                 color = 'bw',
                                                 style = 'style=feature:administrative.country|element:labels|visibility:off&style=feature:water|color:0xffffff|visibility:simplified')

## plot
png('Figure3H_DFT1_E_map.png', 
    width = 1100, height = 1100)
ggmap(tasmania.map.googlemaps_terrain,
      extent = 'panel',
      padding = 0)  +
  scale_x_continuous(limits = c(144.3, 148.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-43.7, -40.55), expand = c(0, 0)) +
  geom_point(mapping = aes(x = Longitude, y = Latitude, group = NULL), 
             data = summary.sf,
             size = 40,
             shape = 'triangle',
             show.legend = F,
             colour = 'cornflowerblue') +
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"), # get rid of legend panel bg
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())
dev.off()


# DFT1-E substitution rates #
#############################

## import data
counts <- as.matrix(read_xlsx('/Users/mstammnitz/Desktop/DFT_evolution/doc/manuscripts/The Evolutionary History of Two Transmissible Cancers in Tasmanian Devils/Tables/v6/Table-S3_v6.xlsx', sheet = 3))
colnames(counts) <- as.character(counts[2,])
counts <- counts[-c(1:2),]
dft1.counts <- counts[counts[,'LINEAGE'] == 'DFT1',]

## metadata
samples <- as.matrix(read_xlsx('/Users/mstammnitz/Desktop/DFT_evolution/doc/manuscripts/The Evolutionary History of Two Transmissible Cancers in Tasmanian Devils/Tables/v6/Table-S2_v6.xlsx', sheet = 1))
dft1.samples <- samples[samples[,9] == 'DFT1',]
dft1.samples <- dft1.samples[which(is.na(dft1.samples[,9]) == F),]
dft1.samples <- dft1.samples[,c(8,4,13)]
dft1.samples[,1] <- str_split_fixed(dft1.samples[,1], ' ', 2)[,1]
dft1.samples[,2] <- gsub('[*]', '', dft1.samples[,2])
colnames(dft1.samples) <- c('Tumour', 'Collection Date', 'Purity')

## combine
DFT1.SNVs.counts <- dft1.samples
DFT1.SNVs.counts <- cbind(DFT1.SNVs.counts[match(dft1.counts[,'TUMOUR ID'], DFT1.SNVs.counts[,1]),],
                          dft1.counts[,"WGD-NORMALISED SUBSTITUTIONS"])
colnames(DFT1.SNVs.counts)[4] <- 'SNVs'
rownames(DFT1.SNVs.counts) <- DFT1.SNVs.counts[,1]
DFT1.SNVs.counts <- DFT1.SNVs.counts[,-1]
DFT1.SNVs.counts[,'Collection Date'] <- decimal_date(as.Date(as.character(unlist(DFT1.SNVs.counts[,'Collection Date'])), format="%d.%m.%Y"))
DFT1.SNVs.counts <- DFT1.SNVs.counts[!is.na(DFT1.SNVs.counts[,'Collection Date']),]
DFT1.SNVs.counts <- as.data.frame(DFT1.SNVs.counts)
DFT1.SNVs.counts[,'Collection Date'] <- as.numeric(as.character(DFT1.SNVs.counts[,'Collection Date']))
DFT1.SNVs.counts[,'Purity'] <- as.numeric(as.character(DFT1.SNVs.counts[,'Purity']))
DFT1.SNVs.counts[,'SNVs'] <- as.numeric(as.character(DFT1.SNVs.counts[,'SNVs']))

pdf("Figure3H_substitution_rates.pdf", 
    height = 9, width = 6)
ggplot(DFT1.SNVs.counts, aes(x = `Collection Date`, y = `SNVs`)) + 
  scale_x_continuous(breaks = seq(from = 2000, to = 2020, by = 10), limits = c(1985, 2035)) + 
  scale_y_continuous(breaks = seq(from = 0, to = 30000, by = 6000), limits = c(0, 35000)) + 
  coord_cartesian(xlim = c(2000, 2020), ylim = c(0, 35000), expand = F) +
  geom_point(data = DFT1.SNVs.counts[-grep('377T1', rownames(DFT1.SNVs.counts)),], color = "cornflowerblue", 
             size = 2, shape = 'circle') + 
  geom_point(data = DFT1.SNVs.counts[grep('377T1', rownames(DFT1.SNVs.counts)),], color = "cornflowerblue",
             size = 18, shape = 'triangle') + 
  geom_smooth(data = DFT1.SNVs.counts[-grep('377T1', rownames(DFT1.SNVs.counts)),],
              method = 'lm', color = "cornflowerblue", fullrange = T) +
  labs(y = "Substitutions", x = "Year") +
  theme_classic(base_size = 30) +
  theme(axis.text = element_text(size = 30),
        axis.line.x = element_line(size = 1, color = 'black'),
        axis.line.y = element_line(size = 1, color = 'black'),
        axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = 50, vjust = -1),
        axis.title.y = element_text(family = 'Helvetica', colour = 'black', size = 50, vjust = 3),
        legend.position = "none",
        text = element_text(family="Helvetica"),
        plot.margin = unit(c(2, 2, 2, 2),"cm"))
dev.off()


# DFT1-E indel rates #
######################

## import data
counts <- as.matrix(read_xlsx('/Users/mstammnitz/Desktop/DFT_evolution/doc/manuscripts/The Evolutionary History of Two Transmissible Cancers in Tasmanian Devils/Tables/v6/Table-S3_v6.xlsx', sheet = 4))
colnames(counts) <- as.character(counts[2,])
counts <- counts[-c(1:2),]
dft1.counts <- counts[counts[,'LINEAGE'] == 'DFT1',]

## metadata
samples <- as.matrix(read_xlsx('/Users/mstammnitz/Desktop/DFT_evolution/doc/manuscripts/The Evolutionary History of Two Transmissible Cancers in Tasmanian Devils/Tables/v6/Table-S2_v6.xlsx', sheet = 1))
dft1.samples <- samples[samples[,9] == 'DFT1',]
dft1.samples <- dft1.samples[which(is.na(dft1.samples[,9]) == F),]
dft1.samples <- dft1.samples[,c(8,4,13)]
dft1.samples[,1] <- str_split_fixed(dft1.samples[,1], ' ', 2)[,1]
dft1.samples[,2] <- gsub('[*]', '', dft1.samples[,2])
colnames(dft1.samples) <- c('Tumour', 'Collection Date', 'Purity')

## combine
DFT1.Indels.counts <- dft1.samples
DFT1.Indels.counts <- cbind(DFT1.Indels.counts[match(dft1.counts[,'TUMOUR ID'], DFT1.Indels.counts[,1]),],
                          dft1.counts[,"WGD-NORMALISED INDELS"])
colnames(DFT1.Indels.counts)[4] <- 'Indels'
rownames(DFT1.Indels.counts) <- DFT1.Indels.counts[,1]
DFT1.Indels.counts <- DFT1.Indels.counts[,-1]
DFT1.Indels.counts[,'Collection Date'] <- decimal_date(as.Date(as.character(unlist(DFT1.Indels.counts[,'Collection Date'])), format="%d.%m.%Y"))
DFT1.Indels.counts <- DFT1.Indels.counts[!is.na(DFT1.Indels.counts[,'Collection Date']),]
DFT1.Indels.counts <- as.data.frame(DFT1.Indels.counts)
DFT1.Indels.counts[,'Collection Date'] <- as.numeric(as.character(DFT1.Indels.counts[,'Collection Date']))
DFT1.Indels.counts[,'Purity'] <- as.numeric(as.character(DFT1.Indels.counts[,'Purity']))
DFT1.Indels.counts[,'Indels'] <- as.numeric(as.character(DFT1.Indels.counts[,'Indels']))

pdf("Figure3H_indel_rates.pdf", 
    height = 9, width = 6)
ggplot(DFT1.Indels.counts, aes(x = `Collection Date`, y = `Indels`)) + 
  scale_x_continuous(breaks = seq(from = 2000, to = 2020, by = 10), limits = c(1985, 2035)) + 
  scale_y_continuous(breaks = seq(from = 0, to = 10000, by = 2000), limits = c(0, 11000)) + 
  coord_cartesian(xlim = c(2000, 2020), ylim = c(0, 11000), expand = F) +
  geom_point(data = DFT1.Indels.counts[-grep('377T1', rownames(DFT1.Indels.counts)),], color = "cornflowerblue", 
             size = 2, shape = 'circle') + 
  geom_point(data = DFT1.Indels.counts[grep('377T1', rownames(DFT1.Indels.counts)),], color = "cornflowerblue", 
             size = 18, shape = 'triangle') + 
  geom_smooth(data = DFT1.Indels.counts[-grep('377T1', rownames(DFT1.Indels.counts)),],
              method = 'lm', color = "cornflowerblue", fullrange = T) +
  labs(y = "Indels", x = "Year") +
  theme_classic(base_size = 30) +
  theme(axis.text = element_text(size = 30),
        axis.line.x = element_line(size = 1, color = 'black'),
        axis.line.y = element_line(size = 1, color = 'black'),
        axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = 50, vjust = -1),
        axis.title.y = element_text(family = 'Helvetica', colour = 'black', size = 50, vjust = 3),
        legend.position = "none",
        text = element_text(family="Helvetica"),
        plot.margin = unit(c(2, 2, 2, 2),"cm"))
dev.off()


# DFT1-E spectra #
##################

## import data
load("/Users/mstammnitz/Desktop/DFT_evolution/doc/manuscripts/The Evolutionary History of Two Transmissible Cancers in Tasmanian Devils/Tables/v6/Supplementary_data/DFT1_DFT2_340T_somatic_variants.Rdata")
reference <- readDNAStringSet('/Users/mstammnitz/Desktop/DFT_evolution/doc/manuscripts/The Evolutionary History of Two Transmissible Cancers in Tasmanian Devils/Tables/v6/Supplementary_data/sarcophilus_harrisii_toplevel.fa.gz')
reference_trinucleotides <- read.table('/Users/mstammnitz/Desktop/DFT_evolution/doc/manuscripts/The Evolutionary History of Two Transmissible Cancers in Tasmanian Devils/Tables/v6/Supplementary_data/mSarhar1.11_trinucleotides.txt', header = T)

## functions for substitution processing (GITHUB repository: https://github.com/MaximilianStammnitz/SubstitutionSafari)
substitution.spectrum <- function(x, normalised){
  
  # a. Isolate triplets form Platypus VCF
  x[,'TRIPLET'] <- as.character(subseq(x = reference[as.character(x[,'CHROM'])], 
                                       start = as.numeric(x[,'POS']) - 1, 
                                       end = as.numeric(x[,'POS']) + 1))
  
  # b. Isolate base changes
  context <- x[,'TRIPLET']
  context.changes <- matrix(NA, nrow = length(context), ncol = 3)
  colnames(context.changes) <- c("REF", "ALT-5'", "ALT-3'")
  context.changes[,1] <- context
  context.changes[,2:3] <- str_split_fixed(context.changes[,1],"",3)[,c(1,3)]
  context.changes[,2] <- paste(context.changes[,2],as.character(x[,'ALT']),context.changes[,3],sep="")
  context.changes <- context.changes[,c(1,2)]
  colnames(context.changes) <- c("REF", "ALT")
  
  # Pyrimidine-context substitutions
  equivalent.mut <- matrix(c('ACA>AAA', 'TGT>TTT', # AC - AA, or TG - TT ###  C>A or G>T
                             'ACC>AAC', 'GGT>GTT',
                             'ACG>AAG', 'CGT>CTT',
                             'ACT>AAT', 'AGT>ATT',
                             'CCA>CAA', 'TGG>TTG', # CC - CA, or GG - GT
                             'CCC>CAC', 'GGG>GTG',
                             'CCG>CAG', 'CGG>CTG',
                             'CCT>CAT', 'AGG>ATG',
                             'GCA>GAA', 'TGC>TTC', # GC - GA, or CG - CT
                             'GCC>GAC', 'GGC>GTC',
                             'GCG>GAG', 'CGC>CTC',
                             'GCT>GAT', 'AGC>ATC',
                             'TCA>TAA', 'TGA>TTA', # TC - TA, or AG - AT
                             'TCC>TAC', 'GGA>GTA',
                             'TCG>TAG', 'CGA>CTA',
                             'TCT>TAT', 'AGA>ATA',
                             'ACA>AGA', 'TGT>TCT', # AC - AG, or TG - TC ###  C>G or G>C
                             'ACC>AGC', 'GGT>GCT',
                             'ACG>AGG', 'CGT>CCT',
                             'ACT>AGT', 'AGT>ACT',
                             'CCA>CGA', 'TGG>TCG', # CC - CG, or GG - GC
                             'CCC>CGC', 'GGG>GCG',
                             'CCG>CGG', 'CGG>CCG',
                             'CCT>CGT', 'AGG>ACG',
                             'GCA>GGA', 'TGC>TCC', # GC - GG, or CG - CC
                             'GCC>GGC', 'GGC>GCC',
                             'GCG>GGG', 'CGC>CCC',
                             'GCT>GGT', 'AGC>ACC',
                             'TCA>TGA', 'TGA>TCA', # TC - TG, or AG - AC
                             'TCC>TGC', 'GGA>GCA',
                             'TCG>TGG', 'CGA>CCA',
                             'TCT>TGT', 'AGA>ACA',
                             'ACA>ATA', 'TGT>TAT', # AC - AT, or TG - TA ###  C>T or G>A
                             'ACC>ATC', 'GGT>GAT',
                             'ACG>ATG', 'CGT>CAT',
                             'ACT>ATT', 'AGT>AAT',
                             'CCA>CTA', 'TGG>TAG', # CC - CT, or GG - GA
                             'CCC>CTC', 'GGG>GAG',
                             'CCG>CTG', 'CGG>CAG',
                             'CCT>CTT', 'AGG>AAG',
                             'GCA>GTA', 'TGC>TAC', # GC - GT, or CG - CA
                             'GCC>GTC', 'GGC>GAC',
                             'GCG>GTG', 'CGC>CAC',
                             'GCT>GTT', 'AGC>AAC',
                             'TCA>TTA', 'TGA>TAA', # TC - TT, or AG - AA
                             'TCC>TTC', 'GGA>GAA',
                             'TCG>TTG', 'CGA>CAA',
                             'TCT>TTT', 'AGA>AAA',
                             'ATA>AAA', 'TAT>TTT', # AT - AA, or TA - TT ###  T>A or A>T
                             'ATC>AAC', 'GAT>GTT',
                             'ATG>AAG', 'CAT>CTT',
                             'ATT>AAT', 'AAT>ATT',
                             'CTA>CAA', 'TAG>TTG', # CT - CA, or GA - GT
                             'CTC>CAC', 'GAG>GTG',
                             'CTG>CAG', 'CAG>CTG',
                             'CTT>CAT', 'AAG>ATG',
                             'GTA>GAA', 'TAC>TTC', # GT - GA, or CA - CT
                             'GTC>GAC', 'GAC>GTC',
                             'GTG>GAG', 'CAC>CTC',
                             'GTT>GAT', 'AAC>ATC',
                             'TTA>TAA', 'TAA>TTA', # TT - TA, or AA - AT
                             'TTC>TAC', 'GAA>GTA',
                             'TTG>TAG', 'CAA>CTA',
                             'TTT>TAT', 'AAA>ATA',
                             'ATA>ACA', 'TAT>TGT', # AT - AC, or TA - TG ### T>C or A>G
                             'ATC>ACC', 'GAT>GGT',
                             'ATG>ACG', 'CAT>CGT',
                             'ATT>ACT', 'AAT>AGT',
                             'CTA>CCA', 'TAG>TGG', # CT - CC, or GA - GG
                             'CTC>CCC', 'GAG>GGG',
                             'CTG>CCG', 'CAG>CGG',
                             'CTT>CCT', 'AAG>AGG',
                             'GTA>GCA', 'TAC>TGC', # GT - GC, or CA - CG
                             'GTC>GCC', 'GAC>GGC',
                             'GTG>GCG', 'CAC>CGC',
                             'GTT>GCT', 'AAC>AGC',
                             'TTA>TCA', 'TAA>TGA', # TT - TC, or AA - AG
                             'TTC>TCC', 'GAA>GGA',
                             'TTG>TCG', 'CAA>CGA',
                             'TTT>TCT', 'AAA>AGA',
                             'ATA>AGA', 'TAT>TCT', # AT - AG, or TA - TC ### T>G or A>C
                             'ATC>AGC', 'GAT>GCT',
                             'ATG>AGG', 'CAT>CCT',
                             'ATT>AGT', 'AAT>ACT',
                             'CTA>CGA', 'TAG>TCG', # CT - CG, or GA - GC
                             'CTC>CGC', 'GAG>GCG',
                             'CTG>CGG', 'CAG>CCG',
                             'CTT>CGT', 'AAG>ACG',
                             'GTA>GGA', 'TAC>TCC', # GT - GG, or CA - CC
                             'GTC>GGC', 'GAC>GCC',
                             'GTG>GGG', 'CAC>CCC',
                             'GTT>GGT', 'AAC>ACC',
                             'TTA>TGA', 'TAA>TCA', # TT - TG, or AA - AC
                             'TTC>TGC', 'GAA>GCA',
                             'TTG>TGG', 'CAA>CCA',
                             'TTT>TGT', 'AAA>ACA'),
                           nrow=96, ncol=2, byrow=TRUE)
  
  # c. Define variants from input table and convert them to pyrimidine context
  triplett_snvs <- paste0(context.changes[,1], '>', context.changes[,2])
  ind.convert <- match(triplett_snvs,equivalent.mut[,2])
  triplett_snvs[!is.na(ind.convert)] <- equivalent.mut[ind.convert[!is.na(ind.convert)],1]
  
  # d. Count triplets
  counts <- table(triplett_snvs)
  consensus.mut <- equivalent.mut[,1]
  if(length(counts)<96){
    add.names <- consensus.mut[which(is.na(match(consensus.mut,names(counts)))==T)]
    names.takeover <- names(counts)
    counts <- c(counts,rep(0,length(add.names)))
    names(counts) <- c(names.takeover,add.names)
  }
  counts <- counts[match(consensus.mut,names(counts))]
  
  # e. normalisation to genome triplet counts
  if(normalised == T){
    
    reference_trinucleotides.div <- c(rep(as.numeric(reference_trinucleotides[,5])[1:16],3),
                                      rep(as.numeric(reference_trinucleotides[,5])[17:32],3))
    counts.norm <- c(counts/reference_trinucleotides.div)/sum(c(counts/reference_trinucleotides.div))
    names(counts) <- names(counts.norm) <- consensus.mut
    out <- list("counts" = counts, "counts.normalised" = counts.norm)
    
  } else if (normalised == F){
    
    names(counts) <- names(counts.norm) <- consensus.mut
    out <- list("counts" = counts, "counts.normalised" = counts.norm)
    
  }
  
  # f. Output
  return(out)
}
plot.substition.spectrum <- function(x, title, peak.colour){
  
  ## setup plot background colours and blocks
  mar.default <- c(2,4,2,2) + 0.1
  par(mar = mar.default + c(3, 11, 3, -2))
  colors = c("grey40", "grey60", "grey40", "grey60", "grey40", "grey60")
  mut <- c('C > A','C > G','C > T','T > A','T > C','T > G')
  y.top = 0.35
  borders = c(0, 115.4*1/6, 115.4*2/6, 115.4*3/6, 115.4*4/6, 115.4*5/6, 115.4)
  alphas = rep(0.2,6)
  plot(1, type="n", ylim=c(0,y.top), xlim=c(0, 115.4), xlab="", ylab="", axes=F,
       main = title,
       cex.main = 4)
  for (i in 1:6) {
    rect(xleft=borders[i], xright=borders[i+1], ybottom=0, ytop=y.top, col=alpha(colors[i], alphas[i]), border="white")
    rect(xleft=borders[i], xright=borders[i+1], ybottom=y.top-c(y.top*0.12), ytop=y.top, col=colors[i], border="white")
    text(x=(borders[i]+borders[i+1])/2, y=y.top-c(y.top*0.12)/2, labels=mut[i], cex=5, col="white")
  }
  
  ## add spectrum
  out <- barplot(x/sum(x),
                 col = rep(peak.colour,96),
                 border = NA,
                 ylab = "",
                 yaxt = 'n',
                 ylim = c(0, 0.35),
                 add = T,
                 names.arg = NA)
  
  ## add title
  title(ylab="Substitutions [%]", line = 9, cex.lab = 6)
  
  # add axes
  axis(side = 1, 
       las = 2,
       at = out[,1],
       labels = c("A[C>A]A", "A[C>A]C", "A[C>A]G", "A[C>A]T", "C[C>A]A",
                  "C[C>A]C", "C[C>A]G", "C[C>A]T", "G[C>A]A", "G[C>A]C",
                  "G[C>A]G", "G[C>A]T", "T[C>A]A", "T[C>A]C", "T[C>A]G",
                  "T[C>A]T", "A[C>G]A", "A[C>G]C", "A[C>G]G", "A[C>G]T",
                  "C[C>G]A", "C[C>G]C", "C[C>G]G", "C[C>G]T", "G[C>G]A",
                  "G[C>G]C", "G[C>G]G", "G[C>G]T", "T[C>G]A", "T[C>G]C",
                  "T[C>G]G", "T[C>G]T", "A[C>T]A", "A[C>T]C", "A[C>T]G",
                  "A[C>T]T", "C[C>T]A", "C[C>T]C", "C[C>T]G", "C[C>T]T",
                  "G[C>T]A", "G[C>T]C", "G[C>T]G", "G[C>T]T", "T[C>T]A",
                  "T[C>T]C", "T[C>T]G", "T[C>T]T", "A[T>A]A", "A[T>A]C",
                  "A[T>A]G", "A[T>A]T", "C[T>A]A", "C[T>A]C", "C[T>A]G",
                  "C[T>A]T", "G[T>A]A", "G[T>A]C", "G[T>A]G", "G[T>A]T",
                  "T[T>A]A", "T[T>A]C", "T[T>A]G", "T[T>A]T", "A[T>C]A",
                  "A[T>C]C", "A[T>C]G", "A[T>C]T", "C[T>C]A", "C[T>C]C",
                  "C[T>C]G", "C[T>C]T", "G[T>C]A", "G[T>C]C", "G[T>C]G",
                  "G[T>C]T", "T[T>C]A", "T[T>C]C", "T[T>C]G", "T[T>C]T",
                  "A[T>G]A", "A[T>G]C", "A[T>G]G", "A[T>G]T", "C[T>G]A",
                  "C[T>G]C", "C[T>G]G", "C[T>G]T", "G[T>G]A", "G[T>G]C",
                  "G[T>G]G", "G[T>G]T", "T[T>G]A", "T[T>G]C", "T[T>G]G",
                  "T[T>G]T"),
       tick = F, line = -1.6, font = 1, family = 'mono', cex.axis = 1.2)
  axis(side = 2, las = 2, cex.axis = 4, col = 'black', col.axis = 'black',
       at = c(0, 0.1, 0.2, 0.3), labels = c('0', '10', '20', '30'), lwd = 5, pos=-2, hadj = 1.4)
  
}

## functions for indel processing (GITHUB repository: https://github.com/MaximilianStammnitz/Indelwald)
indel.spectrum <- function(x){
  
  ## 1. split VCF indels into:
  # (i) 1-bp deletions
  dels.1bp <- x[nchar(as.character(x[,'ALT'])) == nchar(as.character(x[,'REF']))-1,,drop=F]
  
  # (ii) 1-bp insertions
  ins.1bp <- x[nchar(as.character(x[,'REF'])) == nchar(as.character(x[,'ALT']))-1,,drop=F]
  
  # (iii) >1 bp deletions
  
  ## 2
  dels.2bp <- x[nchar(as.character(x[,'ALT'])) == nchar(as.character(x[,'REF']))-2,,drop=F]
  
  ## 3
  dels.3bp <- x[nchar(as.character(x[,'ALT'])) == nchar(as.character(x[,'REF']))-3,,drop=F]
  
  ## 4
  dels.4bp <- x[nchar(as.character(x[,'ALT'])) == nchar(as.character(x[,'REF']))-4,,drop=F]
  
  ## 5+
  dels.5bp <- x[nchar(as.character(x[,'ALT'])) <= nchar(as.character(x[,'REF']))-5,,drop=F]
  
  # (iv) >1 bp insertions
  
  ## 2
  ins.2bp <- x[nchar(as.character(x[,'REF'])) == nchar(as.character(x[,'ALT']))-2,,drop=F]
  
  ## 3
  ins.3bp <- x[nchar(as.character(x[,'REF'])) == nchar(as.character(x[,'ALT']))-3,,drop=F]
  
  ## 4
  ins.4bp <- x[nchar(as.character(x[,'REF'])) == nchar(as.character(x[,'ALT']))-4,,drop=F]
  
  ## 5+
  ins.5bp <- x[nchar(as.character(x[,'REF'])) <= nchar(as.character(x[,'ALT']))-5,,drop=F]
  
  ## 2. classify 1 bp events into:
  
  # (i) 1 bp deletions at homopolymers (length 1 == "no neighbouring homopolymers")
  if(nrow(dels.1bp) > 0){
    
    ## extract 10 bp upstream/downstream sequence context from reference
    dels.1bp.context <- as.character(subseq(x = reference[as.character(dels.1bp[,'CHROM'])], 
                                            start = as.numeric(dels.1bp[,'POS']) - 9, 
                                            end = as.numeric(dels.1bp[,'POS']) + 11))
    dels.1bp.context.middle <- paste0('[', str_split_fixed(dels.1bp.context, '', 21)[,11,drop=F], ']')
    dels.1bp.context.start <- str_split_fixed(dels.1bp.context, '', 11)[,1:10,drop=F]
    dels.1bp.context.start <- paste(dels.1bp.context.start[,1], dels.1bp.context.start[,2], dels.1bp.context.start[,3],
                                    dels.1bp.context.start[,4], dels.1bp.context.start[,5], dels.1bp.context.start[,6],
                                    dels.1bp.context.start[,7], dels.1bp.context.start[,8], dels.1bp.context.start[,9],
                                    dels.1bp.context.start[,10], sep = '')
    dels.1bp.context.end <- str_split_fixed(dels.1bp.context, '', 12)[,12,drop=F]
    dels.1bp[,'TRIPLET'] <- paste(dels.1bp.context.start, dels.1bp.context.middle, dels.1bp.context.end, sep = '')
    colnames(dels.1bp)[5] <- 'CONTEXT FW' 
    
    ## need the central base to be pyrimidine-centred, i.e. C or T, hence also run reverseComplements on full contexts
    dels.1bp.context.rc <- as.character(reverseComplement(subseq(x = reference[as.character(dels.1bp[,'CHROM'])], 
                                                                 start = as.numeric(dels.1bp[,'POS']) - 9, 
                                                                 end = as.numeric(dels.1bp[,'POS']) + 11)))
    dels.1bp.context.rc.middle <- paste0('[', str_split_fixed(dels.1bp.context.rc, '', 21)[,11,drop=F], ']')
    dels.1bp.context.rc.start <- str_split_fixed(dels.1bp.context.rc, '', 11)[,1:10,drop=F]
    dels.1bp.context.rc.start <- paste(dels.1bp.context.rc.start[,1], dels.1bp.context.rc.start[,2], dels.1bp.context.rc.start[,3],
                                       dels.1bp.context.rc.start[,4], dels.1bp.context.rc.start[,5], dels.1bp.context.rc.start[,6],
                                       dels.1bp.context.rc.start[,7], dels.1bp.context.rc.start[,8], dels.1bp.context.rc.start[,9],
                                       dels.1bp.context.rc.start[,10], sep = '')
    dels.1bp.context.rc.end <- str_split_fixed(dels.1bp.context.rc, '', 12)[,12,drop=F]
    dels.1bp <- cbind(dels.1bp[,1:5,drop=F], paste(dels.1bp.context.rc.start, dels.1bp.context.rc.middle, dels.1bp.context.rc.end, sep = ''))
    colnames(dels.1bp)[6] <- 'CONTEXT RC'
    
    ## summarise 1 bp deletions in matrix format
    dels.1bp.summary <- matrix(0, ncol = 2, nrow = 6)
    colnames(dels.1bp.summary) <- c('C', 'T') ## pyrimidine-centred deleted base
    rownames(dels.1bp.summary) <- c('0 bp', '1 bp', '2 bp', '3 bp', '4 bp', '5+ bp') ## contextual homopolymer-length
    fw <- str_split_fixed(str_split_fixed(dels.1bp[,'CONTEXT FW'], '\\[', 2)[,2,drop=F], '\\]', 2)[,1,drop=F]
    rc <- str_split_fixed(str_split_fixed(dels.1bp[,'CONTEXT RC'], '\\[', 2)[,2,drop=F], '\\]', 2)[,1,drop=F]
    fw.if <- fw == 'C' | fw == 'T'
    for (i in 1:nrow(dels.1bp)){
      
      if(fw.if[i] == T){
        
        ## look at forward context
        upstream.tmp <- str_split_fixed(str_split_fixed(as.character(dels.1bp[i,'CONTEXT FW']), '\\[', 2)[,1], '', 10)
        upstream.tmp <- which(upstream.tmp == fw[i])
        
        ### check all 6 categories for upstream bases
        if(any(upstream.tmp %in% 10)){
          
          if(any(upstream.tmp %in% 9)){
            
            if(any(upstream.tmp %in% 8)){
              
              if(any(upstream.tmp %in% 7)){
                
                if(any(upstream.tmp %in% 6)){
                  
                  upstream.tmp <- 6
                  
                }else{
                  upstream.tmp <- 5
                }
                
              }else{
                upstream.tmp <- 4
              }
              
            }else{
              upstream.tmp <- 3
            }
            
          }else{
            upstream.tmp <- 2
          }
          
        }else{
          upstream.tmp <- 1
        }
        
        downstream.tmp <- str_split_fixed(str_split_fixed(as.character(dels.1bp[i,'CONTEXT FW']), '\\]', 2)[,2], '', 10)
        downstream.tmp <- which(downstream.tmp == fw[i])
        
        ### check all 6 categories for downstream bases
        if(any(downstream.tmp %in% 1)){
          
          if(any(downstream.tmp %in% 2)){
            
            if(any(downstream.tmp %in% 3)){
              
              if(any(downstream.tmp %in% 4)){
                
                if(any(downstream.tmp %in% 5)){
                  
                  downstream.tmp <- 6
                  
                }else{
                  downstream.tmp <- 5
                }
                
              }else{
                downstream.tmp <- 4
              }
              
            }else{
              downstream.tmp <- 3
            }
            
          }else{
            downstream.tmp <- 2
          }
          
        }else{
          downstream.tmp <- 1
        }
        
        ## summarise, which homopolymer (upstream vs. downstream context) is longer
        dels.1bp.summary[max(c(upstream.tmp, downstream.tmp)),fw[i]] <- dels.1bp.summary[max(c(upstream.tmp, downstream.tmp)),fw[i]] + 1
        
      } else {
        
        ## look at reverse complement context
        upstream.tmp <- str_split_fixed(str_split_fixed(as.character(dels.1bp[i,'CONTEXT RC']), '\\[', 2)[,1], '', 10)
        upstream.tmp <- which(upstream.tmp == rc[i])
        
        ### check all 6 categories for upstream bases
        if(any(upstream.tmp %in% 10)){
          
          if(any(upstream.tmp %in% 9)){
            
            if(any(upstream.tmp %in% 8)){
              
              if(any(upstream.tmp %in% 7)){
                
                if(any(upstream.tmp %in% 6)){
                  
                  upstream.tmp <- 6
                  
                }else{
                  upstream.tmp <- 5
                }
                
              }else{
                upstream.tmp <- 4
              }
              
            }else{
              upstream.tmp <- 3
            }
            
          }else{
            upstream.tmp <- 2
          }
          
        }else{
          upstream.tmp <- 1
        }
        
        downstream.tmp <- str_split_fixed(str_split_fixed(as.character(dels.1bp[i,'CONTEXT RC']), '\\]', 2)[,2], '', 10)
        downstream.tmp <- which(downstream.tmp == rc[i])
        
        ### check all 6 categories for downstream bases
        if(any(downstream.tmp %in% 1)){
          
          if(any(downstream.tmp %in% 2)){
            
            if(any(downstream.tmp %in% 3)){
              
              if(any(downstream.tmp %in% 4)){
                
                if(any(downstream.tmp %in% 5)){
                  
                  downstream.tmp <- 6
                  
                }else{
                  downstream.tmp <- 5
                }
                
              }else{
                downstream.tmp <- 4
              }
              
            }else{
              downstream.tmp <- 3
            }
            
          }else{
            downstream.tmp <- 2
          }
          
        }else{
          downstream.tmp <- 1
        }
        
        ## summarise, which homopolymer (upstream vs. downstream context) is longer
        dels.1bp.summary[max(c(upstream.tmp, downstream.tmp)),rc[i]] <- dels.1bp.summary[max(c(upstream.tmp, downstream.tmp)),rc[i]] + 1
      }
      
    }
    
  }else{
    
    ## summarise 1 bp deletions in matrix format
    dels.1bp.summary <- matrix(0, ncol = 2, nrow = 6)
    colnames(dels.1bp.summary) <- c('C', 'T') ## pyrimidine-centred deleted base
    rownames(dels.1bp.summary) <- c('0 bp', '1 bp', '2 bp', '3 bp', '4 bp', '5+ bp') ## contextual homopolymer-length
    
  }
  
  # (ii) 1 bp insertions at homopolymers (length 0 == "no neighbouring homopolymers")
  if(nrow(ins.1bp) > 0){
    
    ## extract 10 bp upstream/downstream sequence context from reference
    ins.1bp.context <- as.character(subseq(x = reference[as.character(ins.1bp[,'CHROM'])], 
                                           start = as.numeric(ins.1bp[,'POS']) - 9, 
                                           end = as.numeric(ins.1bp[,'POS']) + 10))
    ## insertion after base pos. 10
    ins.1bp.context.middle <- paste0('[', str_split_fixed(as.character(ins.1bp[,'ALT']), '', 2)[,2,drop=F], ']')
    ins.1bp.context.start <- str_split_fixed(ins.1bp.context, '', 11)[,1:10,drop=F]
    ins.1bp.context.start <- paste(ins.1bp.context.start[,1], ins.1bp.context.start[,2], ins.1bp.context.start[,3],
                                   ins.1bp.context.start[,4], ins.1bp.context.start[,5], ins.1bp.context.start[,6],
                                   ins.1bp.context.start[,7], ins.1bp.context.start[,8], ins.1bp.context.start[,9],
                                   ins.1bp.context.start[,10], sep = '')
    ins.1bp.context.end <- str_split_fixed(ins.1bp.context, '', 11)[,11,drop=F]
    ins.1bp[,'TRIPLET'] <- paste(ins.1bp.context.start, ins.1bp.context.middle, ins.1bp.context.end, sep = '')
    colnames(ins.1bp)[5] <- 'CONTEXT FW'
    
    ## need the central base to be pyrimidine-centred, i.e. C or T, hence also run reverseComplements on full contexts
    ins.1bp.context.rc <- as.character(reverseComplement(subseq(x = reference[as.character(ins.1bp[,'CHROM'])], 
                                                                start = as.numeric(ins.1bp[,'POS']) - 9, 
                                                                end = as.numeric(ins.1bp[,'POS']) + 10)))
    ins.1bp.context.rc.middle <- paste0('[', as.character(reverseComplement(DNAStringSet(str_split_fixed(as.character(ins.1bp[,'ALT']), '', 2)[,2,drop=F]))), ']')
    ins.1bp.context.rc.start <- str_split_fixed(ins.1bp.context.rc, '', 11)[,1:10,drop=F]
    ins.1bp.context.rc.start <- paste(ins.1bp.context.rc.start[,1], ins.1bp.context.rc.start[,2], ins.1bp.context.rc.start[,3],
                                      ins.1bp.context.rc.start[,4], ins.1bp.context.rc.start[,5], ins.1bp.context.rc.start[,6],
                                      ins.1bp.context.rc.start[,7], ins.1bp.context.rc.start[,8], ins.1bp.context.rc.start[,9],
                                      ins.1bp.context.rc.start[,10], sep = '')
    ins.1bp.context.rc.end <- str_split_fixed(ins.1bp.context.rc, '', 11)[,11,drop=F]
    ins.1bp <- cbind(ins.1bp[,1:5,drop=F], paste(ins.1bp.context.rc.start, ins.1bp.context.rc.middle, ins.1bp.context.rc.end, sep = ''))
    colnames(ins.1bp)[6] <- 'CONTEXT RC'
    
    ## summarise 1 bp insertions in matrix format
    ins.1bp.summary <- matrix(0, ncol = 2, nrow = 6)
    colnames(ins.1bp.summary) <- c('C', 'T') ## pyrimidine-centred inserted base
    rownames(ins.1bp.summary) <- c('0 bp', '1 bp', '2 bp', '3 bp', '4 bp', '5+ bp') ## contextual homopolymer-length
    fw <- str_split_fixed(str_split_fixed(ins.1bp[,'CONTEXT FW'], '\\[', 2)[,2,drop=F], '\\]', 2)[,1,drop=F]
    rc <- str_split_fixed(str_split_fixed(ins.1bp[,'CONTEXT RC'], '\\[', 2)[,2,drop=F], '\\]', 2)[,1,drop=F]
    fw.if <- fw == 'C' | fw == 'T'
    for (i in 1:nrow(ins.1bp)){
      
      if(fw.if[i] == T){
        
        ## look at forward context
        upstream.tmp <- str_split_fixed(str_split_fixed(as.character(ins.1bp[i,'CONTEXT FW']), '\\[', 2)[,1], '', 10)
        upstream.tmp <- which(upstream.tmp == fw[i])
        
        ### check all 6 categories for upstream bases
        if(any(upstream.tmp %in% 10)){
          
          if(any(upstream.tmp %in% 9)){
            
            if(any(upstream.tmp %in% 8)){
              
              if(any(upstream.tmp %in% 7)){
                
                if(any(upstream.tmp %in% 6)){
                  
                  upstream.tmp <- 6
                  
                }else{
                  upstream.tmp <- 5
                }
                
              }else{
                upstream.tmp <- 4
              }
              
            }else{
              upstream.tmp <- 3
            }
            
          }else{
            upstream.tmp <- 2
          }
          
        }else{
          upstream.tmp <- 1
        }
        
        downstream.tmp <- str_split_fixed(str_split_fixed(as.character(ins.1bp[i,'CONTEXT FW']), '\\]', 2)[,2], '', 10)
        downstream.tmp <- which(downstream.tmp == fw[i])
        
        ### check all 6 categories for downstream bases
        if(any(downstream.tmp %in% 1)){
          
          if(any(downstream.tmp %in% 2)){
            
            if(any(downstream.tmp %in% 3)){
              
              if(any(downstream.tmp %in% 4)){
                
                if(any(downstream.tmp %in% 5)){
                  
                  downstream.tmp <- 6
                  
                }else{
                  downstream.tmp <- 5
                }
                
              }else{
                downstream.tmp <- 4
              }
              
            }else{
              downstream.tmp <- 3
            }
            
          }else{
            downstream.tmp <- 2
          }
          
        }else{
          downstream.tmp <- 1
        }
        
        ## summarise, which homopolymer is longer
        ins.1bp.summary[max(c(upstream.tmp, downstream.tmp)),fw[i]] <- ins.1bp.summary[max(c(upstream.tmp, downstream.tmp)),fw[i]] + 1
        
      } else {
        
        ## look at reverse complement context
        upstream.tmp <- str_split_fixed(str_split_fixed(as.character(ins.1bp[i,'CONTEXT RC']), '\\[', 2)[,1], '', 10)
        upstream.tmp <- which(upstream.tmp == rc[i])
        
        ### check all 6 categories for upstream bases
        if(any(upstream.tmp %in% 10)){
          
          if(any(upstream.tmp %in% 9)){
            
            if(any(upstream.tmp %in% 8)){
              
              if(any(upstream.tmp %in% 7)){
                
                if(any(upstream.tmp %in% 6)){
                  
                  upstream.tmp <- 6
                  
                }else{
                  upstream.tmp <- 5
                }
                
              }else{
                upstream.tmp <- 4
              }
              
            }else{
              upstream.tmp <- 3
            }
            
          }else{
            upstream.tmp <- 2
          }
          
        }else{
          upstream.tmp <- 1
        }
        
        downstream.tmp <- str_split_fixed(str_split_fixed(as.character(ins.1bp[i,'CONTEXT RC']), '\\]', 2)[,2], '', 10)
        downstream.tmp <- which(downstream.tmp == rc[i])
        
        ### check all 6 categories for downstream bases
        if(any(downstream.tmp %in% 1)){
          
          if(any(downstream.tmp %in% 2)){
            
            if(any(downstream.tmp %in% 3)){
              
              if(any(downstream.tmp %in% 4)){
                
                if(any(downstream.tmp %in% 5)){
                  
                  downstream.tmp <- 6
                  
                }else{
                  downstream.tmp <- 5
                }
                
              }else{
                downstream.tmp <- 4
              }
              
            }else{
              downstream.tmp <- 3
            }
            
          }else{
            downstream.tmp <- 2
          }
          
        }else{
          downstream.tmp <- 1
        }
        
        ## summarise, which homopolymer is longer
        ins.1bp.summary[max(c(upstream.tmp, downstream.tmp)),rc[i]] <- ins.1bp.summary[max(c(upstream.tmp, downstream.tmp)),rc[i]] + 1
      }
      
    } 
    
  }else{
    
    ## summarise 1 bp insertions in matrix format
    ins.1bp.summary <- matrix(0, ncol = 2, nrow = 6)
    colnames(ins.1bp.summary) <- c('C', 'T') ## pyrimidine-centred inserted base
    rownames(ins.1bp.summary) <- c('0 bp', '1 bp', '2 bp', '3 bp', '4 bp', '5+ bp') ## contextual homopolymer-length
    
  }
  
  ## 3. classify >=2 bp deletions into:
  
  # (i) 2 bp deletions at simple repeats (length 1 == "no neighbouring simple repeat")
  if(nrow(dels.2bp) > 0){
    
    ## extract 5 x 2 bp upstream/downstream sequence context from reference
    dels.2bp.context <- as.character(subseq(x = reference[as.character(dels.2bp[,'CHROM'])], 
                                            start = as.numeric(dels.2bp[,'POS']) - 9, 
                                            end = as.numeric(dels.2bp[,'POS']) + 2 + 10))
    dels.2bp.context.middle <- str_split_fixed(dels.2bp.context, '', 22)[,11:12,drop=F]
    dels.2bp.context.middle <- paste0('[', dels.2bp.context.middle[,1], dels.2bp.context.middle[,2], ']')
    dels.2bp.context.start <- str_split_fixed(dels.2bp.context, '', 11)[,1:10,drop=F]
    dels.2bp.context.start <- paste(dels.2bp.context.start[,1], dels.2bp.context.start[,2], dels.2bp.context.start[,3],
                                    dels.2bp.context.start[,4], dels.2bp.context.start[,5], dels.2bp.context.start[,6],
                                    dels.2bp.context.start[,7], dels.2bp.context.start[,8], dels.2bp.context.start[,9],
                                    dels.2bp.context.start[,10], sep = '')
    dels.2bp.context.end <- str_split_fixed(dels.2bp.context, '', 13)[,13,drop=F]
    dels.2bp[,'TRIPLET'] <- paste(dels.2bp.context.start, dels.2bp.context.middle, dels.2bp.context.end, sep = '')
    colnames(dels.2bp)[5] <- 'CONTEXT' 
    
  }
  
  # (ii) 3 bp deletions at simple repeats (length 1 == "no neighbouring simple repeat")
  if(nrow(dels.3bp) > 0){
    
    ## extract 5 x 3 bp upstream/downstream sequence context from reference
    dels.3bp.context <- as.character(subseq(x = reference[as.character(dels.3bp[,'CHROM'])], 
                                            start = as.numeric(dels.3bp[,'POS']) - 14, 
                                            end = as.numeric(dels.3bp[,'POS']) + 3 + 15))
    dels.3bp.context.middle <- str_split_fixed(dels.3bp.context, '', 33)[,16:18,drop=F]
    dels.3bp.context.middle <- paste0('[', dels.3bp.context.middle[,1], dels.3bp.context.middle[,2], dels.3bp.context.middle[,3], ']')
    dels.3bp.context.start <- str_split_fixed(dels.3bp.context, '', 33)[,1:15,drop=F]
    dels.3bp.context.start <- paste(dels.3bp.context.start[,1], dels.3bp.context.start[,2], dels.3bp.context.start[,3],
                                    dels.3bp.context.start[,4], dels.3bp.context.start[,5], dels.3bp.context.start[,6],
                                    dels.3bp.context.start[,7], dels.3bp.context.start[,8], dels.3bp.context.start[,9],
                                    dels.3bp.context.start[,10], dels.3bp.context.start[,11], dels.3bp.context.start[,12], 
                                    dels.3bp.context.start[,13], dels.3bp.context.start[,14], dels.3bp.context.start[,15], sep = '')
    dels.3bp.context.end <- str_split_fixed(dels.3bp.context, '', 19)[,19,drop=F]
    dels.3bp[,'TRIPLET'] <- paste(dels.3bp.context.start, dels.3bp.context.middle, dels.3bp.context.end, sep = '')
    colnames(dels.3bp)[5] <- 'CONTEXT' 
    
  }
  
  # (iii) 4 bp deletions at simple repeats (length 1 == "no neighbouring simple repeat")
  if(nrow(dels.4bp) > 0){
    
    ## extract 5 x 4 bp upstream/downstream sequence context from reference
    dels.4bp.context <- as.character(subseq(x = reference[as.character(dels.4bp[,'CHROM'])], 
                                            start = as.numeric(dels.4bp[,'POS']) - 19, 
                                            end = as.numeric(dels.4bp[,'POS']) + 4 + 20))
    dels.4bp.context.middle <- str_split_fixed(dels.4bp.context, '', 41)[,21:24,drop=F]
    dels.4bp.context.middle <- paste0('[', dels.4bp.context.middle[,1], 
                                      dels.4bp.context.middle[,2], 
                                      dels.4bp.context.middle[,3], 
                                      dels.4bp.context.middle[,4], ']')
    dels.4bp.context.start <- str_split_fixed(dels.4bp.context, '', 21)[,1:20,drop=F]
    dels.4bp.context.start <- paste(dels.4bp.context.start[,1], dels.4bp.context.start[,2], dels.4bp.context.start[,3],
                                    dels.4bp.context.start[,4], dels.4bp.context.start[,5], dels.4bp.context.start[,6],
                                    dels.4bp.context.start[,7], dels.4bp.context.start[,8], dels.4bp.context.start[,9],
                                    dels.4bp.context.start[,10], dels.4bp.context.start[,11], dels.4bp.context.start[,12], 
                                    dels.4bp.context.start[,13], dels.4bp.context.start[,14], dels.4bp.context.start[,15], 
                                    dels.4bp.context.start[,16], dels.4bp.context.start[,17], dels.4bp.context.start[,18], 
                                    dels.4bp.context.start[,19], dels.4bp.context.start[,20], sep = '')
    dels.4bp.context.end <- str_split_fixed(dels.4bp.context, '', 25)[,25,drop=F]
    dels.4bp[,'TRIPLET'] <- paste(dels.4bp.context.start, dels.4bp.context.middle, dels.4bp.context.end, sep = '')
    colnames(dels.4bp)[5] <- 'CONTEXT' 
    
  }
  
  # (iv) 5+ bp deletions at simple repeats (length 1 == "no neighbouring simple repeat")
  if(nrow(dels.5bp) > 0){
    
    ## extract 5 x 100 bp upstream/downstream sequence context from reference
    ## i.e. max. 100 bp repeat motif
    dels.5bp.context <- as.character(subseq(x = reference[as.character(dels.5bp[,'CHROM'])], 
                                            start = as.numeric(dels.5bp[,'POS']) - 499, 
                                            end = as.numeric(dels.5bp[,'POS']) + 5 + 500))
    
    ### account for the different repeat lengths: iterate
    dels.5bp.context.middle <- as.character(dels.5bp[,'REF'])
    dels.5bp.context.middle <- str_split_fixed(dels.5bp.context.middle, '', 2)[,2,drop=F]
    dels.5bp.context.middle <- paste0('[', dels.5bp.context.middle, ']')
    dels.5bp.context.middle.lengths <- nchar(as.character(dels.5bp[,'REF'])) - 1
    dels.5bp.context.start <- rep(NA, nrow(dels.5bp))
    for (i in 1:length(dels.5bp.context.start)){
      tmp.dels.5bp.context.start <- str_split_fixed(dels.5bp.context[i], '', 501)[,c(1 + 500-c(5*dels.5bp.context.middle.lengths[i])):500,drop=F]
      dels.5bp.context.start[i] <- paste(tmp.dels.5bp.context.start, collapse = '')
    }
    dels.5bp.context.end <- rep(NA, nrow(dels.5bp))
    for (i in 1:length(dels.5bp.context.end)){
      tmp.dels.5bp.context.end <- str_split_fixed(dels.5bp.context[i], '', 1006)[,c(501+dels.5bp.context.middle.lengths[i]):c(501 + dels.5bp.context.middle.lengths[i]*6 - 1),drop=F]
      dels.5bp.context.end[i] <- paste(tmp.dels.5bp.context.end, collapse = '')
    }
    dels.5bp[,'TRIPLET'] <- paste(dels.5bp.context.start, dels.5bp.context.middle, dels.5bp.context.end, sep = '')
    colnames(dels.5bp)[5] <- 'CONTEXT'
    
  }
  
  ## summarise >=2 bp deletions at simple repeats in matrix format
  ## also create sub-tables for repeat length = 1 hits, to later assess these for microhomologies
  dels.greater.2bp.summary <- matrix(0, ncol = 4, nrow = 6)
  colnames(dels.greater.2bp.summary) <- c('2 bp', '3 bp', '4 bp', '5+ bp') ## deletion size
  rownames(dels.greater.2bp.summary) <- c('1 or MH', '2', '3', '4', '5', '6+') ## number of repeats
  
  ## 2 bp
  if(nrow(dels.2bp) > 0){
    
    repeat.nts <- str_split_fixed(str_split_fixed(dels.2bp[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,1,drop=F]
    downstream.context <- str_split_fixed(str_split_fixed(dels.2bp[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,2,drop=F]
    dels.2bp.pot.MH.ind <- c()
    for (i in 1:nrow(dels.2bp)){
      
      ## look at repeat
      tmp.repeat.nts <- repeat.nts[i]
      tmp.repeat.length <- nchar(tmp.repeat.nts)
      
      ## how often does it match consecutively in the immediate downstream context?
      tmp.downstream.context <- downstream.context[i]
      tmp.downstream.context <- strsplit(tmp.downstream.context, '')[[1]]
      
      ## group
      tmp.downstream.context <- c(paste(tmp.downstream.context[c(tmp.repeat.length-c(tmp.repeat.length - 1)):tmp.repeat.length], collapse = ''),
                                  paste(tmp.downstream.context[c(2*tmp.repeat.length-c(tmp.repeat.length - 1)):c(2*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(3*tmp.repeat.length-c(tmp.repeat.length - 1)):c(3*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(4*tmp.repeat.length-c(tmp.repeat.length - 1)):c(4*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(5*tmp.repeat.length-c(tmp.repeat.length - 1)):c(5*tmp.repeat.length)], collapse = ''))
      
      ### check all 6 categories for upstream bases
      if(tmp.downstream.context[1] %in% tmp.repeat.nts){
        
        if(tmp.downstream.context[2] %in% tmp.repeat.nts){
          
          if(tmp.downstream.context[3] %in% tmp.repeat.nts){
            
            if(tmp.downstream.context[4] %in% tmp.repeat.nts){
              
              if(tmp.downstream.context[5] %in% tmp.repeat.nts){
                
                dels.greater.2bp.summary['6+', '2 bp'] <- dels.greater.2bp.summary['6+', '2 bp'] + 1
                
              } else{
                
                dels.greater.2bp.summary['5', '2 bp'] <- dels.greater.2bp.summary['5', '2 bp'] + 1
                
              }
              
            }else{
              
              dels.greater.2bp.summary['4', '2 bp'] <- dels.greater.2bp.summary['4', '2 bp'] + 1
              
            }
            
          }else{
            
            dels.greater.2bp.summary['3', '2 bp'] <- dels.greater.2bp.summary['3', '2 bp'] + 1
            
          }
          
        }else{
          
          dels.greater.2bp.summary['2', '2 bp'] <- dels.greater.2bp.summary['2', '2 bp'] + 1
          
        }
        
      }else{
        
        dels.greater.2bp.summary['1 or MH', '2 bp'] <- dels.greater.2bp.summary['1 or MH', '2 bp'] + 1
        dels.2bp.pot.MH.ind <- c(dels.2bp.pot.MH.ind, i)
        
      }
      
    }
    dels.2bp.pot.MH <- dels.2bp[dels.2bp.pot.MH.ind,,drop=F] 
    
  }
  
  ## 3 bp
  if(nrow(dels.3bp) > 0){
    
    repeat.nts <- str_split_fixed(str_split_fixed(dels.3bp[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,1,drop=F]
    downstream.context <- str_split_fixed(str_split_fixed(dels.3bp[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,2,drop=F]
    dels.3bp.pot.MH.ind <- c()
    for (i in 1:nrow(dels.3bp)){
      
      ## look at repeat
      tmp.repeat.nts <- repeat.nts[i]
      tmp.repeat.length <- nchar(tmp.repeat.nts)
      
      ## how often does it match consecutively in the immediate downstream context?
      tmp.downstream.context <- downstream.context[i]
      tmp.downstream.context <- strsplit(tmp.downstream.context, '')[[1]]
      
      ## group
      tmp.downstream.context <- c(paste(tmp.downstream.context[c(tmp.repeat.length-c(tmp.repeat.length - 1)):tmp.repeat.length], collapse = ''),
                                  paste(tmp.downstream.context[c(2*tmp.repeat.length-c(tmp.repeat.length - 1)):c(2*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(3*tmp.repeat.length-c(tmp.repeat.length - 1)):c(3*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(4*tmp.repeat.length-c(tmp.repeat.length - 1)):c(4*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(5*tmp.repeat.length-c(tmp.repeat.length - 1)):c(5*tmp.repeat.length)], collapse = ''))
      
      ### check all 6 categories for upstream bases
      if(tmp.downstream.context[1] %in% tmp.repeat.nts){
        
        if(tmp.downstream.context[2] %in% tmp.repeat.nts){
          
          if(tmp.downstream.context[3] %in% tmp.repeat.nts){
            
            if(tmp.downstream.context[4] %in% tmp.repeat.nts){
              
              if(tmp.downstream.context[5] %in% tmp.repeat.nts){
                
                dels.greater.2bp.summary['6+', '3 bp'] <- dels.greater.2bp.summary['6+', '3 bp'] + 1
                
              } else{
                
                dels.greater.2bp.summary['5', '3 bp'] <- dels.greater.2bp.summary['5', '3 bp'] + 1
                
              }
              
            }else{
              
              dels.greater.2bp.summary['4', '3 bp'] <- dels.greater.2bp.summary['4', '3 bp'] + 1
              
            }
            
          }else{
            
            dels.greater.2bp.summary['3', '3 bp'] <- dels.greater.2bp.summary['3', '3 bp'] + 1
            
          }
          
        }else{
          
          dels.greater.2bp.summary['2', '3 bp'] <- dels.greater.2bp.summary['2', '3 bp'] + 1
          
        }
        
      }else{
        
        dels.greater.2bp.summary['1 or MH', '3 bp'] <- dels.greater.2bp.summary['1 or MH', '3 bp'] + 1
        dels.3bp.pot.MH.ind <- c(dels.3bp.pot.MH.ind, i)
        
      }
      
    }
    dels.3bp.pot.MH <- dels.3bp[dels.3bp.pot.MH.ind,,drop=F] 
    
  }
  
  ## 4 bp
  if(nrow(dels.4bp) > 0){
    
    repeat.nts <- str_split_fixed(str_split_fixed(dels.4bp[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,1,drop=F]
    downstream.context <- str_split_fixed(str_split_fixed(dels.4bp[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,2,drop=F]
    dels.4bp.pot.MH.ind <- c()
    for (i in 1:nrow(dels.4bp)){
      
      ## look at repeat
      tmp.repeat.nts <- repeat.nts[i]
      tmp.repeat.length <- nchar(tmp.repeat.nts)
      
      ## how often does it match consecutively in the immediate downstream context?
      tmp.downstream.context <- downstream.context[i]
      tmp.downstream.context <- strsplit(tmp.downstream.context, '')[[1]]
      
      ## group
      tmp.downstream.context <- c(paste(tmp.downstream.context[c(tmp.repeat.length-c(tmp.repeat.length - 1)):tmp.repeat.length], collapse = ''),
                                  paste(tmp.downstream.context[c(2*tmp.repeat.length-c(tmp.repeat.length - 1)):c(2*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(3*tmp.repeat.length-c(tmp.repeat.length - 1)):c(3*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(4*tmp.repeat.length-c(tmp.repeat.length - 1)):c(4*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(5*tmp.repeat.length-c(tmp.repeat.length - 1)):c(5*tmp.repeat.length)], collapse = ''))
      
      ### check all 6 categories for upstream bases
      if(tmp.downstream.context[1] %in% tmp.repeat.nts){
        
        if(tmp.downstream.context[2] %in% tmp.repeat.nts){
          
          if(tmp.downstream.context[3] %in% tmp.repeat.nts){
            
            if(tmp.downstream.context[4] %in% tmp.repeat.nts){
              
              if(tmp.downstream.context[5] %in% tmp.repeat.nts){
                
                dels.greater.2bp.summary['6+', '4 bp'] <- dels.greater.2bp.summary['6+', '4 bp'] + 1
                
              } else{
                
                dels.greater.2bp.summary['5', '4 bp'] <- dels.greater.2bp.summary['5', '4 bp'] + 1
                
              }
              
            }else{
              
              dels.greater.2bp.summary['4', '4 bp'] <- dels.greater.2bp.summary['4', '4 bp'] + 1
              
            }
            
          }else{
            
            dels.greater.2bp.summary['3', '4 bp'] <- dels.greater.2bp.summary['3', '4 bp'] + 1
            
          }
          
        }else{
          
          dels.greater.2bp.summary['2', '4 bp'] <- dels.greater.2bp.summary['2', '4 bp'] + 1
          
        }
        
      }else{
        
        dels.greater.2bp.summary['1 or MH', '4 bp'] <- dels.greater.2bp.summary['1 or MH', '4 bp'] + 1
        dels.4bp.pot.MH.ind <- c(dels.4bp.pot.MH.ind, i)
        
      }
      
    }
    dels.4bp.pot.MH <- dels.4bp[dels.4bp.pot.MH.ind,,drop=F]
    
  }
  
  ## 5+ bp
  if(nrow(dels.5bp) > 0){
    
    repeat.nts <- str_split_fixed(str_split_fixed(dels.5bp[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,1,drop=F]
    downstream.context <- str_split_fixed(str_split_fixed(dels.5bp[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,2,drop=F]
    dels.5bp.pot.MH.ind <- c()
    for (i in 1:nrow(dels.5bp)){
      
      ## look at repeat
      tmp.repeat.nts <- repeat.nts[i]
      tmp.repeat.length <- nchar(tmp.repeat.nts)
      
      ## how often does it match consecutively in the immediate downstream context?
      tmp.downstream.context <- downstream.context[i]
      tmp.downstream.context <- strsplit(tmp.downstream.context, '')[[1]]
      
      ## group
      tmp.downstream.context <- c(paste(tmp.downstream.context[c(tmp.repeat.length-c(tmp.repeat.length - 1)):tmp.repeat.length], collapse = ''),
                                  paste(tmp.downstream.context[c(2*tmp.repeat.length-c(tmp.repeat.length - 1)):c(2*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(3*tmp.repeat.length-c(tmp.repeat.length - 1)):c(3*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(4*tmp.repeat.length-c(tmp.repeat.length - 1)):c(4*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(5*tmp.repeat.length-c(tmp.repeat.length - 1)):c(5*tmp.repeat.length)], collapse = ''))
      
      ### check all 6 categories for upstream bases
      if(tmp.downstream.context[1] %in% tmp.repeat.nts){
        
        if(tmp.downstream.context[2] %in% tmp.repeat.nts){
          
          if(tmp.downstream.context[3] %in% tmp.repeat.nts){
            
            if(tmp.downstream.context[4] %in% tmp.repeat.nts){
              
              if(tmp.downstream.context[5] %in% tmp.repeat.nts){
                
                dels.greater.2bp.summary['6+', '5+ bp'] <- dels.greater.2bp.summary['6+', '5+ bp'] + 1
                
              } else{
                
                dels.greater.2bp.summary['5', '5+ bp'] <- dels.greater.2bp.summary['5', '5+ bp'] + 1
                
              }
              
            }else{
              
              dels.greater.2bp.summary['4', '5+ bp'] <- dels.greater.2bp.summary['4', '5+ bp'] + 1
              
            }
            
          }else{
            
            dels.greater.2bp.summary['3', '5+ bp'] <- dels.greater.2bp.summary['3', '5+ bp'] + 1
            
          }
          
        }else{
          
          dels.greater.2bp.summary['2', '5+ bp'] <- dels.greater.2bp.summary['2', '5+ bp'] + 1
          
        }
        
      }else{
        
        dels.greater.2bp.summary['1 or MH', '5+ bp'] <- dels.greater.2bp.summary['1 or MH', '5+ bp'] + 1
        dels.5bp.pot.MH.ind <- c(dels.5bp.pot.MH.ind, i)
        
      }
      
    }
    dels.5bp.pot.MH <- dels.5bp[dels.5bp.pot.MH.ind,,drop=F] 
    
  }
  
  # (v) >=2 bp deletions at microhomologies
  ## go directly into in matrix format, thereby also re-classify repeat length = 1 deletions
  dels.greater.2bp.MH.summary <- matrix(0, ncol = 4, nrow = 5)
  colnames(dels.greater.2bp.MH.summary) <- c('2 bp', '3 bp', '4 bp', '5+ bp') ## deletion size
  rownames(dels.greater.2bp.MH.summary) <- c('1 bp MH', '2 bp MH', '3 bp MH', '4 bp MH', '5+ bp MH') ## microhomology length
  dels.greater.2bp.MH.summary[2:5, '2 bp'] <- NA
  dels.greater.2bp.MH.summary[3:5, '3 bp'] <- NA
  dels.greater.2bp.MH.summary[4:5, '4 bp'] <- NA
  
  ## 2 bp
  if(nrow(dels.2bp) > 0){
    
    if(nrow(dels.2bp.pot.MH) > 0){
      
      repeat.nts <- str_split_fixed(str_split_fixed(dels.2bp.pot.MH[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,1,drop=F]
      upstream.context <- str_split_fixed(dels.2bp.pot.MH[,'CONTEXT'], '\\[', 2)[,1,drop=F]
      downstream.context <- str_split_fixed(str_split_fixed(dels.2bp.pot.MH[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,2,drop=F]
      for (i in 1:nrow(dels.2bp.pot.MH)){
        
        ## look at repeat
        tmp.repeat.nts <- repeat.nts[i]
        
        ## isolate first nucleotide in the immediate upstream context
        tmp.upstream.context <- upstream.context[i]
        tmp.upstream.context <- strsplit(tmp.upstream.context, '')[[1]][length(strsplit(tmp.upstream.context, '')[[1]])]
        
        ## isolate first nucleotide in the immediate downstream context
        tmp.downstream.context <- downstream.context[i]
        tmp.downstream.context <- strsplit(tmp.downstream.context, '')[[1]][1]
        
        if(strsplit(tmp.repeat.nts, '')[[1]][2] == tmp.upstream.context | strsplit(tmp.repeat.nts, '')[[1]][1] == tmp.downstream.context){
          
          dels.greater.2bp.MH.summary['1 bp MH', '2 bp'] <- dels.greater.2bp.MH.summary['1 bp MH', '2 bp'] + 1
          dels.greater.2bp.summary['1 or MH', '2 bp'] <- dels.greater.2bp.summary['1 or MH', '2 bp'] - 1
          
        }
      } 
      
    }
    
  }
  
  ## 3 bp
  if(nrow(dels.3bp) > 0){
    
    if(nrow(dels.3bp.pot.MH) > 0){
      
      repeat.nts <- str_split_fixed(str_split_fixed(dels.3bp.pot.MH[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,1,drop=F]
      upstream.context <- str_split_fixed(dels.3bp.pot.MH[,'CONTEXT'], '\\[', 2)[,1,drop=F]
      downstream.context <- str_split_fixed(str_split_fixed(dels.3bp.pot.MH[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,2,drop=F]
      for (i in 1:nrow(dels.3bp.pot.MH)){
        
        ## look at repeat
        tmp.repeat.nts <- repeat.nts[i]
        
        ## isolate first two nucleotides in the immediate upstream context
        tmp.upstream.context <- upstream.context[i]
        tmp.upstream.context <- strsplit(tmp.upstream.context, '')[[1]][c(length(strsplit(tmp.upstream.context, '')[[1]]) - 1):length(strsplit(tmp.upstream.context, '')[[1]])]
        
        ## isolate first two nucleotides in the immediate downstream context
        tmp.downstream.context <- downstream.context[i]
        tmp.downstream.context <- strsplit(tmp.downstream.context, '')[[1]][1:2]
        
        ### MH length = 2 (important to start with the highest possible MH length)
        
        if(all(c(strsplit(tmp.repeat.nts, '')[[1]][2:3] == tmp.upstream.context) == T) | all(c(strsplit(tmp.repeat.nts, '')[[1]][1:2] == tmp.downstream.context) == T)){
          
          dels.greater.2bp.MH.summary['2 bp MH', '3 bp'] <- dels.greater.2bp.MH.summary['2 bp MH', '3 bp'] + 1
          dels.greater.2bp.summary['1 or MH', '3 bp'] <- dels.greater.2bp.summary['1 or MH', '3 bp'] - 1
          
        }else{
          
          ### MH length = 1
          
          if(strsplit(tmp.repeat.nts, '')[[1]][3] == tmp.upstream.context[2] | strsplit(tmp.repeat.nts, '')[[1]][1] == tmp.downstream.context[1]){
            
            dels.greater.2bp.MH.summary['1 bp MH', '3 bp'] <- dels.greater.2bp.MH.summary['1 bp MH', '3 bp'] + 1
            dels.greater.2bp.summary['1 or MH', '3 bp'] <- dels.greater.2bp.summary['1 or MH', '3 bp'] - 1
            
          }
          
        }
        
      } 
      
    }
    
  }
  
  ## 4 bp
  if(nrow(dels.4bp) > 0){
    
    if(nrow(dels.4bp.pot.MH) > 0){
      
      repeat.nts <- str_split_fixed(str_split_fixed(dels.4bp.pot.MH[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,1,drop=F]
      upstream.context <- str_split_fixed(dels.4bp.pot.MH[,'CONTEXT'], '\\[', 2)[,1,drop=F]
      downstream.context <- str_split_fixed(str_split_fixed(dels.4bp.pot.MH[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,2,drop=F]
      for (i in 1:nrow(dels.4bp.pot.MH)){
        
        ## look at repeat
        tmp.repeat.nts <- repeat.nts[i]
        
        ## isolate first three nucleotides in the immediate upstream context
        tmp.upstream.context <- upstream.context[i]
        tmp.upstream.context <- strsplit(tmp.upstream.context, '')[[1]][c(length(strsplit(tmp.upstream.context, '')[[1]]) - 2):length(strsplit(tmp.upstream.context, '')[[1]])]
        
        ## isolate first three nucleotides in the immediate downstream context
        tmp.downstream.context <- downstream.context[i]
        tmp.downstream.context <- strsplit(tmp.downstream.context, '')[[1]][1:3]
        
        ### MH length = 3 (important to start with the highest possible MH length)
        
        if(all(c(strsplit(tmp.repeat.nts, '')[[1]][2:4] == tmp.upstream.context) == T) | all(c(strsplit(tmp.repeat.nts, '')[[1]][1:3] == tmp.downstream.context) == T)){
          
          dels.greater.2bp.MH.summary['3 bp MH', '4 bp'] <- dels.greater.2bp.MH.summary['3 bp MH', '4 bp'] + 1
          dels.greater.2bp.summary['1 or MH', '4 bp'] <- dels.greater.2bp.summary['1 or MH', '4 bp'] - 1
          
        }else{
          
          ### MH length = 2
          
          if(all(c(strsplit(tmp.repeat.nts, '')[[1]][3:4] == tmp.upstream.context[2:3]) ==T) | all(c(strsplit(tmp.repeat.nts, '')[[1]][1:2] == tmp.downstream.context[1:2]) == T)){
            
            dels.greater.2bp.MH.summary['2 bp MH', '4 bp'] <- dels.greater.2bp.MH.summary['2 bp MH', '4 bp'] + 1
            dels.greater.2bp.summary['1 or MH', '4 bp'] <- dels.greater.2bp.summary['1 or MH', '4 bp'] - 1
            
          }else{
            
            ### MH length = 1
            
            if(strsplit(tmp.repeat.nts, '')[[1]][4] == tmp.upstream.context[3] | strsplit(tmp.repeat.nts, '')[[1]][1] == tmp.downstream.context[1]){
              
              dels.greater.2bp.MH.summary['1 bp MH', '4 bp'] <- dels.greater.2bp.MH.summary['1 bp MH', '4 bp'] + 1
              dels.greater.2bp.summary['1 or MH', '4 bp'] <- dels.greater.2bp.summary['1 or MH', '4 bp'] - 1
              
            }
            
          }
          
        }
        
      } 
      
    }
    
  }
  
  ## 5+ bp
  if(nrow(dels.5bp) > 0){
    
    if(nrow(dels.5bp.pot.MH) > 0){
      
      repeat.nts <- str_split_fixed(str_split_fixed(dels.5bp.pot.MH[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,1,drop=F]
      upstream.context <- str_split_fixed(dels.5bp.pot.MH[,'CONTEXT'], '\\[', 2)[,1,drop=F]
      downstream.context <- str_split_fixed(str_split_fixed(dels.5bp.pot.MH[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,2,drop=F]
      for (i in 1:nrow(dels.5bp.pot.MH)){
        
        ## look at repeat
        tmp.repeat.nts <- repeat.nts[i]
        tmp.repeat.length <- nchar(tmp.repeat.nts)
        
        ## isolate first X (= repeat length - 1) nucleotides in the immediate upstream context
        tmp.upstream.context <- upstream.context[i]
        tmp.upstream.context <- strsplit(tmp.upstream.context, '')[[1]][c(length(strsplit(tmp.upstream.context, '')[[1]]) - 4):length(strsplit(tmp.upstream.context, '')[[1]])]
        
        ## isolate first X (= repeat length - 1) nucleotides in the immediate downstream context
        tmp.downstream.context <- downstream.context[i]
        tmp.downstream.context <- strsplit(tmp.downstream.context, '')[[1]][1:5]
        
        ### MH length = 5+ (important to start with the highest possible MH length); only need to look at first 5 up/downstream NTs
        
        if(all(c(strsplit(tmp.repeat.nts, '')[[1]][c(tmp.repeat.length-4):tmp.repeat.length] == tmp.upstream.context[c(length(tmp.upstream.context) - 4):length(tmp.upstream.context)]) == T) | all(c(strsplit(tmp.repeat.nts, '')[[1]][1:5] == tmp.downstream.context[1:5]) == T)){
          
          dels.greater.2bp.MH.summary['5+ bp MH', '5+ bp'] <- dels.greater.2bp.MH.summary['5+ bp MH', '5+ bp'] + 1
          dels.greater.2bp.summary['1 or MH', '5+ bp'] <- dels.greater.2bp.summary['1 or MH', '5+ bp'] - 1
          
        }else{
          
          ### MH length = 4
          
          if(all(c(strsplit(tmp.repeat.nts, '')[[1]][c(tmp.repeat.length-3):tmp.repeat.length] == tmp.upstream.context[c(length(tmp.upstream.context) - 3):length(tmp.upstream.context)]) == T) | all(c(strsplit(tmp.repeat.nts, '')[[1]][1:4] == tmp.downstream.context[1:4]) == T)){
            
            dels.greater.2bp.MH.summary['4 bp MH', '5+ bp'] <- dels.greater.2bp.MH.summary['4 bp MH', '5+ bp'] + 1
            dels.greater.2bp.summary['1 or MH', '5+ bp'] <- dels.greater.2bp.summary['1 or MH', '5+ bp'] - 1
            
          }else{
            
            ### MH length = 3
            
            if(all(c(strsplit(tmp.repeat.nts, '')[[1]][c(tmp.repeat.length-2):tmp.repeat.length] == tmp.upstream.context[c(length(tmp.upstream.context) - 2):length(tmp.upstream.context)]) == T) | all(c(strsplit(tmp.repeat.nts, '')[[1]][1:3] == tmp.downstream.context[1:3]) == T)){
              
              dels.greater.2bp.MH.summary['3 bp MH', '5+ bp'] <- dels.greater.2bp.MH.summary['3 bp MH', '5+ bp'] + 1
              dels.greater.2bp.summary['1 or MH', '5+ bp'] <- dels.greater.2bp.summary['1 or MH', '5+ bp'] - 1
              
            }else{
              
              ### MH length = 2
              
              if(all(c(strsplit(tmp.repeat.nts, '')[[1]][c(tmp.repeat.length-1):tmp.repeat.length] == tmp.upstream.context[c(length(tmp.upstream.context) - 1):length(tmp.upstream.context)]) == T) | all(c(strsplit(tmp.repeat.nts, '')[[1]][1:2] == tmp.downstream.context[1:2]) == T)){
                
                dels.greater.2bp.MH.summary['2 bp MH', '5+ bp'] <- dels.greater.2bp.MH.summary['2 bp MH', '5+ bp'] + 1
                dels.greater.2bp.summary['1 or MH', '5+ bp'] <- dels.greater.2bp.summary['1 or MH', '5+ bp'] - 1
                
              }else{
                
                if(strsplit(tmp.repeat.nts, '')[[1]][tmp.repeat.length] == tmp.upstream.context[length(tmp.upstream.context)] | strsplit(tmp.repeat.nts, '')[[1]][1] == tmp.downstream.context[1]){
                  
                  dels.greater.2bp.MH.summary['1 bp MH', '5+ bp'] <- dels.greater.2bp.MH.summary['1 bp MH', '5+ bp'] + 1
                  dels.greater.2bp.summary['1 or MH', '5+ bp'] <- dels.greater.2bp.summary['1 or MH', '5+ bp'] - 1
                  
                }
                
              }
              
            }
            
          }
          
        }
        
      } 
      
    }
    
  }
  
  ## after MH cases have been taken out, rename row of deletions at 1-unit repeats in original table
  rownames(dels.greater.2bp.summary)[1] <- '1'
  
  ## 4. classify >=2 bp insertions into:
  
  # (i) 2 bp insertions at simple repeats (length 0 == "no neighbouring simple repeat")
  if(nrow(ins.2bp) > 0){
    
    ins.2bp.context <- as.character(subseq(x = reference[as.character(ins.2bp[,'CHROM'])], 
                                           start = as.numeric(ins.2bp[,'POS']) - 9, 
                                           end = as.numeric(ins.2bp[,'POS']) + 10))
    ins.2bp.context.middle <- as.character(ins.2bp[,'ALT'])
    ins.2bp.context.middle <- paste0('[', str_split_fixed(ins.2bp.context.middle, '', 2)[,2,drop=F], ']')
    ins.2bp.context.start <- str_split_fixed(ins.2bp.context, '', 11)[,1:10,drop=F]
    ins.2bp.context.start <- paste(ins.2bp.context.start[,1], ins.2bp.context.start[,2], ins.2bp.context.start[,3],
                                   ins.2bp.context.start[,4], ins.2bp.context.start[,5], ins.2bp.context.start[,6],
                                   ins.2bp.context.start[,7], ins.2bp.context.start[,8], ins.2bp.context.start[,9],
                                   ins.2bp.context.start[,10], sep = '')
    ins.2bp.context.end <- str_split_fixed(ins.2bp.context, '', 11)[,11,drop=F]
    ins.2bp[,'TRIPLET'] <- paste(ins.2bp.context.start, ins.2bp.context.middle, ins.2bp.context.end, sep = '')
    colnames(ins.2bp)[5] <- 'CONTEXT' 
    
  }
  
  # (ii) 3 bp insertions at simple repeats (length 0 == "no neighbouring simple repeat")
  if(nrow(ins.3bp) > 0){
    
    ins.3bp.context <- as.character(subseq(x = reference[as.character(ins.3bp[,'CHROM'])], 
                                           start = as.numeric(ins.3bp[,'POS']) - 14, 
                                           end = as.numeric(ins.3bp[,'POS']) + 15))
    ins.3bp.context.middle <- as.character(ins.3bp[,'ALT'])
    ins.3bp.context.middle <- paste0('[', str_split_fixed(ins.3bp.context.middle, '', 2)[,2,drop=F], ']')
    ins.3bp.context.start <- str_split_fixed(ins.3bp.context, '', 16)[,1:15,drop=F]
    ins.3bp.context.start <- paste(ins.3bp.context.start[,1], ins.3bp.context.start[,2], ins.3bp.context.start[,3],
                                   ins.3bp.context.start[,4], ins.3bp.context.start[,5], ins.3bp.context.start[,6],
                                   ins.3bp.context.start[,7], ins.3bp.context.start[,8], ins.3bp.context.start[,9],
                                   ins.3bp.context.start[,10], ins.3bp.context.start[,11], ins.3bp.context.start[,12], 
                                   ins.3bp.context.start[,13], ins.3bp.context.start[,14], ins.3bp.context.start[,15], sep = '')
    ins.3bp.context.end <- str_split_fixed(ins.3bp.context, '', 16)[,16,drop=F]
    ins.3bp[,'TRIPLET'] <- paste(ins.3bp.context.start, ins.3bp.context.middle, ins.3bp.context.end, sep = '')
    colnames(ins.3bp)[5] <- 'CONTEXT' 
    
  }
  
  # (iii) 4 bp insertions at simple repeats (length 0 == "no neighbouring simple repeat")
  if(nrow(ins.4bp) > 0){
    
    ins.4bp.context <- as.character(subseq(x = reference[as.character(ins.4bp[,'CHROM'])], 
                                           start = as.numeric(ins.4bp[,'POS']) - 19, 
                                           end = as.numeric(ins.4bp[,'POS']) + 20))
    ins.4bp.context.middle <- as.character(ins.4bp[,'ALT'])
    ins.4bp.context.middle <- paste0('[', str_split_fixed(ins.4bp.context.middle, '', 2)[,2,drop=F], ']')
    ins.4bp.context.start <- str_split_fixed(ins.4bp.context, '', 21)[,1:20,drop=F]
    ins.4bp.context.start <- paste(ins.4bp.context.start[,1], ins.4bp.context.start[,2], ins.4bp.context.start[,3],
                                   ins.4bp.context.start[,4], ins.4bp.context.start[,5], ins.4bp.context.start[,6],
                                   ins.4bp.context.start[,7], ins.4bp.context.start[,8], ins.4bp.context.start[,9],
                                   ins.4bp.context.start[,10], ins.4bp.context.start[,11], ins.4bp.context.start[,12], 
                                   ins.4bp.context.start[,13], ins.4bp.context.start[,14], ins.4bp.context.start[,15], 
                                   ins.4bp.context.start[,16], ins.4bp.context.start[,17], ins.4bp.context.start[,18], 
                                   ins.4bp.context.start[,19], ins.4bp.context.start[,20], sep = '')
    ins.4bp.context.end <- str_split_fixed(ins.4bp.context, '', 21)[,21,drop=F]
    ins.4bp[,'TRIPLET'] <- paste(ins.4bp.context.start, ins.4bp.context.middle, ins.4bp.context.end, sep = '')
    colnames(ins.4bp)[5] <- 'CONTEXT' 
    
  }
  
  # (iv) 5+ bp insertions at simple repeats (length 0 == "no neighbouring simple repeat")
  if(nrow(ins.5bp) > 0){
    
    ins.5bp.context <- as.character(subseq(x = reference[as.character(ins.5bp[,'CHROM'])], 
                                           start = as.numeric(ins.5bp[,'POS']) - 499, 
                                           end = as.numeric(ins.5bp[,'POS']) + 500))
    
    ### account for the different repeat lengths: iterate
    ins.5bp.context.middle <- as.character(ins.5bp[,'ALT'])
    ins.5bp.context.middle <- str_split_fixed(ins.5bp.context.middle, '', 2)[,2,drop=F]
    ins.5bp.context.middle <- paste0('[', ins.5bp.context.middle, ']')
    ins.5bp.context.middle.lengths <- nchar(as.character(ins.5bp[,'ALT'])) - 1
    ins.5bp.context.start <- rep(NA, nrow(ins.5bp))
    for (i in 1:length(ins.5bp.context.start)){
      tmp.ins.5bp.context.start <- str_split_fixed(ins.5bp.context[i], '', 501)[,c(1 + 500-c(5*ins.5bp.context.middle.lengths[i])):500,drop=F]
      ins.5bp.context.start[i] <- paste(tmp.ins.5bp.context.start, collapse = '')
    }
    ins.5bp.context.end <- rep(NA, nrow(ins.5bp))
    for (i in 1:length(ins.5bp.context.end)){
      tmp.ins.5bp.context.end <- str_split_fixed(ins.5bp.context[i], '', 1000)[,501:c(501 + ins.5bp.context.middle.lengths[i]*6 - 1),drop=F]
      ins.5bp.context.end[i] <- paste(tmp.ins.5bp.context.end, collapse = '')
    }
    ins.5bp[,'TRIPLET'] <- paste(ins.5bp.context.start, ins.5bp.context.middle, ins.5bp.context.end, sep = '')
    colnames(ins.5bp)[5] <- 'CONTEXT'
    
  }
  
  ## summarise >=2 bp insertions at simple repeats in matrix format
  ins.greater.2bp.summary <- matrix(0, ncol = 4, nrow = 6)
  colnames(ins.greater.2bp.summary) <- c('2 bp', '3 bp', '4 bp', '5+ bp') ## insertion size
  rownames(ins.greater.2bp.summary) <- c('0', '1', '2', '3', '4', '5+') ## number of repeats
  
  ## 2 bp
  if(nrow(ins.2bp) > 0){
    
    repeat.nts <- str_split_fixed(str_split_fixed(ins.2bp[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,1,drop=F]
    downstream.context <- str_split_fixed(str_split_fixed(ins.2bp[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,2,drop=F]
    for (i in 1:nrow(ins.2bp)){
      
      ## look at repeat
      tmp.repeat.nts <- repeat.nts[i]
      tmp.repeat.length <- nchar(tmp.repeat.nts)
      
      ## how often does it match consecutively in the immediate downstream context?
      tmp.downstream.context <- downstream.context[i]
      tmp.downstream.context <- strsplit(tmp.downstream.context, '')[[1]]
      
      ## group
      tmp.downstream.context <- c(paste(tmp.downstream.context[c(tmp.repeat.length-c(tmp.repeat.length - 1)):tmp.repeat.length], collapse = ''),
                                  paste(tmp.downstream.context[c(2*tmp.repeat.length-c(tmp.repeat.length - 1)):c(2*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(3*tmp.repeat.length-c(tmp.repeat.length - 1)):c(3*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(4*tmp.repeat.length-c(tmp.repeat.length - 1)):c(4*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(5*tmp.repeat.length-c(tmp.repeat.length - 1)):c(5*tmp.repeat.length)], collapse = ''))
      
      ### check all 6 categories for upstream bases
      if(tmp.downstream.context[1] %in% tmp.repeat.nts){
        
        if(tmp.downstream.context[2] %in% tmp.repeat.nts){
          
          if(tmp.downstream.context[3] %in% tmp.repeat.nts){
            
            if(tmp.downstream.context[4] %in% tmp.repeat.nts){
              
              if(tmp.downstream.context[5] %in% tmp.repeat.nts){
                
                ins.greater.2bp.summary['5+', '2 bp'] <- ins.greater.2bp.summary['5+', '2 bp'] + 1
                
              } else{
                
                ins.greater.2bp.summary['4', '2 bp'] <- ins.greater.2bp.summary['4', '2 bp'] + 1
                
              }
              
            }else{
              
              ins.greater.2bp.summary['3', '2 bp'] <- ins.greater.2bp.summary['3', '2 bp'] + 1
              
            }
            
          }else{
            
            ins.greater.2bp.summary['2', '2 bp'] <- ins.greater.2bp.summary['2', '2 bp'] + 1
            
          }
          
        }else{
          
          ins.greater.2bp.summary['1', '2 bp'] <- ins.greater.2bp.summary['1', '2 bp'] + 1
          
        }
        
      }else{
        
        ins.greater.2bp.summary['0', '2 bp'] <- ins.greater.2bp.summary['0', '2 bp'] + 1
        
      }
      
    } 
    
  }
  
  ## 3 bp
  if(nrow(ins.3bp) > 0){
    
    repeat.nts <- str_split_fixed(str_split_fixed(ins.3bp[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,1,drop=F]
    downstream.context <- str_split_fixed(str_split_fixed(ins.3bp[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,2,drop=F]
    for (i in 1:nrow(ins.3bp)){
      
      ## look at repeat
      tmp.repeat.nts <- repeat.nts[i]
      tmp.repeat.length <- nchar(tmp.repeat.nts)
      
      ## how often does it match consecutively in the immediate downstream context?
      tmp.downstream.context <- downstream.context[i]
      tmp.downstream.context <- strsplit(tmp.downstream.context, '')[[1]]
      
      ## group
      tmp.downstream.context <- c(paste(tmp.downstream.context[c(tmp.repeat.length-c(tmp.repeat.length - 1)):tmp.repeat.length], collapse = ''),
                                  paste(tmp.downstream.context[c(2*tmp.repeat.length-c(tmp.repeat.length - 1)):c(2*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(3*tmp.repeat.length-c(tmp.repeat.length - 1)):c(3*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(4*tmp.repeat.length-c(tmp.repeat.length - 1)):c(4*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(5*tmp.repeat.length-c(tmp.repeat.length - 1)):c(5*tmp.repeat.length)], collapse = ''))
      
      ### check all 6 categories for upstream bases
      if(tmp.downstream.context[1] %in% tmp.repeat.nts){
        
        if(tmp.downstream.context[2] %in% tmp.repeat.nts){
          
          if(tmp.downstream.context[3] %in% tmp.repeat.nts){
            
            if(tmp.downstream.context[4] %in% tmp.repeat.nts){
              
              if(tmp.downstream.context[5] %in% tmp.repeat.nts){
                
                ins.greater.2bp.summary['5+', '3 bp'] <- ins.greater.2bp.summary['5+', '3 bp'] + 1
                
              } else{
                
                ins.greater.2bp.summary['4', '3 bp'] <- ins.greater.2bp.summary['4', '3 bp'] + 1
                
              }
              
            }else{
              
              ins.greater.2bp.summary['3', '3 bp'] <- ins.greater.2bp.summary['3', '3 bp'] + 1
              
            }
            
          }else{
            
            ins.greater.2bp.summary['2', '3 bp'] <- ins.greater.2bp.summary['2', '3 bp'] + 1
            
          }
          
        }else{
          
          ins.greater.2bp.summary['1', '3 bp'] <- ins.greater.2bp.summary['1', '3 bp'] + 1
          
        }
        
      }else{
        
        ins.greater.2bp.summary['0', '3 bp'] <- ins.greater.2bp.summary['0', '3 bp'] + 1
        
      }
      
    } 
    
  }
  
  ## 4 bp
  if(nrow(ins.4bp) > 0){
    
    repeat.nts <- str_split_fixed(str_split_fixed(ins.4bp[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,1,drop=F]
    downstream.context <- str_split_fixed(str_split_fixed(ins.4bp[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,2,drop=F]
    for (i in 1:nrow(ins.4bp)){
      
      ## look at repeat
      tmp.repeat.nts <- repeat.nts[i]
      tmp.repeat.length <- nchar(tmp.repeat.nts)
      
      ## how often does it match consecutively in the immediate downstream context?
      tmp.downstream.context <- downstream.context[i]
      tmp.downstream.context <- strsplit(tmp.downstream.context, '')[[1]]
      
      ## group
      tmp.downstream.context <- c(paste(tmp.downstream.context[c(tmp.repeat.length-c(tmp.repeat.length - 1)):tmp.repeat.length], collapse = ''),
                                  paste(tmp.downstream.context[c(2*tmp.repeat.length-c(tmp.repeat.length - 1)):c(2*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(3*tmp.repeat.length-c(tmp.repeat.length - 1)):c(3*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(4*tmp.repeat.length-c(tmp.repeat.length - 1)):c(4*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(5*tmp.repeat.length-c(tmp.repeat.length - 1)):c(5*tmp.repeat.length)], collapse = ''))
      
      ### check all 6 categories for upstream bases
      if(tmp.downstream.context[1] %in% tmp.repeat.nts){
        
        if(tmp.downstream.context[2] %in% tmp.repeat.nts){
          
          if(tmp.downstream.context[3] %in% tmp.repeat.nts){
            
            if(tmp.downstream.context[4] %in% tmp.repeat.nts){
              
              if(tmp.downstream.context[5] %in% tmp.repeat.nts){
                
                ins.greater.2bp.summary['5+', '4 bp'] <- ins.greater.2bp.summary['5+', '4 bp'] + 1
                
              } else{
                
                ins.greater.2bp.summary['4', '4 bp'] <- ins.greater.2bp.summary['4', '4 bp'] + 1
                
              }
              
            }else{
              
              ins.greater.2bp.summary['3', '4 bp'] <- ins.greater.2bp.summary['3', '4 bp'] + 1
              
            }
            
          }else{
            
            ins.greater.2bp.summary['2', '4 bp'] <- ins.greater.2bp.summary['2', '4 bp'] + 1
            
          }
          
        }else{
          
          ins.greater.2bp.summary['1', '4 bp'] <- ins.greater.2bp.summary['1', '4 bp'] + 1
          
        }
        
      }else{
        
        ins.greater.2bp.summary['0', '4 bp'] <- ins.greater.2bp.summary['0', '4 bp'] + 1
        
      }
      
    } 
    
  }
  
  ## 5+ bp
  if(nrow(ins.5bp) > 0){
    
    repeat.nts <- str_split_fixed(str_split_fixed(ins.5bp[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,1,drop=F]
    downstream.context <- str_split_fixed(str_split_fixed(ins.5bp[,'CONTEXT'], '\\[', 2)[,2,drop=F], '\\]', 2)[,2,drop=F]
    for (i in 1:nrow(ins.5bp)){
      
      ## look at repeat
      tmp.repeat.nts <- repeat.nts[i]
      tmp.repeat.length <- nchar(tmp.repeat.nts)
      
      ## how often does it match consecutively in the immediate downstream context?
      tmp.downstream.context <- downstream.context[i]
      tmp.downstream.context <- strsplit(tmp.downstream.context, '')[[1]]
      
      ## group
      tmp.downstream.context <- c(paste(tmp.downstream.context[c(tmp.repeat.length-c(tmp.repeat.length - 1)):tmp.repeat.length], collapse = ''),
                                  paste(tmp.downstream.context[c(2*tmp.repeat.length-c(tmp.repeat.length - 1)):c(2*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(3*tmp.repeat.length-c(tmp.repeat.length - 1)):c(3*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(4*tmp.repeat.length-c(tmp.repeat.length - 1)):c(4*tmp.repeat.length)], collapse = ''),
                                  paste(tmp.downstream.context[c(5*tmp.repeat.length-c(tmp.repeat.length - 1)):c(5*tmp.repeat.length)], collapse = ''))
      
      ### check all 6 categories for upstream bases
      if(tmp.downstream.context[1] %in% tmp.repeat.nts){
        
        if(tmp.downstream.context[2] %in% tmp.repeat.nts){
          
          if(tmp.downstream.context[3] %in% tmp.repeat.nts){
            
            if(tmp.downstream.context[4] %in% tmp.repeat.nts){
              
              if(tmp.downstream.context[5] %in% tmp.repeat.nts){
                
                ins.greater.2bp.summary['5+', '5+ bp'] <- ins.greater.2bp.summary['5+', '5+ bp'] + 1
                
              } else{
                
                ins.greater.2bp.summary['4', '5+ bp'] <- ins.greater.2bp.summary['4', '5+ bp'] + 1
                
              }
              
            }else{
              
              ins.greater.2bp.summary['3', '5+ bp'] <- ins.greater.2bp.summary['3', '5+ bp'] + 1
              
            }
            
          }else{
            
            ins.greater.2bp.summary['2', '5+ bp'] <- ins.greater.2bp.summary['2', '5+ bp'] + 1
            
          }
          
        }else{
          
          ins.greater.2bp.summary['1', '5+ bp'] <- ins.greater.2bp.summary['1', '5+ bp'] + 1
          
        }
        
      }else{
        
        ins.greater.2bp.summary['0', '5+ bp'] <- ins.greater.2bp.summary['0', '5+ bp'] + 1
        
      }
      
    } 
    
  }
  
  ## 5. summarise outputs as five lists, each featuring one table per major ID category
  out <- list('1bp del' = dels.1bp.summary,
              '1bp ins' = ins.1bp.summary,
              '>2bp del' = dels.greater.2bp.summary,
              '>2bp ins' = ins.greater.2bp.summary,
              'MH dels' = dels.greater.2bp.MH.summary)
  return(out)
  
}
plot.indel.spectrum <- function(x, title, colors.bars){
  
  # a. Convert spectrum from list to vector
  x.vec <- c(x$`1bp del`[,1], x$`1bp del`[,2],
             x$`1bp ins`[,1], x$`1bp ins`[,2],
             x$`>2bp del`[,1], x$`>2bp del`[,2], x$`>2bp del`[,3], x$`>2bp del`[,4],
             x$`>2bp ins`[,1], x$`>2bp ins`[,2], x$`>2bp ins`[,3], x$`>2bp ins`[,4],
             x$`MH dels`[1,1], x$`MH dels`[1:2,2], x$`MH dels`[1:3,3], x$`MH dels`[1:5,4])
  
  # b. setup plot background colours and blocks
  mar.default <- c(2,4,2,2) + 0.1
  par(mar = mar.default + c(-2, 11, -2, -2))
  
  ## colors
  colors <- rep(c("grey40", "grey60"), 8)
  alphas <- rep(0.2,16)
  mut <- c("1 DEL", "1 INS", ">1 DEL", ">1 INS", "MH")
  mut2 <- c('C', 'T', 'C', 'T',
            '2','3', '4','5+',
            '2','3', '4','5+',
            '2-5+')
  
  y.top = 0.55
  
  borders <- c(0,
               100*12/83,
               100*24/83,
               100*48/83,
               100*72/83,
               100)
  
  borders2 = c(0,
               100*6/83, 
               100*12/83, 
               100*18/83, 
               100*24/83, 
               100*30/83, 
               100*36/83, 
               100*42/83, 
               100*48/83, 
               100*54/83, 
               100*60/83, 
               100*66/83, 
               100*72/83, 
               100)
  
  plot(1, type="n", ylim=c(0, y.top), xlim=c(0, 100), xlab="", ylab="", axes=F)
  for (i in 1:5) {
    rect(xleft=borders[i], xright=borders[i+1], ybottom=y.top*0.88, ytop=y.top, 
         col=alpha(colors[1], c(0.8,0.7,0.8,0.7,0.8)[i]), border="white")
    text(x = (borders[i]+borders[i+1])/2, 
         y = y.top*0.94, 
         labels = mut[i], 
         cex = 4.5,
         col = 'white')
  }
  
  for (i in 1:16) {
    rect(xleft=borders2[i], xright=borders2[i+1], ybottom=0, ytop=y.top*0.88, col=alpha(colors[i], alphas[i]), border="white")
    rect(xleft=borders2[i], xright=borders2[i+1], ybottom=y.top*0.76, ytop=y.top*0.88, col=colors[i], border="white")
    text(x=(borders2[i]+borders2[i+1])/2, 
         y=y.top*0.82, 
         labels=mut2[i], 
         cex=5, 
         col='white')
  }
  
  # c. add spectrum in %
  perc <- x.vec/sum(x.vec)
  out <- barplot(perc, 
                 col = colors.bars,
                 cex.names = 0.4,
                 names = rep('', 83),
                 las = 2,
                 border = NA,
                 main = paste0(title, " (", sum(x.vec), " Indels)"), 
                 cex.main = 2,
                 ylab = "", 
                 yaxt = 'n', 
                 ylim = c(0, 0.4),
                 add = T)
  
  # Axes
  axis(side = 2, las = 2, cex.axis = 4,
       at = c(0, 0.1, 0.2, 0.3, 0.4), labels = c('0', '10', '20', '30', '40'), lwd = 5, pos=-2, hadj = 1.4)
  title(ylab = "Indels [%]", line = 9, cex.lab = 6)
}

## DFT1-377T1 only SNV spectrum
DFT1.SNVs.377Tunique <- DFT1.SNVs[which(rowSums(DFT1.SNVs.bin) == 1 & DFT1.SNVs.bin[,'377T1'] == 1),]
DFT1.SNVs.377Tunique.spectrum <- substitution.spectrum(DFT1.SNVs.377Tunique, normalised = T)
pdf("Figure3H_377T1_SNV_spectrum.pdf",
    height = 10, width = 24)
plot.substition.spectrum(DFT1.SNVs.377Tunique.spectrum$counts.normalised,
                  title = '', 
                  peak.colour = 'cornflowerblue')
dev.off()

## DFT1-377T1 only Indel spectrum
DFT1.Indels.377Tunique <- DFT1.Indels[which(rowSums(DFT1.Indels.bin) == 1 & DFT1.Indels.bin[,'377T1'] == 1),]
DFT1.Indels.377Tunique.spectrum <- indel.spectrum(DFT1.Indels.377Tunique)
pdf("Figure3H_377T1_Indel_spectrum.pdf", 
    height = 10, width = 24)
plot.indel.spectrum(DFT1.Indels.377Tunique.spectrum,
                    title = '', 
                    colors.bars = 'cornflowerblue')
dev.off()


# MLH1 deletion - Campbellgram #
################################

### load and subset 377T1 SVs
DFT1.SVs <- as.matrix(read_xlsx('/Users/mstammnitz/Desktop/DFT_evolution/doc/manuscripts/The Evolutionary History of Two Transmissible Cancers in Tasmanian Devils/Tables/v6/Table-S5_v6.xlsx', sheet = 2))
colnames(DFT1.SVs) <- as.character(DFT1.SVs[2,])
DFT1.SVs <- DFT1.SVs[-c(1:2),]
DFT1.SVs.377T1.MSG <- DFT1.SVs[grep('MSG', DFT1.SVs[,'CALLER']),,drop=F]
DFT1.SVs.377T1.MSG <- DFT1.SVs.377T1.MSG[which(apply(DFT1.SVs.377T1.MSG[,grep('377T1', colnames(DFT1.SVs.377T1.MSG)),drop=F], 1, function(x){out <- as.numeric(str_split_fixed(x, '/', 2)[,1]); all(out >= 1)}) == T),,drop=F]
DFT1.SVs.377T1.SvABA <- DFT1.SVs[-grep('MSG', DFT1.SVs[,'CALLER']),,drop=F]
DFT1.SVs.377T1.SvABA <- DFT1.SVs.377T1.SvABA[which(apply(DFT1.SVs.377T1.SvABA[,grep('377T1', colnames(DFT1.SVs.377T1.SvABA)),drop=F], 1, function(x){out <- as.numeric(x); all(out >= 3)}) == T),,drop=F]

### 377T1 SVs (all non-unique)
DFT1.SVs.377T1.MSG.non.unique <- DFT1.SVs.377T1.MSG[which(apply(DFT1.SVs.377T1.MSG[,c(42:105)], 1, function(x){out <- as.numeric(str_split_fixed(x, '/', 2)[,1]); any(out >= 1)}) == T),,drop=F]
DFT1.SVs.377T1.SvABA.non.unique <- DFT1.SVs.377T1.SvABA[which(apply(DFT1.SVs.377T1.SvABA[,c(42:105)], 1, function(x){out <- as.numeric(x); any(out >= 3)}) == T),,drop=F]
DFT1.SVs.377T1.non.unique <- rbind(DFT1.SVs.377T1.MSG.non.unique, DFT1.SVs.377T1.SvABA.non.unique)
DFT1.SVs.377T1.non.unique <- DFT1.SVs.377T1.non.unique[order(as.character(DFT1.SVs.377T1.non.unique[,"CHROM1"]), as.numeric(DFT1.SVs.377T1.non.unique[,"ORIGINAL POS1"])),,drop=F]
DFT1.SVs.377T1.non.unique.left.GR <- GRanges(seqnames = as.character(DFT1.SVs.377T1.non.unique[,'CHROM1']), 
                                             ranges = IRanges::IRanges(start = as.numeric(as.character(DFT1.SVs.377T1.non.unique[,'ORIGINAL POS1'])),
                                                                       end = as.numeric(as.character(DFT1.SVs.377T1.non.unique[,'ORIGINAL POS1']))))
DFT1.SVs.377T1.non.unique.right.GR <- GRanges(seqnames = as.character(DFT1.SVs.377T1.non.unique[,'CHROM2']), 
                                              ranges = IRanges::IRanges(start = as.numeric(as.character(DFT1.SVs.377T1.non.unique[,'ORIGINAL POS2'])),
                                                                        end = as.numeric(as.character(DFT1.SVs.377T1.non.unique[,'ORIGINAL POS2']))))
OL.left <- findOverlaps(DFT1.SVs.377T1.non.unique.left.GR, GRanges(seqnames = 1, 
                                                                   ranges = IRanges::IRanges(start = 1,
                                                                                             end = 50000000)))
OL.right <- findOverlaps(DFT1.SVs.377T1.non.unique.right.GR, GRanges(seqnames = 1, 
                                                                     ranges = IRanges::IRanges(start = 1,
                                                                                               end = 50000000)))
OL.left <- as.matrix(OL.left); OL.right <- as.matrix(OL.right)
DFT1.SVs.377T1.non.unique.local <- DFT1.SVs.377T1.non.unique[OL.left[OL.left[,1] %in% OL.right[,1],1],,drop=F]
DFT1.SVs.377T1.non.unique.inter.right <- DFT1.SVs.377T1.non.unique[OL.left[!OL.left[,1] %in% OL.right[,1],1],,drop=F]
DFT1.SVs.377T1.non.unique.inter.left <- DFT1.SVs.377T1.non.unique[OL.right[!OL.right[,1] %in% OL.left[,1],1],,drop=F]

### 377T1 SVs (all unique)
DFT1.SVs.377T1.MSG.unique <- DFT1.SVs.377T1.MSG[which(apply(DFT1.SVs.377T1.MSG[,c(42:105)], 1, function(x){out <- as.numeric(str_split_fixed(x, '/', 2)[,1]); all(out < 1)}) == T),,drop=F]
DFT1.SVs.377T1.SvABA.unique <- DFT1.SVs.377T1.SvABA[which(apply(DFT1.SVs.377T1.SvABA[,c(42:105)], 1, function(x){out <- as.numeric(x); all(out < 3)}) == T),,drop=F]
DFT1.SVs.377T1.unique <- rbind(DFT1.SVs.377T1.MSG.unique, DFT1.SVs.377T1.SvABA.unique)
DFT1.SVs.377T1.unique <- DFT1.SVs.377T1.unique[order(as.character(DFT1.SVs.377T1.unique[,"CHROM1"]), as.numeric(DFT1.SVs.377T1.unique[,"ORIGINAL POS1"])),,drop=F]
DFT1.SVs.377T1.unique.left.GR <- GRanges(seqnames = as.character(DFT1.SVs.377T1.unique[,'CHROM1']), 
                                         ranges = IRanges::IRanges(start = as.numeric(as.character(DFT1.SVs.377T1.unique[,'ORIGINAL POS1'])),
                                                                   end = as.numeric(as.character(DFT1.SVs.377T1.unique[,'ORIGINAL POS1']))))
DFT1.SVs.377T1.unique.right.GR <- GRanges(seqnames = as.character(DFT1.SVs.377T1.unique[,'CHROM2']), 
                                          ranges = IRanges::IRanges(start = as.numeric(as.character(DFT1.SVs.377T1.unique[,'ORIGINAL POS2'])),
                                                                    end = as.numeric(as.character(DFT1.SVs.377T1.unique[,'ORIGINAL POS2']))))
OL.left <- findOverlaps(DFT1.SVs.377T1.unique.left.GR, GRanges(seqnames = 1, 
                                                               ranges = IRanges::IRanges(start = 1,
                                                                                         end = 50000000)))
OL.right <- findOverlaps(DFT1.SVs.377T1.unique.right.GR, GRanges(seqnames = 1, 
                                                                 ranges = IRanges::IRanges(start = 1,
                                                                                           end = 50000000)))
OL.left <- as.matrix(OL.left); OL.right <- as.matrix(OL.right)
DFT1.SVs.377T1.unique.local <- DFT1.SVs.377T1.unique[OL.left[OL.left[,1] %in% OL.right[,1],1],,drop=F]
DFT1.SVs.377T1.unique.inter.right <- DFT1.SVs.377T1.unique[OL.left[!OL.left[,1] %in% OL.right[,1],1],,drop=F]
DFT1.SVs.377T1.unique.inter.left <- DFT1.SVs.377T1.unique[OL.right[!OL.right[,1] %in% OL.left[,1],1],,drop=F]

## load 377T1 (normalised LogR or CNV calls from Kevin) 
copynumber.377T1 <- read.table("/Users/mstammnitz/Desktop/DFT_evolution/doc/manuscripts/The Evolutionary History of Two Transmissible Cancers in Tasmanian Devils/Tables/v6/Supplementary_data/Copynumber_377T1_chr1.txt", header = T)
copynumber.377T1 <- copynumber.377T1[which(copynumber.377T1[,"excluded_from_segmentation"] == F),]

## display CN estimates for SNPs & SNVs
ylims <- c(0.5, 3.5)
lwds <- 2
transp <- 0.6

pdf("Figure3H_Campbellgram_MLH1.pdf", 
    height = 16, width = 13)
par(mar = c(15, 12, 2, 2))
plot(x = c(0,50000000),
     y = c(1,1),
     pch = 16, 
     col = 'white',
     cex = 1,
     ylab = '', 
     xlab = '',
     yaxt = 'n', 
     xaxt = 'n',
     bty = 'n',
     xlim = c(0,50000000),
     ylim = ylims)

### local non-unique SVs
fit.func <- function(x){
  y <- fit$coefficients[3]*x^2+c(fit$coefficients[2]*x)+fit$coefficients[1]
  return(y)
}

### distant non-unique SVs
for (i in 1:nrow(DFT1.SVs.377T1.non.unique.inter.right)){
  
  ## arcs
  bp <- DFT1.SVs.377T1.non.unique.inter.right[i,]
  bp <- bp[c(grep('ORIGINAL POS1', colnames(DFT1.SVs.377T1.non.unique.inter.right)),
             grep('CHR2', colnames(DFT1.SVs.377T1.non.unique.inter.right)))]
  end <- as.numeric(bp[1]) + 10*c(50000000 - as.numeric(bp[1]))
  
  if(bp[2] != 'X' & bp[2] != '6' & bp[2] != '5'){
    
    bp.x <- c(as.numeric(bp[1]), 50000000, 50000000 + c(50000000 - as.numeric(bp[1])))
    
  }else{
    
    bp.x <- c(0 - c(as.numeric(bp[1]) - 0), 0, as.numeric(bp[1]))
    
  }
  
  fit <- glm(c(3, 3.5, 3) ~ poly(bp.x, 2, raw=TRUE))
  x.fit.out <- seq(f = bp.x[1], t = bp.x[3], length.out = 5000)
  y.fit.out <- fit.func(seq(f = bp.x[1], t = bp.x[3], length.out = 5000))
  lines(x = x.fit.out[which(x.fit.out > 0 & x.fit.out < 50000000)],
        y = y.fit.out[which(x.fit.out > 0 & x.fit.out < 50000000)],
        col = 'black',
        pch = 16,
        lwd = lwds)
  
  ## borders
  lines(x = rep(bp[1],2), y = c(0.5,3), col = 'black', lwd = lwds)
  
}

### local unique SVs
fit.func <- function(x){
  y <- fit$coefficients[3]*x^2+c(fit$coefficients[2]*x)+fit$coefficients[1]
  return(y)
}
for (i in 1:nrow(DFT1.SVs.377T1.unique.local)){
  
  ## arcs
  bp <- DFT1.SVs.377T1.unique.local[i,]
  bp <- bp[grep('ORIGINAL', colnames(DFT1.SVs.377T1.unique.local))]
  bp.x <- c(as.numeric(bp[1]),c(as.numeric(bp[1])+as.numeric(bp[2]))/2,as.numeric(bp[2]))
  bp.y <- c(3, 3.5, 3)
  fit <- glm(bp.y ~ poly(bp.x, 2, raw = T))
  x.fit.out <- seq(f = bp.x[1], t = bp.x[3], length.out = 5000)
  y.fit.out <- fit.func(x.fit.out)
  lines(x = x.fit.out[which(x.fit.out > 0 & x.fit.out < 50000000)],
        y = y.fit.out[which(x.fit.out > 0 & x.fit.out < 50000000)],
        col = 'cornflowerblue', 
        pch = 16, 
        lwd = lwds)
  
  ## borders
  lines(x = rep(bp[1],2), y = c(0.5,3), col = 'cornflowerblue', lwd = lwds)
  lines(x = rep(bp[2],2), y = c(0.5,3), col = 'cornflowerblue', lwd = lwds)
}

### data points (besides MLH1)
points(x = as.numeric(copynumber.377T1[which(as.numeric(copynumber.377T1[,'START']) < 13952508 | as.numeric(copynumber.377T1[,'END']) > 14001592),'START']),
       y = as.numeric(copynumber.377T1[which(as.numeric(copynumber.377T1[,'START']) < 13952508 | as.numeric(copynumber.377T1[,'END']) > 14001592),'total_cn']),
       pch = 16, 
       col = alpha('cornflowerblue', transp),
       cex = 0.6)

### data points (only MLH1)
points(x = as.numeric(copynumber.377T1[which(as.numeric(copynumber.377T1[,'START']) >= 13952508 & as.numeric(copynumber.377T1[,'END']) < 14001592),'START']),
       y = as.numeric(copynumber.377T1[which(as.numeric(copynumber.377T1[,'START']) >= 13952508 & as.numeric(copynumber.377T1[,'END']) < 14001592),'total_cn']),
       pch = 16, 
       col = alpha('black', transp),
       cex = 1.8)

### axes
axis(1, at = seq(from = 0, to = 50000000, length = 6), 
     labels = paste0(seq(from = 0, to = 50000000, length = 6)/1000000),
     cex.axis = 5, 
     padj = 1, 
     lwd = 5)
mtext(side = 1, 
      text = 'Chromosome 1 [MBP]', 
      line = 12, 
      cex = 6)

axis(2, 
     at = c(1, 2, 3), 
     labels = c(1, 2, 3), 
     cex.axis = 5, 
     las = 2, 
     lwd
     = 5, 
     line = -1)
mtext(side = 2, 
      text = 'Copy number', 
      line = 7, 
      cex = 6)

dev.off()