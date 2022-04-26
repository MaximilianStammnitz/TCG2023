## Figure 3G
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

## set input path(s)
setwd('/Tables')


# DFT1 tree #
#############

## clade C2/3 separation
clades <- as.matrix(read_xlsx('Table-S2.xlsx', sheet = 1))
colnames(clades) <- clades[2,]
clades.dft1 <- clades[grep('DFT1', clades[,9]),]
samples.cladeC <-  clades.dft1[grep('DFT1-C', clades.dft1[,10]),]
samples.cladeC <- samples.cladeC[-grep('DFT1-C1', samples.cladeC[,10]),]
samples.cladeC <- samples.cladeC[-grep('1439T7', samples.cladeC[,8]),]

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

pdf("Figure3G_DFT1_C2-3_tree.pdf",
    width = 13, height = 11)

## hack in the edge colors and width
edge.sizes <- rep(1, DFT1.beast.tree.hpd@phylo$Nnode+length(DFT1.beast.tree.hpd@phylo$tip.label))

## highlight branch
edge.sizes[clade_cnot1] <- 4
d <- data.frame(node=1:(DFT1.beast.tree.hpd@phylo$Nnode+length(DFT1.beast.tree.hpd@phylo$tip.label)), 
                size = edge.sizes)

p <- ggtree(DFT1.beast.tree.hpd, 
            root.position = 2018.13424657534 - max(DFT1.beast.tree.hpd@data[,'height']), 
            layout = 'rectangular', 
            col = 'black',
            ladderize = T) %<+% d + aes(lwd = I(size)) + 
  scale_x_continuous(breaks = seq(f=1990, t=2020, by=10),
                     labels = seq(f=1990, t=2020, by=10),
                     limits = c(1986,2020)) +
  geom_rootedge(rootedge = 2018.13424657534 - max(DFT1.beast.tree.hpd@data[,'height']) - mean(DFT1.beast.trees.MRCAs.rates[,'age at 818 singletons']),
                lwd = 1, col = 'black') +
  geom_point2(aes(subset = (node %in% 79),
                  x = x - c(2018.13424657534 - max(DFT1.beast.tree.hpd@data[,'height']) - mean(DFT1.beast.trees.MRCAs.rates[,'age at 818 singletons']))),
              shape = 16, size = 3, col = 'black') +
  geom_hilight(node = clade_cnot1, fill = "cornflowerblue", alpha = 0.3, extend = 0.5) +
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

## clean up environment
rm(list=ls())


# DFT1-C2/3 sampling map #
##########################

## clade C2/3 separation
clades <- as.matrix(read_xlsx('Table-S2.xlsx', sheet = 1))
colnames(clades) <- clades[2,]
clades.dft1 <- clades[grep('DFT1', clades[,9]),]
samples.cladeC <-  clades.dft1[grep('DFT1-C', clades.dft1[,10]),]
samples.cladeC <- samples.cladeC[-grep('DFT1-C1', samples.cladeC[,10]),]
samples.cladeC[23,"LOCATION OF SAMPLING"] <- str_split_fixed(samples.cladeC[23,"LOCATION OF SAMPLING"], " ", 2)[,1]
set.seed(100)
add.noise <- names(table(samples.cladeC[,"LOCATION OF SAMPLING"])[table(samples.cladeC[,"LOCATION OF SAMPLING"]) > 1])
mult <- samples.cladeC[samples.cladeC[,"LOCATION OF SAMPLING"] %in% add.noise,]
sing <- samples.cladeC[!samples.cladeC[,"LOCATION OF SAMPLING"] %in% add.noise,]
for (i in 1:length(add.noise)){
  print(i)
  add.noise.tmp <- add.noise[i]
  add.noise.tmp.id <- which(mult[,"LOCATION OF SAMPLING"] == add.noise.tmp)
  mult[add.noise.tmp.id,"LOCATION LATITUDE"] <- as.numeric(mult[add.noise.tmp.id,"LOCATION LATITUDE"]) + runif(n = length(add.noise.tmp.id), min = -0.1, max = 0.1)
  mult[add.noise.tmp.id,"LOCATION LONGITUDE"] <- as.numeric(mult[add.noise.tmp.id,"LOCATION LONGITUDE"]) + runif(n = length(add.noise.tmp.id), min = -0.1, max = 0.1)
}
summary <- rbind(mult,sing)
rownames(summary) <- summary[,8]
summary.sf <- summary[,6:7]
class(summary.sf) <- 'numeric'
summary.sf <- as.data.frame(summary.sf)

## downlaod Tasmania map via ggmap package, requires a Google API key:
## (example video: https://www.youtube.com/watch?v=ewYGC2JjKuE)
## register_google(key = '') fill in
tasmania.map.googlemaps_terrain <- get_googlemap(center = c(146.8087, -42.0409), 
                                                 zoom = 7,
                                                 size = c(640, 640),
                                                 scale = 2,
                                                 maptype = 'terrain',
                                                 color = 'bw',
                                                 style = 'style=feature:administrative.country|element:labels|visibility:off&style=feature:water|color:0xffffff|visibility:simplified')

## plot
png('Figure3G_DFT1_C2-3_map.png', 
    width = 1100, height = 1100)
ggmap(tasmania.map.googlemaps_terrain,
      extent = 'panel',
      padding = 0)  +
  scale_x_continuous(limits = c(144.3, 148.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-43.7, -40.55), expand = c(0, 0)) +
  geom_point(mapping = aes(x = Longitude, y = Latitude, group = NULL), 
             data = summary.sf,
             size = 30,
             shape = 'triangle',
             show.legend = F,
             colour = 'cornflowerblue',
             alpha = 0.6) +
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

## clean up environment
rm(list=ls())


# DFT1-C2/3 substitution rates #
################################

## import data
counts <- as.matrix(read_xlsx('Table-S3.xlsx', sheet = 3))
colnames(counts) <- as.character(counts[2,])
counts <- counts[-c(1:2),]
dft1.counts <- counts[counts[,'LINEAGE'] == 'DFT1',]
dft1.counts <- dft1.counts[which(dft1.counts[,'TUMOUR ID'] != '377T1'),]

## metadata
samples <- as.matrix(read_xlsx('Table-S2.xlsx', sheet = 1))
dft1.samples <- samples[samples[,9] == 'DFT1',]
dft1.samples <- dft1.samples[which(is.na(dft1.samples[,9]) == F),]
dft1.samples <- dft1.samples[,c(8,4,13)]
dft1.samples[,1] <- str_split_fixed(dft1.samples[,1], ' ', 2)[,1]
dft1.samples[,2] <- gsub('[*]', '', dft1.samples[,2])
dft1.samples <- dft1.samples[which(dft1.samples[,1] != '377T1'),]
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

## clade C2/3 separation
clades <- as.matrix(read_xlsx('Table-S2.xlsx', sheet = 1))
colnames(clades) <- clades[2,]
clades.dft1 <- clades[grep('DFT1', clades[,9]),]
samples.cladeC <-  clades.dft1[grep('DFT1-C', clades.dft1[,10]),]
samples.cladeC <- samples.cladeC[-grep('DFT1-C1', samples.cladeC[,10]),]
samples.cladeC <- samples.cladeC[-grep('1439T7', samples.cladeC[,8]),]

pdf("Figure3G_substitution_rates.pdf", 
    height = 12, width = 18)
ggplot(DFT1.SNVs.counts, aes(x = `Collection Date`, y = `SNVs`)) + 
  scale_x_continuous(breaks = seq(from = 1990, to = 2030, by = 10), limits = c(1980, 2040)) + 
  scale_y_continuous(breaks = seq(from = 0, to = 10000, by = 2000), limits = c(-15000, 15000)) + 
  coord_cartesian(xlim = c(1985, 2035), ylim = c(0, 10000), expand = F) +
  geom_point(data = DFT1.SNVs.counts[-match(samples.cladeC[,8], rownames(DFT1.SNVs.counts)),], 
             color = "cornflowerblue",
             size = 5, shape = 'circle') + 
  geom_smooth(data = DFT1.SNVs.counts[-match(samples.cladeC[,8], rownames(DFT1.SNVs.counts)),],
              method = 'lm', color = "cornflowerblue", fullrange = T) +
  geom_point(data = DFT1.SNVs.counts[match(samples.cladeC[,8], rownames(DFT1.SNVs.counts)),], 
             color = alpha("cornflowerblue", 0.6), 
             size = 13, shape = 'triangle') + 
  geom_smooth(data = DFT1.SNVs.counts[match(samples.cladeC[,8], rownames(DFT1.SNVs.counts)),],
              method = 'lm', color = "cornflowerblue", fullrange = T, fill = alpha('cornflowerblue', 0.3)) +
  labs(y = "Substitutions", x = "Year") +
  theme_classic(base_size = 50) +
  theme(axis.text = element_text(size = 50),
        axis.line.x = element_line(size = 1, color = 'black'),
        axis.line.y = element_line(size = 1, color = 'black'),
        axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = 70, vjust = -1),
        axis.title.y = element_text(family = 'Helvetica', colour = 'black', size = 70, vjust = 3),
        legend.position = "none",
        text = element_text(family="Helvetica"),
        plot.margin = unit(c(2, 2, 2, 2),"cm"))
dev.off()

## clean up environment
rm(list=ls())


# DFT1-C2/3 SBS1 rates #
########################

## import data
counts <- as.matrix(read_xlsx('Table-S3_v6.xlsx', sheet = 3))
colnames(counts) <- as.character(counts[2,])
counts <- counts[-c(1:2),]
dft1.counts <- counts[counts[,'LINEAGE'] == 'DFT1',]
dft1.counts <- dft1.counts[which(dft1.counts[,'TUMOUR ID'] != '377T1'),]

## metadata
samples <- as.matrix(read_xlsx('Table-S2.xlsx', sheet = 1))
dft1.samples <- samples[samples[,9] == 'DFT1',]
dft1.samples <- dft1.samples[which(is.na(dft1.samples[,9]) == F),]
dft1.samples <- dft1.samples[,c(8,4,13)]
dft1.samples[,1] <- str_split_fixed(dft1.samples[,1], ' ', 2)[,1]
dft1.samples[,2] <- gsub('[*]', '', dft1.samples[,2])
dft1.samples <- dft1.samples[which(dft1.samples[,1] != '377T1'),]
colnames(dft1.samples) <- c('Tumour', 'Collection Date', 'Purity')

## combine
DFT1.SNVs.counts.SBS1 <- dft1.samples
DFT1.SNVs.counts.SBS1 <- cbind(DFT1.SNVs.counts.SBS1[match(dft1.counts[,'TUMOUR ID'], DFT1.SNVs.counts.SBS1[,1]),],
                          dft1.counts[,"WGD-NORMALISED SBS1 SUBSTITUTIONS"])
colnames(DFT1.SNVs.counts.SBS1)[4] <- 'SNVs'
rownames(DFT1.SNVs.counts.SBS1) <- DFT1.SNVs.counts.SBS1[,1]
DFT1.SNVs.counts.SBS1 <- DFT1.SNVs.counts.SBS1[,-1]
DFT1.SNVs.counts.SBS1[,'Collection Date'] <- decimal_date(as.Date(as.character(unlist(DFT1.SNVs.counts.SBS1[,'Collection Date'])), format="%d.%m.%Y"))
DFT1.SNVs.counts.SBS1 <- DFT1.SNVs.counts.SBS1[!is.na(DFT1.SNVs.counts.SBS1[,'Collection Date']),]
DFT1.SNVs.counts.SBS1 <- as.data.frame(DFT1.SNVs.counts.SBS1)
DFT1.SNVs.counts.SBS1[,'Collection Date'] <- as.numeric(as.character(DFT1.SNVs.counts.SBS1[,'Collection Date']))
DFT1.SNVs.counts.SBS1[,'Purity'] <- as.numeric(as.character(DFT1.SNVs.counts.SBS1[,'Purity']))
DFT1.SNVs.counts.SBS1[,'SNVs'] <- as.numeric(as.character(DFT1.SNVs.counts.SBS1[,'SNVs']))

## clade C2/3 separation
clades <- as.matrix(read_xlsx('Table-S2.xlsx', sheet = 1))
colnames(clades) <- clades[2,]
clades.dft1 <- clades[grep('DFT1', clades[,9]),]
samples.cladeC <-  clades.dft1[grep('DFT1-C', clades.dft1[,10]),]
samples.cladeC <- samples.cladeC[-grep('DFT1-C1', samples.cladeC[,10]),]
samples.cladeC <- samples.cladeC[-grep('1439T7', samples.cladeC[,8]),]

pdf("Figure3G_SBS1_rates.pdf", 
    height = 9, width = 9.3)
ggplot(DFT1.SNVs.counts.SBS1, aes(x = `Collection Date`, y = `SNVs`)) + 
  scale_x_continuous(breaks = seq(from = 2000, to = 2020, by = 10), limits = c(1985, 2035)) + 
  scale_y_continuous(breaks = seq(from = 0, to = 750, by = 150), limits = c(-1000, 1000)) + 
  coord_cartesian(xlim = c(1995, 2025), ylim = c(0, 750), expand = F) +
  geom_point(data = DFT1.SNVs.counts.SBS1[-match(samples.cladeC[,8], rownames(DFT1.SNVs.counts.SBS1)),], 
             color = "cornflowerblue",
             size = 3, shape = 'circle') + 
  geom_smooth(data = DFT1.SNVs.counts.SBS1[-match(samples.cladeC[,8], rownames(DFT1.SNVs.counts.SBS1)),],
              method = 'lm', color = "cornflowerblue", fullrange = T) +
  geom_point(data = DFT1.SNVs.counts.SBS1[match(samples.cladeC[,8], rownames(DFT1.SNVs.counts.SBS1)),], 
             color = alpha("cornflowerblue", 0.6), 
             size = 9, shape = 'triangle') + 
  geom_smooth(data = DFT1.SNVs.counts.SBS1[match(samples.cladeC[,8], rownames(DFT1.SNVs.counts.SBS1)),],
              method = 'lm', color = "cornflowerblue", fullrange = T, fill = alpha('cornflowerblue', 0.3)) +
  labs(y = "SBS1 Substitutions", x = "Year") +
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

## clean up environment
rm(list=ls())


# DFT1-C2/3 SBS5 rates #
########################

## import data
counts <- as.matrix(read_xlsx('Table-S3.xlsx', sheet = 3))
colnames(counts) <- as.character(counts[2,])
counts <- counts[-c(1:2),]
dft1.counts <- counts[counts[,'LINEAGE'] == 'DFT1',]
dft1.counts <- dft1.counts[which(dft1.counts[,'TUMOUR ID'] != '377T1'),]

## metadata
samples <- as.matrix(read_xlsx('Table-S2.xlsx', sheet = 1))
dft1.samples <- samples[samples[,9] == 'DFT1',]
dft1.samples <- dft1.samples[which(is.na(dft1.samples[,9]) == F),]
dft1.samples <- dft1.samples[,c(8,4,13)]
dft1.samples[,1] <- str_split_fixed(dft1.samples[,1], ' ', 2)[,1]
dft1.samples[,2] <- gsub('[*]', '', dft1.samples[,2])
dft1.samples <- dft1.samples[which(dft1.samples[,1] != '377T1'),]
colnames(dft1.samples) <- c('Tumour', 'Collection Date', 'Purity')

## combine
DFT1.SNVs.counts.SBS5 <- dft1.samples
DFT1.SNVs.counts.SBS5 <- cbind(DFT1.SNVs.counts.SBS5[match(dft1.counts[,'TUMOUR ID'], DFT1.SNVs.counts.SBS5[,1]),],
                               dft1.counts[,"WGD-NORMALISED SBS5 SUBSTITUTIONS"])
colnames(DFT1.SNVs.counts.SBS5)[4] <- 'SNVs'
rownames(DFT1.SNVs.counts.SBS5) <- DFT1.SNVs.counts.SBS5[,1]
DFT1.SNVs.counts.SBS5 <- DFT1.SNVs.counts.SBS5[,-1]
DFT1.SNVs.counts.SBS5[,'Collection Date'] <- decimal_date(as.Date(as.character(unlist(DFT1.SNVs.counts.SBS5[,'Collection Date'])), format="%d.%m.%Y"))
DFT1.SNVs.counts.SBS5 <- DFT1.SNVs.counts.SBS5[!is.na(DFT1.SNVs.counts.SBS5[,'Collection Date']),]
DFT1.SNVs.counts.SBS5 <- as.data.frame(DFT1.SNVs.counts.SBS5)
DFT1.SNVs.counts.SBS5[,'Collection Date'] <- as.numeric(as.character(DFT1.SNVs.counts.SBS5[,'Collection Date']))
DFT1.SNVs.counts.SBS5[,'Purity'] <- as.numeric(as.character(DFT1.SNVs.counts.SBS5[,'Purity']))
DFT1.SNVs.counts.SBS5[,'SNVs'] <- as.numeric(as.character(DFT1.SNVs.counts.SBS5[,'SNVs']))

## clade C2/3 separation
clades <- as.matrix(read_xlsx('Table-S2.xlsx', sheet = 1))
colnames(clades) <- clades[2,]
clades.dft1 <- clades[grep('DFT1', clades[,9]),]
samples.cladeC <-  clades.dft1[grep('DFT1-C', clades.dft1[,10]),]
samples.cladeC <- samples.cladeC[-grep('DFT1-C1', samples.cladeC[,10]),]
samples.cladeC <- samples.cladeC[-grep('1439T7', samples.cladeC[,8]),]

pdf("Figure3G_SBS5_rates.pdf", 
    height = 9, width = 9.3)
ggplot(DFT1.SNVs.counts.SBS5, aes(x = `Collection Date`, y = `SNVs`)) + 
  scale_x_continuous(breaks = seq(from = 2000, to = 2020, by = 10), limits = c(1985, 2035)) + 
  scale_y_continuous(breaks = seq(from = 0, to = 10000, by = 2000), limits = c(-10000, 15000)) + 
  coord_cartesian(xlim = c(1995, 2025), ylim = c(0, 10000), expand = F) +
  geom_point(data = DFT1.SNVs.counts.SBS5[-match(samples.cladeC[,8], rownames(DFT1.SNVs.counts.SBS5)),], 
             color = "cornflowerblue",
             size = 3, shape = 'circle') + 
  geom_smooth(data = DFT1.SNVs.counts.SBS5[-match(samples.cladeC[,8], rownames(DFT1.SNVs.counts.SBS5)),],
              method = 'lm', color = "cornflowerblue", fullrange = T) +
  geom_point(data = DFT1.SNVs.counts.SBS5[match(samples.cladeC[,8], rownames(DFT1.SNVs.counts.SBS5)),], 
             color = alpha("cornflowerblue", 0.6), 
             size = 9, shape = 'triangle') + 
  geom_smooth(data = DFT1.SNVs.counts.SBS5[match(samples.cladeC[,8], rownames(DFT1.SNVs.counts.SBS5)),],
              method = 'lm', color = "cornflowerblue", fullrange = T, fill = alpha('cornflowerblue', 0.3)) +
  labs(y = "SBS5 Substitutions", x = "Year") +
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

## clean up environment
rm(list=ls())
