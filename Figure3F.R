## Figure 3F
## 'The evolution of two transmissible cancers in Tasmanian devils'
## Stammnitz et al., 2022
## maxrupsta@gmail

library(readxl)
library(stringr)
library(lubridate)
library(ggplot2)

## set input path(s)
setwd('/Tables')


# DFT1 and DFT2 ID1 rates #
############################

## import data
counts <- as.matrix(read_xlsx('Table-S3.xlsx', sheet = 4))
colnames(counts) <- as.character(counts[2,])
counts <- counts[-c(1:2),]
dft1.counts <- counts[counts[,'LINEAGE'] == 'DFT1',]
dft1.counts <- dft1.counts[which(dft1.counts[,'TUMOUR ID'] != '377T1'),]
dft2.counts <- counts[counts[,'LINEAGE'] == 'DFT2',]

## metadata
samples <- as.matrix(read_xlsx('Table-S2.xlsx', sheet = 1))

### DFT1
dft1.samples <- samples[samples[,9] == 'DFT1',]
dft1.samples <- dft1.samples[which(is.na(dft1.samples[,9]) == F),]
dft1.samples <- dft1.samples[,c(8,4,13)]
dft1.samples[,1] <- str_split_fixed(dft1.samples[,1], ' ', 2)[,1]
dft1.samples[,2] <- gsub('[*]', '', dft1.samples[,2])
dft1.samples <- dft1.samples[which(dft1.samples[,1] != '377T1'),]
colnames(dft1.samples) <- c('Tumour', 'Collection Date', 'Purity')

### DFT2
dft2.samples <- samples[samples[,9] == 'DFT2',]
dft2.samples <- dft2.samples[which(is.na(dft2.samples[,9]) == F),]
dft2.samples <- dft2.samples[,c(8,4,13)]
dft2.samples[,1] <- str_split_fixed(dft2.samples[,1], ' ', 2)[,1]
dft2.samples[,2] <- gsub('[*]', '', dft2.samples[,2])
colnames(dft2.samples) <- c('Tumour', 'Collection Date', 'Purity')

## combine
DFT1.DFT2.ID1.counts <- rbind(dft1.samples, dft2.samples)
DFT1.DFT2.ID1.counts <- cbind(DFT1.DFT2.ID1.counts[c(match(dft1.counts[,'TUMOUR ID'], DFT1.DFT2.ID1.counts[,1]),
                                                       match(dft2.counts[,'TUMOUR ID'], DFT1.DFT2.ID1.counts[,1])),],
                               c(dft1.counts[,"WGD-NORMALISED ID1 INDELS"], dft2.counts[,"WGD-NORMALISED ID1 INDELS"]))
colnames(DFT1.DFT2.ID1.counts)[4] <- 'ID1'
rownames(DFT1.DFT2.ID1.counts) <- DFT1.DFT2.ID1.counts[,1]
DFT1.DFT2.ID1.counts <- DFT1.DFT2.ID1.counts[,-1]
DFT1.DFT2.ID1.counts[,'Collection Date'] <- decimal_date(as.Date(as.character(unlist(DFT1.DFT2.ID1.counts[,'Collection Date'])), format="%d.%m.%Y"))
DFT1.DFT2.ID1.counts <- DFT1.DFT2.ID1.counts[!is.na(DFT1.DFT2.ID1.counts[,'Collection Date']),]
DFT1.DFT2.ID1.counts <- as.data.frame(DFT1.DFT2.ID1.counts)
DFT1.DFT2.ID1.counts[,'Collection Date'] <- as.numeric(as.character(DFT1.DFT2.ID1.counts[,'Collection Date']))
DFT1.DFT2.ID1.counts[,'Purity'] <- as.numeric(as.character(DFT1.DFT2.ID1.counts[,'Purity']))
DFT1.DFT2.ID1.counts[,'ID1'] <- as.numeric(as.character(DFT1.DFT2.ID1.counts[,'ID1']))

pdf("Figure3F_ID1_rates.pdf", 
    height = 12, width = 18)
ggplot(DFT1.DFT2.ID1.counts, aes(x = `Collection Date`, y = `ID1`)) + 
  scale_x_continuous(breaks = seq(from = 1990, to = 2030, by = 10), limits = c(1985, 2035)) + 
  scale_y_continuous(breaks = seq(from = 0, to = 1000, by = 200), limits = c(-1000, 1000)) + 
  coord_cartesian(xlim = c(1985, 2035), ylim = c(0, 1000), expand = F) +
  geom_point(color = c(rep("cornflowerblue", 75), rep("red", 39)), size = 5) + 
  geom_smooth(data = DFT1.DFT2.ID1.counts[1:75,],
              method = 'lm', color = "cornflowerblue", fullrange = T) +
  geom_smooth(data = DFT1.DFT2.ID1.counts[76:nrow(DFT1.DFT2.ID1.counts),], 
              method = 'lm', color = "red", fullrange = T) +
  labs(y = "ID1 Indels", x = "Year") +
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


# DFT1 and DFT2 ID2 rates #
############################

## import data
counts <- as.matrix(read_xlsx('Table-S3.xlsx', sheet = 4))
colnames(counts) <- as.character(counts[2,])
counts <- counts[-c(1:2),]
dft1.counts <- counts[counts[,'LINEAGE'] == 'DFT1',]
dft1.counts <- dft1.counts[which(dft1.counts[,'TUMOUR ID'] != '377T1'),]
dft2.counts <- counts[counts[,'LINEAGE'] == 'DFT2',]

## metadata
samples <- as.matrix(read_xlsx('Table-S2.xlsx', sheet = 1))

### DFT1
dft1.samples <- samples[samples[,9] == 'DFT1',]
dft1.samples <- dft1.samples[which(is.na(dft1.samples[,9]) == F),]
dft1.samples <- dft1.samples[,c(8,4,13)]
dft1.samples[,1] <- str_split_fixed(dft1.samples[,1], ' ', 2)[,1]
dft1.samples[,2] <- gsub('[*]', '', dft1.samples[,2])
dft1.samples <- dft1.samples[which(dft1.samples[,1] != '377T1'),]
colnames(dft1.samples) <- c('Tumour', 'Collection Date', 'Purity')

### DFT2
dft2.samples <- samples[samples[,9] == 'DFT2',]
dft2.samples <- dft2.samples[which(is.na(dft2.samples[,9]) == F),]
dft2.samples <- dft2.samples[,c(8,4,13)]
dft2.samples[,1] <- str_split_fixed(dft2.samples[,1], ' ', 2)[,1]
dft2.samples[,2] <- gsub('[*]', '', dft2.samples[,2])
colnames(dft2.samples) <- c('Tumour', 'Collection Date', 'Purity')

## combine
DFT1.DFT2.ID2.counts <- rbind(dft1.samples, dft2.samples)
DFT1.DFT2.ID2.counts <- cbind(DFT1.DFT2.ID2.counts[c(match(dft1.counts[,'TUMOUR ID'], DFT1.DFT2.ID2.counts[,1]),
                                                       match(dft2.counts[,'TUMOUR ID'], DFT1.DFT2.ID2.counts[,1])),],
                               c(dft1.counts[,"WGD-NORMALISED ID2 INDELS"], dft2.counts[,"WGD-NORMALISED ID2 INDELS"]))
colnames(DFT1.DFT2.ID2.counts)[4] <- 'ID2'
rownames(DFT1.DFT2.ID2.counts) <- DFT1.DFT2.ID2.counts[,1]
DFT1.DFT2.ID2.counts <- DFT1.DFT2.ID2.counts[,-1]
DFT1.DFT2.ID2.counts[,'Collection Date'] <- decimal_date(as.Date(as.character(unlist(DFT1.DFT2.ID2.counts[,'Collection Date'])), format="%d.%m.%Y"))
DFT1.DFT2.ID2.counts <- DFT1.DFT2.ID2.counts[!is.na(DFT1.DFT2.ID2.counts[,'Collection Date']),]
DFT1.DFT2.ID2.counts <- as.data.frame(DFT1.DFT2.ID2.counts)
DFT1.DFT2.ID2.counts[,'Collection Date'] <- as.numeric(as.character(DFT1.DFT2.ID2.counts[,'Collection Date']))
DFT1.DFT2.ID2.counts[,'Purity'] <- as.numeric(as.character(DFT1.DFT2.ID2.counts[,'Purity']))
DFT1.DFT2.ID2.counts[,'ID2'] <- as.numeric(as.character(DFT1.DFT2.ID2.counts[,'ID2']))

pdf("Figure3F_ID2_rates.pdf", 
    height = 12, width = 18)
ggplot(DFT1.DFT2.ID2.counts, aes(x = `Collection Date`, y = `ID2`)) + 
  scale_x_continuous(breaks = seq(from = 1990, to = 2030, by = 10), limits = c(1985, 2035)) + 
  scale_y_continuous(breaks = seq(from = 0, to = 500, by = 100), limits = c(-750, 1000)) + 
  coord_cartesian(xlim = c(1985, 2035), ylim = c(0, 500), expand = F) +
  geom_point(color = c(rep("cornflowerblue", 75), rep("red", 39)), size = 5) + 
  geom_smooth(data = DFT1.DFT2.ID2.counts[1:75,],
              method = 'lm', color = "cornflowerblue", fullrange = T) +
  geom_smooth(data = DFT1.DFT2.ID2.counts[76:nrow(DFT1.DFT2.ID2.counts),], 
              method = 'lm', color = "red", fullrange = T) +
  labs(y = "ID2 Indels", x = "Year") +
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
