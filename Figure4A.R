## Figure 4A
## 'The evolution of two transmissible cancers in Tasmanian devils'
## Stammnitz et al., 2022
## maxrupsta@gmail

library(readxl)
library(ggplot2)
library(lubridate)
library(stringr)

## set input path(s)
setwd('/Tables')


# LINE-1 rates #
################

## import LINE-1 lists
DFT.LINE1.rates <- as.matrix(read_xlsx('Table-S4.xlsx', sheet = 1))
colnames(DFT.LINE1.rates) <- DFT.LINE1.rates[2,]
DFT.LINE1.rates <- DFT.LINE1.rates[-c(1:2),]
DFT1.LINE1.counts <- DFT.LINE1.rates[DFT.LINE1.rates[,1] == 'DFT1',]
DFT1.LINE1.counts[1,] <- gsub('\\*', '', DFT1.LINE1.counts[1,])
DFT2.LINE1.counts <- DFT.LINE1.rates[DFT.LINE1.rates[,1] == 'DFT2',]
DFT2.LINE1.counts[34,] <- gsub('\\*', '', DFT2.LINE1.counts[34,])

## load metadata
samples <- as.matrix(read_xlsx('Table-S2.xlsx', sheet = 1))

### DFT1
dft1.samples <- samples[samples[,9] == 'DFT1',]
dft1.samples <- dft1.samples[which(is.na(dft1.samples[,9]) == F),]
dft1.samples <- dft1.samples[,c(8,4,12)]
dft1.samples[,1] <- str_split_fixed(dft1.samples[,1], ' ', 2)[,1]
dft1.samples[,2] <- gsub('[*]', '', dft1.samples[,2])
colnames(dft1.samples) <- c('Tumour', 'Collection Date', 'Purity')

### DFT2
dft2.samples <- samples[samples[,9] == 'DFT2',]
dft2.samples <- dft2.samples[which(is.na(dft2.samples[,9]) == F),]
dft2.samples <- dft2.samples[,c(8,4,12)]
dft2.samples[,1] <- str_split_fixed(dft2.samples[,1], ' ', 2)[,1]
dft2.samples[,2] <- gsub('[*]', '', dft2.samples[,2])
colnames(dft2.samples) <- c('Tumour', 'Collection Date', 'Purity')

## combine
DFT1.DFT2.L1.counts <- rbind(dft1.samples, dft2.samples)
DFT1.DFT2.L1.counts <- cbind(DFT1.DFT2.L1.counts[c(match(DFT1.LINE1.counts[,"TUMOUR ID"], DFT1.DFT2.L1.counts[,1]),
                                                   match(DFT2.LINE1.counts[,"TUMOUR ID"], DFT1.DFT2.L1.counts[,1])),],
                             c(DFT1.LINE1.counts[,"LINE-1 INSERTIONS"], DFT2.LINE1.counts[,"LINE-1 INSERTIONS"]))
colnames(DFT1.DFT2.L1.counts)[4] <- 'L1s'
DFT1.DFT2.L1.counts <- DFT1.DFT2.L1.counts[,-1]
DFT1.DFT2.L1.counts[,'Collection Date'] <- decimal_date(as.Date(as.character(unlist(DFT1.DFT2.L1.counts[,'Collection Date'])), format="%d.%m.%Y"))
DFT1.DFT2.L1.counts <- DFT1.DFT2.L1.counts[!is.na(DFT1.DFT2.L1.counts[,'Collection Date']),]
DFT1.DFT2.L1.counts <- as.data.frame(DFT1.DFT2.L1.counts)
DFT1.DFT2.L1.counts[,'Collection Date'] <- as.numeric(as.character(DFT1.DFT2.L1.counts[,'Collection Date']))
DFT1.DFT2.L1.counts[,'Purity'] <- as.numeric(as.character(DFT1.DFT2.L1.counts[,'Purity']))
DFT1.DFT2.L1.counts[,'L1s'] <- as.numeric(as.character(DFT1.DFT2.L1.counts[,'L1s']))

## plot
pdf("Figure4A_L1_rates.pdf", 
    height = 12, width = 10)
ggplot(DFT1.DFT2.L1.counts, aes(x = `Collection Date`, y = `L1s`)) + 
  scale_x_continuous(breaks = seq(from = 1990, to = 2030, by = 20), limits = c(1985, 2035)) + 
  scale_y_continuous(breaks = seq(from = 0, to = 200, by = 40), limits = c(-500, 500)) + 
  coord_cartesian(xlim = c(1985, 2035), ylim = c(0, 200), expand = F) +
  geom_point(color = c(rep("cornflowerblue", 64), rep("red", 39)), size = 4) + 
  geom_smooth(data = DFT1.DFT2.L1.counts[1:64,],
              method = 'lm', color = "cornflowerblue", fullrange = T) +
  geom_smooth(data = DFT1.DFT2.L1.counts[65:nrow(DFT1.DFT2.L1.counts),], 
              method = 'lm', color = "red", fullrange = T) +
  theme_classic(base_size = 20) +
  labs(y = "LINE-1 insertions", x = "Year") +
  theme(axis.text = element_text(size = 35),
        axis.line.x = element_line(size = 1, color = 'black'),
        axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = 50, vjust = -1),
        axis.title.y = element_text(family = 'Helvetica', colour = 'black', size = 50, vjust = 3),
        legend.position = "none",
        text = element_text(family="Helvetica"),
        plot.margin = unit(c(2, 2, 2, 2),"cm"))
dev.off()

## clean up environment
rm(list=ls())
