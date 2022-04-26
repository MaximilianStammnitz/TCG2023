## Figure 5B
## 'The evolution of two transmissible cancers in Tasmanian devils'
## Stammnitz et al., 2022
## maxrupsta@gmail

library(readxl)
library(stringr)
library(lubridate)
library(ggplot2)

## set input path(s)
setwd('/Users/ms37/Desktop/Labwork/DFT_evolution/doc/manuscripts/The Evolutionary History of Two Transmissible Cancers in Tasmanian Devils/Tables/v6/')


# Rearrangement event rates #
#############################

## import rearrangement events
DFT.SVs <- as.matrix(read_xlsx('/Users/mstammnitz/Desktop/DFT_evolution/doc/manuscripts/The Evolutionary History of Two Transmissible Cancers in Tasmanian Devils/Tables/v6/Table-S5_v6.xlsx', sheet = 1))
colnames(DFT.SVs) <- as.character(DFT.SVs[2,])
DFT.SVs <- DFT.SVs[-c(1:2),]
DFT1.SVs <- DFT.SVs[DFT.SVs[,'LINEAGE'] == 'DFT1',]
DFT2.SVs <- DFT.SVs[DFT.SVs[,'LINEAGE'] == 'DFT2',]

## load metadata
samples <- as.matrix(read_xlsx('/Users/mstammnitz/Desktop/DFT_evolution/doc/manuscripts/The Evolutionary History of Two Transmissible Cancers in Tasmanian Devils/Tables/v6/Table-S2_v6.xlsx', sheet = 1))

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
DFT1.DFT2.SV.event.counts <- rbind(dft1.samples, dft2.samples)
DFT1.DFT2.SV.event.counts <- cbind(DFT1.DFT2.SV.event.counts[c(match(DFT1.SVs[,"TUMOUR ID"], DFT1.DFT2.SV.event.counts[,1]),
                                                               match(DFT2.SVs[,"TUMOUR ID"], DFT1.DFT2.SV.event.counts[,1])),],
                                   c(DFT1.SVs[,"REARRANGEMENT EVENTS"], DFT2.SVs[,"REARRANGEMENT EVENTS"]))
colnames(DFT1.DFT2.SV.event.counts)[4] <- 'SV events'
DFT1.DFT2.SV.event.counts <- DFT1.DFT2.SV.event.counts[,-1]
DFT1.DFT2.SV.event.counts[,'Collection Date'] <- decimal_date(as.Date(as.character(unlist(DFT1.DFT2.SV.event.counts[,'Collection Date'])), format="%d.%m.%Y"))
DFT1.DFT2.SV.event.counts <- DFT1.DFT2.SV.event.counts[!is.na(DFT1.DFT2.SV.event.counts[,'Collection Date']),]
DFT1.DFT2.SV.event.counts <- as.data.frame(DFT1.DFT2.SV.event.counts)
DFT1.DFT2.SV.event.counts[,'Collection Date'] <- as.numeric(as.character(DFT1.DFT2.SV.event.counts[,'Collection Date']))
DFT1.DFT2.SV.event.counts[,'Purity'] <- as.numeric(as.character(DFT1.DFT2.SV.event.counts[,'Purity']))
DFT1.DFT2.SV.event.counts[,'SV events'] <- as.numeric(as.character(DFT1.DFT2.SV.event.counts[,'SV events']))

## plot
pdf("Figure5B_Rearrangement_event_rates.pdf", 
    height = 12, width = 10)

ggplot(DFT1.DFT2.SV.event.counts, aes(x = `Collection Date`, y = `SV events`)) + 
  scale_x_continuous(breaks = seq(from = 1990, to = 2030, by = 20), limits = c(1985, 2035)) + 
  scale_y_continuous(breaks = seq(from = 0, to = 100, by = 20), limits = c(-200, 300)) + 
  coord_cartesian(xlim = c(1985, 2035), ylim = c(0, 100), expand = F) +
  geom_point(color = c(rep("cornflowerblue", 64), rep("red", 39)), size = 5) + 
  geom_smooth(data = DFT1.DFT2.SV.event.counts[1:64,],
              method = 'lm', color = "cornflowerblue", fullrange = T) +
  geom_smooth(data = DFT1.DFT2.SV.event.counts[65:nrow(DFT1.DFT2.SV.event.counts),], 
              method = 'lm', color = "red", fullrange = T) +
  labs(y = "Rearrangement events", x = "Year") +
  theme_classic(base_size = 20) +
  theme(axis.text = element_text(size = 45),
        axis.line.x = element_line(size = 1, color = 'black'),
        axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = 60, vjust = -1),
        axis.title.y = element_text(family = 'Helvetica', colour = 'black', size = 60, vjust = 3),
        legend.position = "none",
        text = element_text(family="Helvetica"),
        plot.margin = unit(c(2, 2, 2, 2),"cm"))
dev.off()

## clean up environment
rm(list=ls())


# CNV rates #
#############

## DFT1
DFT1.cnvs <- as.matrix(read_xlsx('/Users/mstammnitz/Desktop/DFT_evolution/doc/manuscripts/The Evolutionary History of Two Transmissible Cancers in Tasmanian Devils/Tables/v6/Table-S6_v6.xlsx', sheet = 2))
colnames(DFT1.cnvs) <- DFT1.cnvs[2,]
DFT1.cnvs <- DFT1.cnvs[-c(1:2),]
DFT1.cnvs <- DFT1.cnvs[-which(DFT1.cnvs[,'MARKER 5'] == 'YES'),] ## M5s
DFT1.cnvs <- DFT1.cnvs[-which(DFT1.cnvs[,'DOUBLE MINUTE'] == 'YES'),] ## DMs

### half the count of post-tetraploidisation events
dft1.tetraploids <- c("18T", "49T1", "209T3", "372T1", "377T3", 
                      "378T1", "398T1", "420T1", "421T1", "477T3", 
                      "528T2", "876T1", "1191T1", "2690T", "2694Ta", "2694Tb")
for (i in 1:nrow(DFT1.cnvs)){

  print(i)
  
  ## fetch data
  if(as.character(DFT1.cnvs[i,'EVENT TYPE']) == 'GAIN' | 
     as.character(DFT1.cnvs[i,'EVENT TYPE'])  == 'GAIN (pre-tetraploidisation)' |
     as.character(DFT1.cnvs[i,'EVENT TYPE']) == 'LOSS' | 
     as.character(DFT1.cnvs[i,'EVENT TYPE'])  == 'LOSS (pre-tetraploidisation)'){
    
    next
    
  }else if (as.character(DFT1.cnvs[i,'EVENT TYPE']) == 'GAIN (post-tetraploidisation)' | 
            as.character(DFT1.cnvs[i,'EVENT TYPE']) == 'LOSS (post-tetraploidisation)'){
    
    for(j in match(dft1.tetraploids, colnames(DFT1.cnvs))){
      
      if(DFT1.cnvs[i,j] == 0){
        
        next
        
      }else if (DFT1.cnvs[i,j] == 1){
        
        DFT1.cnvs[i,j] <- 0.5
        
      }

    }

  }

}

## DFT2
DFT2.cnvs <- as.matrix(read_xlsx('/Users/mstammnitz/Desktop/DFT_evolution/doc/manuscripts/The Evolutionary History of Two Transmissible Cancers in Tasmanian Devils/Tables/v6/Table-S6_v6.xlsx', sheet = 3))
colnames(DFT2.cnvs) <- DFT2.cnvs[2,]
DFT2.cnvs <- DFT2.cnvs[-c(1:2),]

### half the count of post-tetraploidisation events
dft2.tetraploids <- c("1531T1", "1538T1", "1553T1")
for (i in 1:nrow(DFT2.cnvs)){
  
  print(i)
  
  ## fetch data
  if(as.character(DFT2.cnvs[i,'EVENT TYPE']) == 'GAIN' | 
     as.character(DFT2.cnvs[i,'EVENT TYPE'])  == 'GAIN (pre-tetraploidisation)' |
     as.character(DFT2.cnvs[i,'EVENT TYPE']) == 'LOSS' | 
     as.character(DFT2.cnvs[i,'EVENT TYPE'])  == 'LOSS (pre-tetraploidisation)'){
    
    next
    
  }else if (as.character(DFT2.cnvs[i,'EVENT TYPE']) == 'GAIN (post-tetraploidisation)' | 
            as.character(DFT2.cnvs[i,'EVENT TYPE']) == 'LOSS (post-tetraploidisation)'){
    
    for(j in match(dft2.tetraploids, colnames(DFT2.cnvs))){
      
      if(DFT2.cnvs[i,j] == 0){
        
        next
        
      }else if (DFT2.cnvs[i,j] == 1){
        
        DFT2.cnvs[i,j] <- 0.5
        
      }
      
    }
    
  }
  
}

## derive rates
total.cnvs.dft1 <- apply(DFT1.cnvs[,-c(1:12)], 2, function(x){sum(abs(as.numeric(x)))})
total.cnvs.dft2 <- apply(DFT2.cnvs[,-c(1:10)], 2, function(x){sum(abs(as.numeric(x)))})

## load metadata
samples <- as.matrix(read_xlsx('/Users/mstammnitz/Desktop/DFT_evolution/doc/manuscripts/The Evolutionary History of Two Transmissible Cancers in Tasmanian Devils/Tables/v6/Table-S2_v6.xlsx', sheet = 1))

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
DFT1.DFT2.CNVs.counts <- rbind(dft1.samples, dft2.samples)
DFT1.DFT2.CNVs.counts <- cbind(DFT1.DFT2.CNVs.counts[c(match(names(total.cnvs.dft1), DFT1.DFT2.CNVs.counts[,1]),
                                                       match(names(total.cnvs.dft2), DFT1.DFT2.CNVs.counts[,1])),],
                               c(total.cnvs.dft1, total.cnvs.dft2))
colnames(DFT1.DFT2.CNVs.counts)[4] <- 'CNVs'
DFT1.DFT2.CNVs.counts <- DFT1.DFT2.CNVs.counts[,-1]
DFT1.DFT2.CNVs.counts[,'Collection Date'] <- decimal_date(as.Date(as.character(unlist(DFT1.DFT2.CNVs.counts[,'Collection Date'])), format="%d.%m.%Y"))
DFT1.DFT2.CNVs.counts <- DFT1.DFT2.CNVs.counts[!is.na(DFT1.DFT2.CNVs.counts[,'Collection Date']),]
DFT1.DFT2.CNVs.counts <- as.data.frame(DFT1.DFT2.CNVs.counts)
DFT1.DFT2.CNVs.counts[,'Collection Date'] <- as.numeric(as.character(DFT1.DFT2.CNVs.counts[,'Collection Date']))
DFT1.DFT2.CNVs.counts[,'Purity'] <- as.numeric(as.character(DFT1.DFT2.CNVs.counts[,'Purity']))
DFT1.DFT2.CNVs.counts[,'CNVs'] <- as.numeric(as.character(DFT1.DFT2.CNVs.counts[,'CNVs']))

## plot
pdf("Figure5B_CNV_rates.pdf", 
    height = 12, width = 10)

ggplot(DFT1.DFT2.CNVs.counts, aes(x = `Collection Date`, y = `CNVs`)) + 
  scale_x_continuous(breaks = seq(from = 1990, to = 2030, by = 20), limits = c(1985, 2035)) + 
  scale_y_continuous(breaks = seq(from = 0, to = 100, by = 20), limits = c(-200, 300)) + 
  coord_cartesian(xlim = c(1985, 2035), ylim = c(0, 100), expand = F) +
  geom_point(color = c(rep("cornflowerblue", 76), rep("red", 39)), size = 5) + 
  geom_smooth(data = DFT1.DFT2.CNVs.counts[1:76,],
              method = 'lm', color = "cornflowerblue", fullrange = T) +
  geom_smooth(data = DFT1.DFT2.CNVs.counts[77:nrow(DFT1.DFT2.CNVs.counts),], 
              method = 'lm', color = "red", fullrange = T) +
  labs(y = "Copy number variants", x = "Year") +
  theme_classic(base_size = 20) +
  theme(axis.text = element_text(size = 45),
        axis.line.x = element_line(size = 1, color = 'black'),
        axis.title.x = element_text(family = 'Helvetica', colour = 'black', size = 60, vjust = -1),
        axis.title.y = element_text(family = 'Helvetica', colour = 'black', size = 60, vjust = 3),
        legend.position = "none",
        text = element_text(family="Helvetica"),
        plot.margin = unit(c(2, 2, 2, 2),"cm"))
dev.off()

## clean up environment
rm(list=ls())