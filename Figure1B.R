## Figure 1B
## 'The evolution of two transmissible cancers in Tasmanian devils'
## Stammnitz et al., 2022
## maxrupsta@gmail

## libraries
library(readxl)
library(stringr)
library(ggplot2)
library(ggmap)
library(rstudioapi)

## set input path(s)
setwd('/Users/ms37/Desktop/Labwork/DFT_evolution/doc/manuscripts/The Evolutionary History of Two Transmissible Cancers in Tasmanian Devils/Tables/v6/')


# DFT1 sampling map #
#####################

## fetch location of all deep-sequenced DFT1 tumours
samples <- as.matrix(read_xlsx('Table-S2_v6.xlsx', sheet = 1))
dft1.samples <- samples[samples[,9] == 'DFT1',]
dft1.samples <- dft1.samples[which(is.na(dft1.samples[,9]) == F),]
dft1.samples <- dft1.samples[,c(8,5,6,7)]
dft1.samples[,1] <- str_split_fixed(dft1.samples[,1], ' ', 2)[,1]
dft1.samples[,2] <- str_split_fixed(dft1.samples[,2], '\\[', 2)[,1]
colnames(dft1.samples) <- c('Tumour', 'Native Location', 'Latitude', 'Longitude')

## add clade info
clade_assignments <- samples[,c(8,10)]
dft1.samples <- cbind(dft1.samples, 'Clade' = clade_assignments[match(dft1.samples[,'Tumour'],clade_assignments[,1]),2])
dft1.samples[,'Clade'] <- str_split_fixed(dft1.samples[,'Clade'], '-', 2)[,2]

## add a bit of (uniform) spatial noise to locations sampled more than once  
set.seed(100)
add.noise <- names(table(dft1.samples[,'Native Location'])[table(dft1.samples[,'Native Location']) > 1])
mult <- dft1.samples[dft1.samples[,'Native Location'] %in% add.noise,]
sing <- dft1.samples[!dft1.samples[,'Native Location'] %in% add.noise,]
for (i in 1:length(add.noise)){
  print(i)
  add.noise.tmp <- add.noise[i]
  add.noise.tmp.id <- which(mult[,'Native Location'] == add.noise.tmp)
  mult[add.noise.tmp.id,'Latitude'] <- as.numeric(mult[add.noise.tmp.id,'Latitude']) + runif(n = length(add.noise.tmp.id), min = -0.1, max = 0.1)
  mult[add.noise.tmp.id,'Longitude'] <- as.numeric(mult[add.noise.tmp.id,'Longitude']) + runif(n = length(add.noise.tmp.id), min = -0.1, max = 0.1)
}

## add a single DFT2 data point
summary <- rbind(mult,sing, c('DFT2', 'Channel', '-43.155303', '147.172596', NA))
rownames(summary) <- summary[,'Tumour']
summary <- summary[,3:4]
class(summary) <- 'numeric'
summary <- as.data.frame(summary)

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
pdf('Figure1B_map_DFT1.pdf', width = 7, height = 6)
ggmap(tasmania.map.googlemaps_terrain, extent = 'panel', padding = 0)  +
  scale_x_continuous(limits = c(144.3, 148.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-43.7, -40.55), expand = c(0, 0)) +
  geom_point(data = summary,
             mapping = aes(x = Longitude, y = Latitude, group = NULL), 
             size = c(rep(8, 78), 12),
             show.legend = F,
             colour = c(rep('cornflowerblue', 78), 'red'),
             alpha = 0.5) +
  theme(panel.background = element_rect(fill = 'transparent'),
        plot.background = element_rect(fill = 'transparent', color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill = 'transparent'),
        legend.box.background = element_rect(fill = 'transparent'),
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  scalebar(x.min = 144.38, 
           x.max = 148.5, 
           y.min = -43.65, 
           y.max = -40.55,
           group = NA,
           location = 'bottomleft', 
           dist = 50,
           dist_unit = 'km', 
           transform = TRUE,  
           model = 'WGS84',
           st.color = 'black',
           st.dist = 0.04,
           st.size = 8.5,
           st.bottom = F,
           box.fill = c('grey70', 'grey90'),
           box.color = 'black',
           height = 0.02)
dev.off()

## clean up environment
rm(list=ls())


# DFT2 sampling map #
#####################

## fetch location of all deep-sequenced DFT2 tumours
samples <- as.matrix(read_xlsx('Table-S2_v6.xlsx', sheet = 1))
dft2.samples <- samples[samples[,9] == 'DFT2',]
dft2.samples <- dft2.samples[which(is.na(dft2.samples[,9]) == F),]
dft2.samples <- dft2.samples[,c(8,5,6,7)]
dft2.samples <- dft2.samples[-which(dft2.samples[,2] == 'N/A'),] ## no spatial information for 1525T1
dft2.samples[,1] <- str_split_fixed(dft2.samples[,1], ' ', 2)[,1]
colnames(dft2.samples) <- c('Tumour', 'Native Location', 'Latitude', 'Longitude')

## add clade info
clade_assignments <- samples[,c(8,10)]
dft2.samples <- cbind(dft2.samples, 'Clade' = clade_assignments[match(dft2.samples[,'Tumour'],clade_assignments[,1]),2])
dft2.samples[,'Clade'] <- str_split_fixed(dft2.samples[,'Clade'], '-', 2)[,2]

## add a bit of (uniform) spatial noise to locations sampled more than once  
add.noise <- names(table(dft2.samples[,'Native Location'])[table(dft2.samples[,'Native Location']) > 1])
mult <- dft2.samples[dft2.samples[,'Native Location'] %in% add.noise,]
sing <- dft2.samples[!dft2.samples[,'Native Location'] %in% add.noise,]
for (i in 1:length(add.noise)){
  print(i)
  add.noise.tmp <- add.noise[i]
  add.noise.tmp.id <- which(mult[,'Native Location'] == add.noise.tmp)
  mult[add.noise.tmp.id,'Latitude'] <- as.numeric(mult[add.noise.tmp.id,'Latitude']) + runif(n = length(add.noise.tmp.id), min = -0.003, max = 0.003)
  mult[add.noise.tmp.id,'Longitude'] <- as.numeric(mult[add.noise.tmp.id,'Longitude']) + runif(n = length(add.noise.tmp.id), min = -0.003, max = 0.003)
}

## add two DFT1 data points
dft2.samples <- rbind(mult, sing, 
                      c('812T2', 'Lower Snug, Powers Road', '-43.079632933333301','147.24856074459299', NA),
                      c('818T2', 'Lower Snug, Powers Road', '-43.079632933333301','147.24856074459299', NA))
rownames(dft2.samples) <- dft2.samples[,'Tumour']
dft2.samples <- dft2.samples[,3:4]
class(dft2.samples) <- 'numeric'
dft2.samples <- as.data.frame(dft2.samples)

## downlaod Channel map via ggmap package, requires a Google API key:
### video link ###
## register_google(key = '') fill in
channel.map.googlemaps_terrain <- get_googlemap(center = c(147.172596, -43.155303), 
                                                zoom = 11,
                                                size = c(640, 640),
                                                scale = 2,
                                                maptype = 'terrain',
                                                color = 'bw',
                                                style = 'style=feature:administrative.country|element:labels|visibility:off&style=feature:water|color:0xffffff|visibility:simplified&style=feature:road|color:0xffffff|visibility:off')


## plot
pdf('Figure1B_map_DFT2.pdf', width = 7, height = 6)
ggmap(channel.map.googlemaps_terrain,
      extent = 'panel',
      padding = 0)  +
  scale_x_continuous(limits = c(146.87, 147.34), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-43.3, -43.00), expand = c(0, 0)) +
  geom_point(data = dft2.samples,
             mapping = aes(x = Longitude, y = Latitude, group = NULL), 
             size = rep(8, nrow(dft2.samples)),
             show.legend = F,
             colour = c(rep('red',40),
                        rep('cornflowerblue', 2)),
             alpha = 0.5) +
  theme(panel.background = element_rect(fill = 'transparent'),
        plot.background = element_rect(fill = 'transparent', color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill = 'transparent'),
        legend.box.background = element_rect(fill = 'transparent'),
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  scalebar(x.min = 146.965, 
           x.max = 147.34, 
           y.min = -43.295, 
           y.max = -43.00,
           group = NA,
           location = 'bottomleft', 
           dist = 5,
           dist_unit = 'km', 
           transform = TRUE,  
           model = 'WGS84',
           st.color = 'black',
           st.bottom = F,
           st.dist = 0.04,
           st.size = 8.5,
           box.fill = c('grey70', 'grey90'),
           box.color = 'black',
           height = 0.02)
dev.off()

## clean up environment
rm(list=ls())