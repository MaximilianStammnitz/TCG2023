## Figure 6D
## 'The evolution of two transmissible cancers in Tasmanian devils'
## Stammnitz et al., 2022
## maxrupsta@gmail

## set input path(s)
setwd('/Tables')


# MGA truncation variant display #
##################################

## load MGA annotation
load("Sarcophilus_harrisii.mSarHar1.11.102.gtf.Rdata")
MGA.ensembl.out <- ensembl[grep('MGA', ensembl[,9]),]
MGA.ensembl.out <- MGA.ensembl.out[which(MGA.ensembl.out[,1] == '2' & MGA.ensembl.out[,4] < 263000000),]
MGA.ensembl.out.exons <- MGA.ensembl.out[which(MGA.ensembl.out[,3] == 'exon'),]
MGA.ensembl.out.exons <- MGA.ensembl.out.exons[grep('MGA-201', MGA.ensembl.out.exons[,9]),]

pdf("Figure6D_MGA_mutations.pdf",
    width = 16, height = 5)
par(mar = c(8,1,1,1))

cexs <- 2

plot(x = c(262870000, 262980000),
     y = rep(1,2.3), col = 'white',
     yaxt = 'n', frame = F, ylim = c(0,2.3), ylab = '', 
     xaxt = 'n', xlab = '', xlim = c(262870000, 262980000),
     main = '', cex.main = 2)

## add exon polygons and intron zick-zacks
for(i in 4:nrow(MGA.ensembl.out.exons)){
  
  ## exon "squares"
  polygon(x = as.numeric(MGA.ensembl.out.exons[i,c(4,5,5,4)]),
          y = c(0.5,0.5,0.75,0.75), 
          col = 'black',
          border = 'black')
  
  ## intron zick-zack
  if(i != nrow(MGA.ensembl.out.exons)){
    
    tmp.intron.middle <- mean(c(as.numeric(MGA.ensembl.out.exons[i,5]),
                                as.numeric(MGA.ensembl.out.exons[c(i+1),4])))
    
    lines(x = c(as.numeric(MGA.ensembl.out.exons[i,5]),
                tmp.intron.middle),
          y = c(0.625,0.75),
          col = 'black')
    
    lines(x = c(tmp.intron.middle, as.numeric(MGA.ensembl.out.exons[c(i+1),4])),
          y = c(0.75,0.625),
          col = 'black') 
    
  }
  
}

## display 5' UTR
text(x = 262877000,
     y = 0.625, 
     labels = "5' UTR", cex = cexs)
lines(x = c(262883000,262885457),
      y = c(0.625, 0.625),
      col = 'black',
      lty = 2)

## display natural ATG
text(x = 262885535,
     y = 0.15, 
     labels = "ATG", 
     pos = 2,
     cex = cexs,
     offset = -0.1)
lines(x = c(262885725,262885725),
      y = c(0.3,0.5),
      col = 'black')

## display natural TGA
text(x = 262974035,
     y = 0.15, 
     labels = "TGA", 
     pos = 4,
     cex = cexs,
     offset = -0.1)
lines(x = c(262974035,262974035),
      y = c(0.3,0.5),
      col = 'black')

## display mutation #1
lines(x = c(262885725,262885725),
      y = c(0.75,1.25),
      col = 'cornflowerblue',
      lwd = 2)
points(x = 262885725,
       y = 1.25,
       col = 'cornflowerblue',
       bg = 'cornflowerblue',
       cex = cexs*3, 
       pch = 25)
text(x = 262885725,
     y = 1.55, 
     labels = "VK64-65VX", 
     col = 'cornflowerblue',
     cex = cexs)

## display mutation #2
lines(x = c(262901923,262901923),
      y = c(0.75,1.25),
      col = 'cornflowerblue',
      lwd = 2)
lines(x = c(262901971,262901971),
      y = c(0.75,1.6),
      col = 'cornflowerblue',
      lwd = 2)
points(x = 262901923,
       y = 1.25,
       col = 'cornflowerblue', 
       bg = 'cornflowerblue',
       cex = cexs*3, 
       pch = 25)
points(x = 262901971,
       y = 1.60,
       col = 'cornflowerblue', 
       bg = 'cornflowerblue',
       cex = cexs*3, 
       pch = 25)
text(x = 262901923,
     y = 1.90, 
     labels = "D576DX", 
     col = 'cornflowerblue',
     cex = cexs)
text(x = 262901971,
     y = 2.2, 
     labels = "L592LX", 
     col = 'cornflowerblue',
     cex = cexs)

## display mutation #3 - 377T1 (heterozygous!)
lines(x = c(262917104,262917104),
      y = c(0.625,1.25),
      col = 'cornflowerblue',
      lwd = 2)
points(x = 262917104,
       y = 1.25,
       col = 'cornflowerblue', 
       bg = 'cornflowerblue', 
       cex = cexs*3, 
       pch = 25)
text(x = 262917104,
     y = 1.55,
     labels = "T1158*",
     col = 'cornflowerblue',
     cex = cexs)

## display mutation #4
lines(x = c(262953328,262953328),
      y = c(0.625,1.25),
      col = 'cornflowerblue',
      lwd = 2)
points(x = 262953328,
       y = 1.25,
       col = 'cornflowerblue', 
       bg = 'cornflowerblue', 
       cex = cexs*3, 
       pch = 25)
text(x = 262953328,
     y = 1.55,
     labels = "A1541*",
     col = 'cornflowerblue',
     cex = cexs)

## display mutation #5
lines(x = c(262968827,262968827),
      y = c(0.75,1.25),
      col = 'cornflowerblue',
      lwd = 2)
points(x = 262968827,
       y = 1.25,
       col = 'cornflowerblue', 
       bg = 'cornflowerblue', 
       cex = cexs*3, 
       pch = 25)
text(x = 262968827,
     y = 1.55, 
     labels = "Q2548*", 
     col = 'cornflowerblue',
     cex = cexs)

## axis
axis(side = 1, 
     at = seq(f = 262880000, t = 262980000, length.out = 5),
     labels = c("262.880", "262.905", "262.930", "262.955", "262.980"),
     cex.axis = 1.8)
mtext(substitute(paste(italic("MGA"), " locus on Chromosome 2 [MBP]")), side = 1, line = 5, cex = 2.5)
dev.off()

## clean up environment
rm(list=ls())
