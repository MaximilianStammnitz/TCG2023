## Figure 6E
## 'The evolution of two transmissible cancers in Tasmanian devils'
## Stammnitz et al., 2022
## by Kevin Gori: kcg25@cam.ac.uk

library(readxl)
library(data.table)
library(plotrix)
library(bit64)
library(rectanglePacking)

## set input path(s)
setwd('/Tables')


# DFT1/DFT2 copy number event display #
#######################################

### Load data:
## chrlengths:
# A table giving chromosome sizes of the Devil reference genome.
# Also provides genomic start and end coordinates of the chromosomes
# Finally, the GenomicOffset can be added to within-chromosome positions
# to map them to genomic position
chrlengths <- data.table(CHROM=c(1:6, "X", "Y"),
                         LENGTH=c(716413629, 662751787, 611347268, 464895054, 288121652, 254895979, 83081154, 130564))
chrlengths[, GenomeStart := 0]
chrlengths[, GenomeEnd := cumsum(LENGTH)]
chrlengths[, GenomeStart := 1 + c(0, shift(GenomeEnd)[2:.N])]
chrlengths[, GenomeOffset := GenomeStart - 1]

## centromere:
# A table giving the centromeric position for chromosomes 1-X
centromere <- data.table(CHROM=c(1,2,3,4,5,6,"X"),
                         START=c(231517792,331316704,189436634,117838506,98635937,108814554,7949349),
                         END=c(280330846,507217211,229173042,299488986,130085044,172768340,14576901),
                         POS=c(286700000,341400000,192000000,293000000,101800000,122000000,11000000),
                         CHREND = c(716413629, 662751787, 611347268, 464895054, 288121652, 254895979, 83081154),
                         INVERT = c(FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE))

## cnvtable:
# The Copy number events inferred for DFT1 and DFT2
cnvtable.DFT1 <- as.matrix(read_xlsx("Table-S6.xlsx", sheet = 2))
colnames(cnvtable.DFT1) <- as.character(cnvtable.DFT1[2,])
cnvtable.DFT1 <- cnvtable.DFT1[-c(1:2),]
cnvtable.DFT2 <- as.matrix(read_xlsx("Table-S6.xlsx", sheet = 3))
colnames(cnvtable.DFT2) <- as.character(cnvtable.DFT2[2,])
cnvtable.DFT2 <- cnvtable.DFT2[-c(1:2),]
cnvtable <- as.data.table(rbind(cbind("CHROM" = cnvtable.DFT1[,2], 
                                      "clone" = 'DFT1', 
                                      "ID" = as.numeric(cnvtable.DFT1[,1]), 
                                      cnvtable.DFT1[,3:5], 
                                      "EVENT" = cnvtable.DFT1[,10], 
                                      "N" = cnvtable.DFT1[,9]),
                                cbind("CHROM" = cnvtable.DFT2[,2], 
                                      "clone" = 'DFT2', 
                                      "ID" = as.numeric(cnvtable.DFT2[,1]), 
                                      cnvtable.DFT2[,3:5], 
                                      "EVENT" = cnvtable.DFT2[,10],
                                      "N" = cnvtable.DFT2[,9])))

cnvtable[, EventType := ifelse(grepl("GAIN", EVENT), "gain", "loss")]
cnvtable[, N := as.integer(N)]
cnvtable[, START := as.integer(START)]
cnvtable[, END := as.integer(END)]

# Convert positions to genomic coordinates (use integer64 type to accommodate large numbers)
cnvtable[chrlengths, GenomeStartPos := START + GenomeOffset, on = "CHROM"]
cnvtable[chrlengths, GenomeEndPos := END + GenomeOffset, on = "CHROM"]
centromere[chrlengths, GenomePos := POS + GenomeOffset, on = "CHROM"]

# Assign a uniqueID to each CNV
cnvtable[, oldID := ID]
cnvtable[, ID := paste(clone, CHROM, ID, sep=".")]
stopifnot(cnvtable[, uniqueN(ID) == .N])

# Sort the table (shouldn't strictly be necessary, as the layer assignment code sorts its input,
# but can't hurt)
setorder(cnvtable, CHROM, clone, EventType, GenomeStartPos, GenomeEndPos)

# Assign plot layers to each CNV
for (chrom_ in cnvtable[, sort(unique(CHROM))]) {
    for (clone_ in cnvtable[, sort(unique(clone))]) {
        for (event_type_ in cnvtable[, sort(unique(EventType))]) {
            chrom_table <- cnvtable[CHROM == chrom_ & clone == clone_ & EventType == event_type_]
            rectangle_packing_cpp(chrom_table)
            cnvtable[chrom_table, layer := i.layer, on = c("CHROM", "clone", "EventType", "ID")]
        }
    }
}

get_plot_dimensions <- function(cnv_table) {
    ymin <- cnv_table[EventType=="loss", c(-max(layer))]
    ymax <- cnv_table[EventType=="gain", c(max(layer))]
    list(xlim = cnv_table[, c(min(GenomeStartPos), max(GenomeEndPos))],
         ylim = c(ymin, ymax))
}
draw.chromosome <- function(left, right, middle, bottom, top, centrosize, ...) {
    xs <- c(left + centrosize/2, middle - centrosize/2, middle, middle + centrosize/2, right - centrosize/2,
            right - centrosize/2, middle + centrosize/2, middle, middle - centrosize/2, left + centrosize/2)
    ys <- c(top, top, (3 * top + bottom) / 4, top, top,
            bottom, bottom, (top + 3 * bottom) / 4, bottom, bottom)
    polygon(xs, ys, ...)
    plotrix::draw.ellipse(c(left+centrosize/2, right-centrosize/2), rep((top+bottom)/2, 2),
                          centrosize/2, (top-bottom)/2,
                          ...)
}
plot_chrom_map <- function(cnv_table, chr_lengths, centromeres, plot_options = NULL) {
    default_plot_options <- list(
        chromosome_height = 0.9,
        cnv_height = 0.7,
        col = "cornflowerblue",
        border = NA,
        title = "Chromosome Map",
        centromere_col = "white",
        centromere_cex = 2,
        show_cnv_id = FALSE
    )
    
    if (!is.null(plot_options)) {
        stopifnot(is.list(plot_options))
        for (option in names(plot_options)) {
            default_plot_options[[option]] <- plot_options[[option]]
        }
    }
    
    plot_options <- default_plot_options
    
    dims <- get_plot_dimensions(cnv_table)
    if ("ylim" %in% names(plot_options)) {
        dims[["ylim"]] <- plot_options[["ylim"]]
    }
    plot(NA, xlim = dims$xlim, ylim = dims$ylim, main = plot_options$title,
         ylab = NA, xlab = NA, xaxt = "n", yaxt = "n")
    
    # Useful to know the extents of the plotting window
    window_xrange <- par()$usr[1:2]
    window_yrange <- par()$usr[3:4]
    
    # Sort out the y-axis
    y_axis_limits <- round(2*dims$ylim, -1)/2
    y_axis_major_ticks <- seq(from=y_axis_limits[1], to=y_axis_limits[2], by = 5)
    y_axis_minor_ticks <- seq(from=y_axis_limits[1], to=y_axis_limits[2], by = 1)
    axis(2, at = y_axis_major_ticks, labels = abs(y_axis_major_ticks), las = 2)
    print(par()$tck)
    axis(2, at = y_axis_minor_ticks, labels = NA, tck=-0.01)
    
    mtext(c("Loss depth", "Gain depth"), 2, line=3, at = window_yrange/2, adj=c(0.5, 0.5))
    mtext(chr_lengths[CHROM != "Y", CHROM], 1, line = 0,
          at = chr_lengths[CHROM != "Y", (GenomeStart + GenomeEnd) / 2],
          adj = 0.5)
    
    # Plot background rectangles for chromosomes
    rect(window_xrange[1], window_yrange[1], window_xrange[2], window_yrange[2],
         col = "grey95", border = NA, lwd = 2)
    chr_lengths[, rect(GenomeStart, window_yrange[1], GenomeEnd, window_yrange[2],
                       col = c("grey90", "grey95"), border = "NA")]
    # Plot rectangles for the chromosomes
    # chr_lengths[, rect(GenomeStart, -plot_options$chromosome_height/2, GenomeEnd, plot_options$chromosome_height/2, col = c("grey15", "grey20"), border=NA)]
    # centromeres[, segments(GenomePos, -plot_options$chromosome_height/2, GenomePos, plot_options$chromosome_height/2, col = plot_options$centromere_col, lwd = 3, lend = "butt")]
    # chr_lengths[, rect(GenomeStart, -plot_options$chromosome_height/2, GenomeEnd, plot_options$chromosome_height/2, col = NA, border = c("grey90", "grey95"))]
    for (i in 1:7) {
        chrom_ <- chr_lengths[i, CHROM]
        draw.chromosome(left = chr_lengths[CHROM == chrom_, GenomeStart+3e6],
                        right = chr_lengths[CHROM == chrom_, GenomeEnd-3e6],
                        middle = centromeres[CHROM == chrom_, GenomePos],
                        bottom = -0.42, top = 0.42, centrosize = ifelse(chrom_ == "X", 1.0e7, 2.0e7),
                        col = "white", border = "black")
        draw.chromosome(left = chr_lengths[CHROM == chrom_, GenomeStart+3e6],
                        right = chr_lengths[CHROM == chrom_, GenomeEnd-3e6],
                        middle = centromeres[CHROM == chrom_, GenomePos],
                        bottom = -0.4, top = 0.4, centrosize = ifelse(chrom_ == "X", 1.0e7, 2.0e7),
                        col = c("grey60", "grey55")[i %% 2 + 1], border = NA)
    }
    
    cnv_table[EventType == "gain",
              rect(xleft = GenomeStartPos,
                   ybottom = layer - plot_options$cnv_height/2,
                   xright = GenomeEndPos,
                   ytop = layer + plot_options$cnv_height/2,
                   col = plot_options$col, border = plot_options$border, lwd = 1)]
    cnv_table[EventType == "loss",
              rect(xleft = GenomeStartPos,
                   ybottom = -(layer - plot_options$cnv_height/2),
                   xright = GenomeEndPos,
                   ytop = -(layer + plot_options$cnv_height/2),
                   col = plot_options$col, border = plot_options$border, lwd = 1)]
    if (plot_options$show_cnv_id) {
        cnv_table[, text(x = (GenomeStartPos + GenomeEndPos)/2, y = ifelse(EventType=="gain", 1, -1) * (layer),
                         labels = oldID, adj = c(0.5, 0.5), cex = 0.6)]
    }
    centromeres[, points(GenomePos, rep(0, .N), col = plot_options$centromere_col, pch = 20, cex = plot_options$centromere_cex)]
    
}

pdf("Figure6E_DFT1_DFT2_CNV_event_map.pdf", width=11.7, height=8.3)
par(mfrow = c(2,1), mar = c(4,5,3,3))
plot_chrom_map(cnvtable[clone == "DFT1"], chrlengths, centromere,
               plot_options = list(col = "cornflowerblue",
                                   ylim = c(-19, 20),
                                   centromere_col = "black",
                                   centromere_cex = 0.5,
                                   title = "DFT1"))
plot_chrom_map(cnv_table = cnvtable[clone == "DFT2"],
               chr_lengths = chrlengths,
               centromeres = centromere,
               plot_options = list(col = "red",
                                   ylim = c(-19, 20),
                                   centromere_col = "black",
                                   centromere_cex = 0.5,
                                   title = "DFT2"))
dev.off()

## clean up environment
rm(list=ls())
