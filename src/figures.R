theme_slide <- theme_grey(base_size = 16, base_family = "Lato") +
  theme(plot.background = element_rect(fill = "transparent", colour = NA))

###################
### PHENOTYPING ###
###################

rawData <- data.table(read.csv('data/Phenotype_Panicle_corrected.csv',
                              stringsAsFactors = FALSE))

# remove "Echelle" row
rawData <- rawData[!grep("Echelle", file_name)]

rawData[, Accession := unlist(strsplit(file_name, split = "_", fixed = TRUE))[1], by = file_name]
rawData[, Species := plyr::revalue(toupper(Accession),
                        c(B88 = "O. barthii", IR64 = "O. sativa indica",
                          NIP = "O. sativa japonica", TOG5681 = "O. glaberrima",
                          W1654 = "O. rufipogon")),
        by = Accession]
rawData[Accession == "tog5681", Accession := "Tog5681"]
rawData[, Indiv := unlist(strsplit(file_name, split = "_", fixed = TRUE))[2], by = file_name]
rawData[, Panicle := unlist(strsplit(file_name, split = "_", fixed = TRUE))[3], by = file_name]

setkey(rawData, "Species", "Indiv", "Panicle")
pd <- unique(rawData[, .(Species, Indiv, Panicle, TA_nb, Sp_nb)])
ptypes <- ggplot(pd, aes(x = TA_nb, y = Sp_nb)) +
  theme_slide +
  theme(legend.text = element_text(face = "italic")) +
  geom_smooth(method = lm, se = FALSE, colour = "black", size = 1, alpha = 0.5) +
  geom_point(aes(colour = Species), size = 3, alpha = 0.8) +
  scale_color_brewer(palette = "Set1", guide = guide_legend(title = NULL)) +
  ylab("Number of spikelets") + xlab("Number of secondary branches")


rd2 <- as.data.table(read.table("data/OsOgObOrPTRAPdata.txt",
                                sep = "\t",
                                header = TRUE,
                                stringsAsFactors = FALSE,
                                dec = ","))

rd2[, Species := plyr::revalue(Origin, c(Ob = "O. barthii", Os = "O. sativa",
                                     Og = "O. glaberrima", Or = "O. rufipogon"))]

setkey(rd2, "Species", "Bar_Code", "Sowing_nb", "Repet_nb", "Plant_nb", "Panicle_nb")

pd2 <- unique(rd2[, .(Species, Bar_Code, Sowing_nb, Repet_nb, Plant_nb, Panicle_nb, Sp_nb, TA_nb)])
ptypes2 <- ggplot(pd2, aes(x = TA_nb, y = Sp_nb)) +
  theme_slide +
  theme(legend.text = element_text(face = "italic")) +
  geom_point(aes(fill = Species), colour = NA, size = 3, alpha = 0.6, shape = 21) +
  scale_fill_brewer(palette = "Set1", guide = guide_legend(title = NULL)) +
  geom_smooth(method = lm, se = FALSE, colour = "black", size = 0.5) +
  ylab("Number of spikelets") + xlab("Number of secondary branches")


#############
### Mfuzz ###
#############

expressionMatrix <- readRDS('output/mfuzz/expressionMatrix.Rds')
c1 <- readRDS('output/mfuzz/c1.Rds')
memCutoff <- 0.5

# get clusters and membership
cluster <- data.table(id = names(c1$cluster), Cluster = c1$cluster,
                      Membership = apply(c1$membership, 1, max), key = "id")

# get standardised VST counts
exprs <- data.table(Biobase::exprs(expressionMatrix), keep.rownames = TRUE)
exprs[,id := rn][, rn := NULL]
setkey(exprs, "id")

plotData.wide <- exprs[cluster[Membership > memCutoff]]
plotData.wide[, number := length(id), by = Cluster]
plotData.wide[, label := factor(paste0("Cluster ", Cluster, "\n(", number, " genes)"))]

# relevel clusters
centres.wide <- data.table(c1$centers)

# re-order the cluster for the plot
centres.wide[, Cluster := paste("Cluster", 1:nrow(centres.wide))]

# find the changes between RM and PBM and PBM and SM for each cluster
centres.wide[, c("n1n2", "n2n4") :=
               list(PBM - RM,
                    SM - PBM)]
# divide these changes into categories
centres.wide[, c("dn1n2", "dn2n4") :=
               list(c('dec', 'small', 'inc')[cut(n1n2, breaks = c(-Inf, -0.5, 0.5, Inf))],
                    c('dec', 'small', 'inc')[cut(n2n4, breaks = c(-Inf, -1, 1, Inf))])]               

# first, show gradual increase / decrease
centres.wide[dn1n2 == dn2n4, cOrder := c(1,2)[order(RM)]]

# next, big changes in n1n2 but small in n2n4
centres.wide[dn2n4 == 'small', cOrder := c(3,4)[order(RM)]]

# small changes in n1n2, then large change
centres.wide[dn1n2 == 'small', cOrder := c(5,6)[order(RM)]]

# complex patterns 
centres.wide[!dn1n2 == dn2n4 & !dn1n2 == "small" & !dn2n4 == "small",
             cOrder := c(7,8)[order(SM)]]

# order any leftovers on RM
if (any(is.na(centres.wide[, cOrder]))) {
  orderMax <- max(centres.wide[,cOrder], na.rm = TRUE)
  centres.wide[is.na(cOrder), cOrder := c((orderMax + 1):nrow(centres.wide))]
}

# relevel the clusters by cOrder
plotData.wide[, label :=
                factor(label, levels = levels(label)[order(centres.wide[,cOrder])])]

# add label to centres.wide
setkey(centres.wide, 'cOrder')
centres.wide[, label := plotData.wide[, levels(label)]]
centres.wide[, label := factor(label, levels = label)]

# make long
plotData <- reshape2::melt(plotData.wide,
                           id.vars = c("id", "Cluster", "Membership", "label", "number"),
                           variable.name = "Stage",
                           value.name = "Scaled, transformed read counts")

# fix stage label
plotData[, Stage := plyr::mapvalues(Stage, "ePBM.SBM", "ePBM/SBM")]

# set up heatscale
heatscale <- RColorBrewer::brewer.pal(n = 6, name = "YlOrRd")

# add centres to plot
centres <- reshape2::melt(centres.wide, id.vars = 'label',
                          measure.vars = c("RM", "PBM", "ePBM.SBM", "SM"),
                          variable.name = 'Stage',
                          value.name = "Scaled, transformed read counts")
centres[, Stage := plyr::mapvalues(Stage, "ePBM.SBM", "ePBM/SBM")]

# main cluster figure
mfuzz <- ggplot(plotData,
                aes(x = Stage, y = `Scaled, transformed read counts`,
                    colour = Membership, group = id)) +
  theme_slide +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  xlab(NULL) +
  scale_colour_gradientn(colours = heatscale, limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  geom_line(alpha = 0.8) +
  geom_line(data = centres, mapping = aes(group = 1), colour = "black", alpha = 0.5) +
  facet_wrap("label", ncol = 4)

#################
### HYPERGEOM ###
#################

t_hypergeom <- readRDS('output/mfuzz/sigFamilies.Rds')

# fix p-values. 
pFix <- function(x) {
  format.pval(c(0.122, x), digits = 2, eps = 0.001, na.form = "")[-1]
}
# just use a for loop, should use *apply() but too complicated
for (x in colnames(t_hypergeom)) {
  if (grepl("adj", x)) {
    t_hypergeom[,x] <- pFix(t_hypergeom[,x])
  }
}

# remove NAs
t_hypergeom[is.na(t_hypergeom)] <- ""

############
### GSEA ###
############

gsea <- readRDS('output/gsea/gseaTable.Rds')
famCat <- readRDS("data/tfdb/famCat.Rds")

# format some labels
setnames(gsea, old = "Test statistic", new = "Test\nstatistic")
gsea[, Stage := plyr::mapvalues(Stage, "ePBM/SBM", "ePBM/\nSBM")]

# separate by TF / other proteins
setkey(famCat, 'Family')
setkey(gsea, 'rn')
gsea[, Category := famCat[gsea][,Category]]
gsea[, Category := plyr::mapvalues(Category, from = c("TF", "Other"),
                                   to = c("Transcription factors", "Others"))]
gsea[, Category := factor(Category, levels = c("Transcription factors", "Others"))]

heatscale <- rev(RColorBrewer::brewer.pal(6, "RdBu"))

f_gsea <- ggplot(gsea, aes(x = Stage, y = rn, label = padj, fill = `Test\nstatistic`)) +
  xlab(NULL) + ylab(NULL) +
  theme_slide +
  theme(axis.text.x = element_text(vjust = 0.5)) +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  scale_fill_gradientn(colours = heatscale) +
  facet_grid(Category ~ ., scales = "free_y", space = "free_y") +
  geom_raster()
if (gsea[,any(showPval)]){
  f_gsea <- f_gsea + geom_text(data = gsea[which(showPval),], size = 2)
}



############
### ALOG ###
############

# HDZIP4 <- c('LOC_Os02g45250', 'LOC_Os04g48070', 'LOC_Os04g53540', 'LOC_Os08g04190', 'LOC_Os08g08820', 'LOC_Os08g19590', 'LOC_Os09g35760', 'LOC_Os10g42490')
# HDZIP3 <- c('LOC_Os01g48170', 'LOC_Os03g01890', 'LOC_Os03g43930', 'LOC_Os05g48820', 'LOC_Os06g39906', 'LOC_Os10g33960', 'LOC_Os12g41860')
# 
# expTable <- data.table(msuId = HDZIP3, type = 'hb')

library(ggtree)

svpMads <- c("LOC_Os02g52340", "LOC_Os03g08754", "LOC_Os06g11330")
expG1s <- c("LOC_Os02g07030", "LOC_Os06g46030", "LOC_Os10g33780")
ck <- c("LOC_Os01g51210","LOC_Os04g43840","LOC_Os01g40630", "LOC_Os01g10110")

expTable <- rbind(data.table(msuId = svpMads, type = "mads"),
                  data.table(msuId = expG1s, type = "g1l"),
                  data.table(msuId = ck, type = 'ck'))

expTable[, symbol := oryzr::LocToGeneName(msuId)$symbols, by = msuId]

# add expression values
tpm <- data.table(readRDS('output/tpm/tpm.Rds'), keep.rownames = TRUE)
setkey(tpm, "rn")
plotData.wide <- tpm[expTable]

# convert to long
plotData <- reshape2::melt(plotData.wide, id.vars = c("rn", "type", "symbol"),
                           variable.name = "library", value.name = "Expression (TPM)")
setkey(plotData, "rn", "library")

# add expression calls
expGenTT.wide <- data.table(readRDS('output/expressedGenes/expGenTT.Rds'), key = "id")
expGenTT <- reshape2::melt(expGenTT.wide, id.vars = "id", variable = "library", value = "isExpr")
setkey(expGenTT, "id", "library")
plotData <- expGenTT[plotData, .(
  id,
  type,
  library,
  symbol,
  `Expression (TPM)`,
  isExpr)]

# add stage
plotData[, stage := substr(library, start = 1, stop = 2)]
old <- c("n1", "n2", "n3", "n4")
new <- c("RM", "PBM", "ePBM/\nSBM", "SM")
plotData[, stage := factor(plyr::mapvalues(stage, from = old, to = new), levels = new)]

# set up labels
plotData[!is.na(symbol), symbol := paste(symbol, "·", id)]
plotData[is.na(symbol), symbol := id]

# make a plot
cols <- RColorBrewer::brewer.pal(3, "Set1")[c(2,1)]
alogPlot <- function(plotData) {
  return(
    ggplot(plotData, aes(x = stage, y = `Expression (TPM)`, group = symbol,
                         colour = isExpr)) +
      theme_minimal(base_size = 8, base_family = "Helvetica") +
      theme(axis.text.x = element_text(vjust = 0.5),
            strip.text = element_text(face = "italic"),
            plot.title = element_text(hjust = 0, face = "bold"),
            plot.background = element_rect(colour = "black", size = 0.25)) +
      scale_colour_manual(values = cols, guide = FALSE) +
      xlab(NULL) +
      stat_smooth(se = FALSE, colour = "grey", size = 0.5) +
      geom_point(shape = 16, alpha = 0.7, position = position_jitter(height = 0, width = 0.3))
  )
}

f_alogFamily_a <- alogPlot(plotData[type == 'ck']) + ggtitle("a") + facet_wrap(~symbol, ncol = 2)
f_alogFamily_b <- alogPlot(plotData[type == 'g1l']) + ggtitle("b") + facet_wrap(~symbol, ncol = 1)
f_alogFamily_c <- alogPlot(plotData[type == 'mads']) + ggtitle("c") + facet_wrap(~symbol, ncol = 1)

#f_alogFamily <- gridExtra::grid.arrange(f_alogFamily_a, f_alogFamily_b, ncol = 2)

#################
### MADS TREE ###
#################

library(ggtree)

njTree <- readRDS('output/madsComp/clustal/njTree.Rds')
madsPeptides <- readRDS('output/madsComp/clustal/madsPeptides.Rds')
minpcnongap <- readRDS('output/madsComp/clustal/minpcnongap.Rds')
minpcident <- readRDS('output/madsComp/clustal/minpcident.Rds')
minProtLength <- readRDS('output/madsComp/clustal/minProtLength.Rds')
og <- readRDS('output/madsComp/clustal/og.Rds')

# set up expression values for annotation
setkey(madsPeptides, "name")
exprAnnot <- madsPeptides[unique(njTree$tip.label), .(name, log2FoldChange)]

# draw a tree
heatscale <- rev(RColorBrewer::brewer.pal(5, "PuOr"))
sf_madsTree <- ggtree::ggtree(njTree, aes(x = x, y = y, label = label), size = 0.025) +
  xlab(NULL) + ylab(NULL) +
  scale_y_continuous(expand = c(0,1)) +
  theme_minimal(base_size = 8, base_family = "Helvetica") +
  theme(axis.text = element_blank(),
        panel.grid = element_blank(),
        legend.key.size = grid::unit(8, "point"),
        legend.text = element_text(size = 4),
        legend.title = element_text(size = 5),
        legend.position = c(0,0.5),
        legend.justification = c(0,0.5))
# add expression values as an annotation
sf_madsTree <- sf_madsTree %<+% exprAnnot
sf_madsTree <- sf_madsTree +
  scale_fill_gradient2(low = heatscale[1], mid = 'grey90', high = heatscale[5],
                       midpoint = 0, na.value = "white",
                       name = expression(L[2]*"FC")) +
  geom_label(mapping = aes(fill = log2FoldChange), na.rm = TRUE,
             hjust = "left",
             size = 1.5,
             colour = NA,
             label.padding = grid::unit(0.05, "lines"),
             label.r = grid::unit(0.05, "lines")) +
  ggplot2::geom_text(size= 1.5, na.rm = TRUE, hjust = -0.01)

# annotate clades
sf_madsTree <- ggtree::annotation_clade(sf_madsTree, node = 145, "AGL2-like",
                                       offset = 0.15, font.size = 1.5, bar.size = 0.5)
sf_madsTree <- ggtree::annotation_clade(sf_madsTree, node = 159, "AGL6-like",
                                       offset = 0.09, font.size = 1.5, bar.size = 0.5)
sf_madsTree <- ggtree::annotation_clade(sf_madsTree, node = 132, "TM3-like",
                                       offset = 0.1, font.size = 1.5, bar.size = 0.5)
sf_madsTree <- ggtree::annotation_clade(sf_madsTree, node = 168, "AG-like",
                                       offset = 0.09, font.size = 1.5, bar.size = 0.5)
sf_madsTree <- ggtree::annotation_clade(sf_madsTree, node = 165, "AGL12-like",
                                       offset = 0.09, font.size = 1.5, bar.size = 0.5)
sf_madsTree <- ggtree::annotation_clade(sf_madsTree, node = 178, "SQUA-like",
                                       offset = 0.13, font.size = 1.5, bar.size = 0.5)
sf_madsTree <- ggtree::annotation_clade(sf_madsTree, node = 192, "STMADS11-like",
                                       offset = 0.12, font.size = 1.5, bar.size = 0.5)
sf_madsTree <- ggtree::annotation_clade(sf_madsTree, node = 200, "AGL17-like",
                                       offset = 0.14, font.size = 1.5, bar.size = 0.5)
sf_madsTree <- ggtree::annotation_clade(sf_madsTree, node = 208, "FLC-like",
                                       offset = 0.09, font.size = 1.5, bar.size = 0.5)
sf_madsTree <- ggtree::annotation_clade(sf_madsTree, node = 214, "GLO-like",
                                       offset = 0.06, font.size = 1.5, bar.size = 0.5)
sf_madsTree <- ggtree::annotation_clade(sf_madsTree, node = 218, "DEF-like",
                                       offset = 0.08, font.size = 1.5, bar.size = 0.5)
sf_madsTree <- ggtree::annotation_clade(sf_madsTree, node = 114, "MIKC*",
                                       offset = 0.14, font.size = 1.5, bar.size = 0.5)
sf_madsTree <- ggtree::annotation_clade(sf_madsTree, node = 211, "GGM13-like",
                                       offset = 0.025, font.size = 1.5, bar.size = 0.5,
                                       angle = 0, offset.text = 0.04)

# 1-col width=3.150,
# max height=8.661,
# max width = 6.614
# cairo_pdf(filename = "output/madsComp/clustal/tempTree.pdf", width = 3.150,
#    height = 8.661)
#  sf_madsTree
#  dev.off()
# 
# quick print for node labels
# cairo_pdf(filename = "output/madsComp/clustal/nodes.pdf", width = 6.614,
#           height = 8.661)
# ggtree(njTree) + ggplot2::geom_text(aes(label = label), hjust = -0.1, size = 2) + 
#   ggplot2::geom_text(mapping = aes(x = branch, label = node), vjust = -0.3, size = 2)
# 
# dev.off()

################
### HOMEOBOX ###
################

plotData.long <- readRDS("output/homeobox/plotData.long.Rds")
segData <- readRDS("output/homeobox/segData.Rds")

library(ggplot2)
heatscale <- RColorBrewer::brewer.pal(6, "YlOrRd")

f_hb <- ggplot() +
  theme_minimal(base_size = 8, base_family = "Helvetica") +
  theme(
    legend.key.size = grid::unit(8, "points"),
    legend.text	= element_text(size = 5.5),
    axis.text.x = element_text(vjust = 0.5, hjust = 1, angle = 90),
    axis.text.y = element_text(face = "italic", size = 5.5),
    #legend.title = element_text(size = 6),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  xlab(NULL) + ylab(NULL) +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  geom_raster(aes(x = Stage, y = symbol, fill = `Scaled reads`),
              data = plotData.long) +
  scale_fill_gradientn(colours = heatscale) +
  geom_segment(aes(x = 0.255, y = y, xend = 0.255, yend = yend, colour = class),
               data = segData, size = 5) +
  scale_colour_brewer(palette = "Set3", guide = guide_legend(title = NULL))
