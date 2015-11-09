
###################
### STATS TABLE ###
###################

st_libStats <- readRDS('output/quantStats/libStats.Rds')

################
### PCA PLOT ###
################

expressedGenes <- readRDS('output/expressedGenes/expressedGenesAll.Rds')
vst <- readRDS('output/DESeq2/vstAll.Rds')

exprVst <- GenomicRanges::assay(vst)[expressedGenes,]
exprVst <- GenomicRanges::assay(vst)
pca <- prcomp(t(exprVst))
percentVar <- pca$sdev^2/sum(pca$sdev^2)

pcaPlotData <- data.frame(
  label = toupper(rownames(pca$x)),
  PCA1 = pca$x[,1],
  PCA2 = pca$x[,2],
  Stage = GenomicRanges::colData(vst)$stage
)

# row index to help with overplotted labels
pcaPlotData$rIdx <- 1:dim(pcaPlotData)[1]
opIdx <- pcaPlotData[c('n2r3', 'n4r1'),'rIdx']

# format legend
pcaPlotData$Stage <- plyr::mapvalues(pcaPlotData$Stage, from = "ePBM/SBM",
                                     to = "ePBM/\nSBM")

sf_pca <- ggplot(data = pcaPlotData,
                 aes(x = PCA1, y = PCA2, colour = Stage, label = label)) +
  theme_minimal(base_size = 8, base_family = "Helvetica") +
  scale_colour_brewer(palette = "Set1", name = 'Stage') +
  guides(colour=guide_legend(title=NULL)) +
  scale_fill_brewer(palette = "Set1", guide = FALSE) +
  xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
  ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
  xlim(-50, 77) +
  geom_point(shape = 16, size = 2, alpha = 0.7) +
  geom_text(data = pcaPlotData[-opIdx,], size = 2,
            hjust = 1.2, show.legend = FALSE) +
  geom_text(data = pcaPlotData[opIdx,], size = 2,
            hjust = -0.2, show.legend = FALSE)

#######################
### COMPARE IN SITU ###
#######################

compare <- readRDS('output/compare/compare.Rds')
setkey(compare, "zhangRef")

# convert zhangRef to citekey (these papers were mined from zhang paper manually)
citekeys <- structure(c("@Ikeda:2007ja", "@Li:2011es", "@Ren:2013jv", "@Suzaki:2004ib", 
                        "@Chu:2006fw", "@Lee:2012hj", "@Yoshida:2013ff", "@Xue:2008ki", 
                        "@Ashikari:2005eg", "@Yan:2011hw", "@Li:2013iq", "@Yoshida:2012fg", 
                        "@Kurakawa:2007go", "@Lee:2012hj", "@Lee:2007cj", "@Gao:2010iz", 
                        "@IkedaKawakatsu:2012co", "@Horigome:2009gt", "@Lee:2007cj", 
                        "@Jiao:2010ft", "@Miura:2010it", "@Huang:2009kf", "@Kobayashi:2010hn", 
                        "@Komatsu:2003iu", "@Li:2010ks", "@Oikawa:2009cs", "@Tabuchi:2011gx", 
                        "@Zhu:2010je", "@LopezDee:1999ds", "@Yadav:2007im", "@Nagasawa:2003gu", 
                        "@Yun:2013id", "@Duan:2012bt", "@Ohmori:2009kj", "@Dreni:2007jd", 
                        "@Cui:2010kt", "@Fornara:2004bi", "**m.Komatsu2009**"),
                      .Names = c("56", 
                                 "85", "118", "129", "21", "82", "154", "148", "3", "150", "86", 
                                 "153", "78", "82", "83", "48", "59", "54", "83", "68", "100", 
                                 "55", "73", "75", "84", "108", "131", "165", "991", "992", "993", 
                                 "994", "995", "996", "997", "998", "999", "990"))

compare[, Reference := citekeys[as.character(zhangRef)]]

# SI table
st_reviewInSitu <- data.table(reshape2::dcast(compare, msuId + Reference ~ stage,
                                              value.var = 'compare'))

# deal with genes that have the same expression pattern in >1 report
st_reviewInSitu <- st_reviewInSitu[, .(
  Reference = paste(Reference, collapse = ", "))
  , by = c("msuId", "RM", "PBM", "SBM", "SM", "FM")]

st_reviewInSitu[, `Gene symbol` :=
                  oryzr::LocToGeneName(msuId)$symbols,
                by = msuId]
setnames(st_reviewInSitu, old = "msuId", new = "MSU identifier")
setcolorder(st_reviewInSitu, neworder =
              c('Gene symbol', 'MSU identifier', 'RM', 'PBM', 'SBM', 'SM', 'FM',
                'Reference'))
setkey(st_reviewInSitu, "Gene symbol")

# get long tpm data with stage
tpm.wide <- data.table(readRDS('output/tpm/tpm.Rds'), keep.rownames = TRUE)
setnames(tpm.wide, 'rn', "MSU identifier")
tpm <- reshape2::melt(tpm.wide, id.vars = "MSU identifier",
                      variable.name = "library", value.name = "tpm")
colData <- data.table(as.data.frame(GenomicRanges::colData(
  readRDS('output/DESeq2/ddsLrt.Rds'))), keep.rownames = TRUE)
tpm[, stage := colData[rn %in% as.character(library), stage], by = library]

setkey(tpm, "MSU identifier", "stage")

# merge tpm data with compare calls
genes <- compare[, .(msuId, stage, call.z)]
setkey(genes, "msuId", "stage")
genes <- unique(genes)
genes[, stage := plyr::mapvalues(stage, "SBM", "ePBM/SBM")]

# dirty kludge, fix in compareVsInSitu.R
genes <- genes[!stage == "FM"]
setkey(genes, "msuId", "stage")

genTpm <- tpm[genes, .(
  msuId = `MSU identifier`,
  library,
  stage,
  call.z,
  tpm
)]

# add call per library
expGenTT.wide <- data.table(readRDS('output/expressedGenes/expGenTT.Rds'), key = "id")
expGenTT <- reshape2::melt(expGenTT.wide, id.vars = "id",
                           variable.name = "library", value.name = "isExpr")
setkey(expGenTT, "id", "library")
setkey(genTpm, "msuId", "library")

genTpm <- expGenTT[genTpm]
setnames(genTpm, "id", "msuId")

# add replicate number to colour points
genTpm[, Replicate := sub(".*(\\d+)", "\\1", library)]

# format gene name and stage for the plot
genTpm[, name := oryzr::LocToGeneName(msuId)$symbols, by = msuId]
genTpm[, plotLabel := name]
genTpm[!is.na(plotLabel), plotLabel := paste(plotLabel, "·", msuId)]
genTpm[is.na(plotLabel), plotLabel := msuId]
genTpm[, stage := plyr::mapvalues(stage, "ePBM/SBM", "ePBM/\nSBM")]

# order the plots
setkey(genTpm, "name")
genTpm[, plotLabel := factor(plotLabel, levels = unique(plotLabel))]

# add coordinates for shading
genTpm[, c("xmin", "xmax") := .(
  as.numeric(substr(library, 2, 2)) - 0.2,
  as.numeric(substr(library, 2, 2)) + 0.2) ]

# draw it
setkey(genTpm, "plotLabel", "stage")
# sf_isGenesTpm <- ggplot(genTpm, aes(x = stage, y = tpm, colour = Replicate, shape = isExpr)) +
#   theme_minimal(base_size = 8, base_family = "Helvetica") +
#   theme(strip.text = element_text(face = "italic")) +
#   xlab(NULL) +
#   scale_color_brewer(palette = "Set1") +
#   guides(shape = FALSE, size = FALSE) +
#   geom_rect(aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
#             unique(genTpm[call.z == TRUE,]), colour = NA,  fill = "grey90", alpha = 0.75) +
#   stat_smooth(aes(group = plotLabel), method = "loess", se = FALSE, colour = "grey", size = 0.5) + 
#   geom_point(size = 1, alpha = 0.8, position = position_jitter(width = 0.2)) +
#   facet_wrap(~ plotLabel, scales = "free_y", ncol = 4)

# keep a subset of genes
keep <- c("LOC_Os06g45460", "LOC_Os09g26999", "LOC_Os07g42410", "LOC_Os07g47330",
          "LOC_Os10g33780", "LOC_Os03g60430", "LOC_Os01g61480", "LOC_Os03g11614",
          "LOC_Os01g40630", "LOC_Os04g49150", "LOC_Os07g41370", "LOC_Os02g52340",
          "LOC_Os09g32948", "LOC_Os06g11330", "LOC_Os08g41950", "LOC_Os02g45770",
          "LOC_Os03g51690", "LOC_Os03g54170", "LOC_Os04g51000", "LOC_Os07g13170",
          "LOC_Os08g39890")

genTpm <- genTpm[msuId %in% keep]

refTable <- compare[, .(msuId, plotLabel, Reference)]
refTable <- refTable[msuId %in% keep]
setkey(refTable, "msuId", "Reference")
refTable <- unique(refTable)
refTable[, Reference := paste(unlist(Reference), collapse = ", "), by = msuId]
# refTable[, Reference := paste0("\\cite{", paste(Reference, collapse = ","), "}"), by = msuId]
# refTable <- unique(refTable)
# refTable[, Reference := gsub("@", "", Reference, fixed = TRUE)]
setkey(refTable, "plotLabel")
st_refTable <- unique(refTable)
setnames(st_refTable, c("msuId", "plotLabel"), c("MSU name", "CSGNL symbol"))

sf_isGenesTpm <- ggplot() +
  theme_minimal(base_size = 8, base_family = "Helvetica") +
  theme(strip.text = element_text(face = "italic")) +
  xlab(NULL) +
  scale_color_brewer(palette = "Set1") +
  guides(shape = FALSE, size = FALSE) +
  geom_blank(mapping = aes(x = stage, y = tpm), data = genTpm) +
  geom_rect(mapping = aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            data = unique(genTpm[call.z == TRUE,]),
            colour = NA,  fill = "grey90") +
  stat_smooth(mapping = aes(x = stage, y = tpm, group = plotLabel),
              data = genTpm,
              method = "loess", se = FALSE, colour = "black", size = 0.25, alpha = 0.6) + 
  geom_point(mapping = aes(x = stage, y = tpm, colour = Replicate, shape = isExpr),
             data = genTpm,
             size = 1, alpha = 0.6) +
  facet_wrap(~ plotLabel, scales = "free_y", ncol = 2)


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
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  xlab(NULL) +
  scale_colour_gradientn(colours = heatscale, limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  geom_line(alpha = 0.8) +
  geom_line(data = centres, mapping = aes(group = 1), colour = "black", alpha = 0.5) +
  facet_wrap("label", ncol = 4)

# centroid dis vs. c (for SI)
centroids <- readRDS('output/mfuzz/centroids.Rds')
sf_mfuzzCentroids <- ggplot(centroids, aes(x = x, y = y)) +
  theme_minimal(base_size = 8, base_family = "Helvetica") +
  theme(plot.title = element_text(hjust = 0, face = "bold"),
        plot.background = element_rect(colour = "black", size = 0.25)) +
  xlab(expression(Cluster~number~"("*italic(c)*")")) +
  ylab("Minimum centroid distance") +
  ggtitle("a") +
  stat_smooth(method = loess, se = FALSE, size = 0.5) +
  geom_point(size = 1, alpha = 0.5)

# MDS for SI
vg.mds <- readRDS('output/mfuzz/vg.mds.Rds')
sf_mfuzzPca <- ggplot(vg.mds, aes(x = MDS1, y=MDS2, colour = factor(cluster))) +
  theme_minimal(base_size = 8, base_family = "Helvetica") +
  theme(legend.key.size = grid::unit(8, "point"),
        plot.title = element_text(hjust = 0, face = "bold"),
        plot.background = element_rect(colour = "black", size = 0.25)) +
  #coord_fixed(ratio = 1) +
  ggtitle("b") +
  geom_point(aes(size=max.membership),alpha=0.5, shape = 16) +
  scale_colour_brewer(palette = "Set1", name = "Cluster") +
  scale_size_area(guide = FALSE, max_size = 1)

sf_mfuzzPcaCentroids <- gridExtra::arrangeGrob(
  sf_mfuzzCentroids,
  sf_mfuzzPca,
  padding = unit(0, "lines"), ncol = 1)

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
        legend.key.size = unit(8, "point"),
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
             label.padding = unit(0.05, "lines"),
             label.r = unit(0.05, "lines")) +
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

detach("package:ggtree", unload=TRUE)

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