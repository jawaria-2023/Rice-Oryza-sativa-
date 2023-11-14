# Rice-Oryza-sativa-
The analysis of microarrays on laser-microdissected tissues reveals the biosynthesis of suberin in the outer regions of roots during the development of a barrier to radial oxygen loss in rice (Oryza sativa).

Organism:	Oryza sativa
Experiment type: 	Expression profiling by array
Summary: For plants growing in soil with too much water, it's important for their roots to get enough air inside. A special barrier, called the radial oxygen loss (ROL) barrier, helps transport oxygen over long distances from the plant's stem to the tip of its roots. This higher oxygen concentration at the root tip allows the plant to grow even in soil without much oxygen. The ROL barrier forms in the outer part of the roots (OPR). Some substances, like suberin and/or lignin, found in the cell walls of the roots, are believed to help create this barrier, but it's not clear which one is the main contributor.

This study looked at the genes that are active when the ROL barrier is forming in rice roots. Researchers used a precise method called laser microdissection to isolate specific tissues in the outer part of the roots, and then they analyzed the genes in these tissues using a method called microarray. They found that certain genes associated with suberin production were strongly turned on during the barrier formation, while genes linked to lignin production did not show significant changes. An analysis of the promoters of these active genes revealed elements associated with specific transcription factors (WRKY, AP2/ERF, NAC, bZIP, MYB, CBT/DREB, and MADS). These factors were particularly connected to the expression of genes containing WRKY, AP2, and MYB domains.

When the researchers checked the activity of specific genes related to suberin production using a method called semiquantitative reverse-transcription PCR, they confirmed that these genes were highly active during the formation of the ROL barrier. In conclusion, these findings suggest that suberin is a major component of the ROL barrier in rice roots.

############################### R SCRIPT######################################
#   Data plots for selected GEO samples
library(GEOquery)
library(limma)
library(umap)

# load series and platform data from GEO

gset <- getGEO("GSE58804", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL6864", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

ex <- exprs(gset)
# log2 transform
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
          (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
  ex <- log2(ex) }

# box-and-whisker plot
par(mar=c(7,4,2,1))
title <- paste ("GSE58804", "/", annotation(gset), sep ="")
boxplot(ex, boxwex=0.7, notch=T, main=title, outline=FALSE, las=2)

# expression value distribution plot
par(mar=c(4,4,2,1))
title <- paste ("GSE58804", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex, main=title, legend=F)

# mean-variance trend
ex <- na.omit(ex) # eliminate rows with NAs
plotSA(lmFit(ex), main="Mean variance trend, GSE58804")

# UMAP plot (multi-dimensional scaling)
ex <- ex[!duplicated(ex), ]  # remove duplicates
ump <- umap(t(ex), n_neighbors = 8, random_state = 123)
plot(ump$layout, main="UMAP plot, nbrs=8", xlab="", ylab="", pch=20, cex=1.5)
library("maptools")  # point labels without overlaps
pointLabel(ump$layout, labels = rownames(ump$layout), method="SANN", cex=0.6)

