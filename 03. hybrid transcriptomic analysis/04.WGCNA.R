## R: source /GS01/project/pengms_group/pengms20t1/dir.xumm/bin/python3.9/bin/activate SQANTI3.env

setwd("/GS01/project/pengms_group/pengms20t1/dir.xumm/Hybrid_assembly/mule_duck/WGCNA/res1")

set.seed(20230803)
rm(list = ls())

#### load packages
library(dplyr)

#### load raw read counts matrix
raw.matrix <- read.table("Mule.Hinny.ducks.AL.OV.counts.txt",header = T, row.names = 1)
raw.matrix <- select(raw.matrix, -c("Hin.OV.27", "Mul.OV.22", "Hin.AL.103"))

#### filter raw read counts matrix
#### removing all features that have a count of less than 10 in more than 90% of the samples
temp <- apply(raw.matrix[,colnames(raw.matrix)] >= 10, 1, sum) >= 4
table(temp)
filter.mat <- raw.matrix[temp,]

#### normalized filtered read counts matrix by varianceStabilizingTransformation() function in DESeq2 packages
library(DESeq2)    ## version: 1.34.0
species <- factor(c(rep("Hinny", 18), rep("Mule",19)))
feed <- factor(c(rep("AL", 9), rep("OV", 9), c(rep("AL",10), rep("OV", 9))))
colData <- data.frame(row.names=colnames(filter.mat), species, feed)
dds <- DESeqDataSetFromMatrix(floor(filter.mat), colData = colData, design= ~ species + feed )
vsd.mat <- varianceStabilizingTransformation(dds, blind = FALSE)


#### load packages
library(WGCNA)                 ## version: 1.72.1
options(stringsAsFactors = FALSE)
enableWGCNAThreads(nThreads = 3)                 ## enable multiple threads

#### 1.prepare input file for WGCNA
datExpr <- data.frame(t((assay(vsd.mat))))
datTraits <- data.frame(sampleID = rownames(colData), species = colData$species, feed = colData$feed, state = paste(colData$species, colData$feed, sep = "."))
rownames(datTraits) <- rownames(colData)

## match order of samples in datExpr with those in datTraits
sampleNames <- rownames(datExpr)
traitRows <- match(sampleNames, datTraits$sampleID)
rownames(datTraits) <- datTraits[traitRows, 1]

##save data
save(datExpr, datTraits, file = "WGCNA-input-data.RData")


#### 2.choosing the best beta value (soft-thresholding power)
powers <- c(c(1:20), seq(from = 20, to = 60, by = 2))

## Call the network topology analysis function
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, networkType = "signed", corFnc = "bicor")

## Plot the results
pdf(file="beta-value-choose.pdf", height = 7.5, width = 15)
par(mfrow = c(1,2))
cex1 = 0.9

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit,signed R^2",type = "n", main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,cex=cex1,col="red")

# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()


#### 3.one-step network construction and module detection
## gene number in datExpr must be less than maxBlockSize value, otherwise blockwiseModules will split the data set into several blocks
## strongly recommend using the argument maxPOutliers = 0.05 or 0.10 whenever the biweight midcorrelation is used
## this step take much time and require bigger memories
net <- blockwiseModules(datExpr, power = sft$powerEstimate, minModuleSize = 100, maxBlockSize = 30000,
                        deepSplit = 4, networkType = "signed", corType = "bicor", maxPOutliers = 0.1,
                        reassignThreshold = 0, mergeCutHeight = 0.2, detectCutHeight = 0.9999,
                        pamStage = FALSE, pamRespectsDendro = TRUE,
                        numericLabels = TRUE, saveTOMs = FALSE, verbose = 3)

## see identified modules and number of genes in each modules
## the lable 0 is reserved for genes outside of all modules
table(net$colors)

## module visualization
## Plot the dendrogram and the module colors underneath
## assign all of the gene to their corresponding module 
## hclust for the genes
# Convert labels to colors for plotting
mergedColors <- labels2colors(net$colors)
table(mergedColors)
pdf(file = "clustering-dendrogram-genes-modules.pdf", height = 10, width = 15)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()

## save data
save(sft, net, file = "networkConstruction-auto.RData")


#### (optional step): cluster samples
#### visualize the relationship between samples and assigned modules
#### from the plot, outliers (few sample cant cluster with other samples under same condition) can be removed
## samples clustering
datExpr_tree <- hclust(dist(datExpr), method = "average")
par(mar = c(0,5,2,0))
plot(datExpr_tree, main = "Sample clustering", sub="", xlab="", cex.lab = 2, cex.axis = 1, cex.main = 1, cex.lab = 1)

pdf(file = "clustering-samples-modules.pdf", height = 9, width = 10)
## module visualization
sample_colors <- numbers2colors(as.numeric(factor(datTraits$state)), colors = c("blue","red","green","yellow"), signed = FALSE)
par(mar = c(1,4,3,1),cex=0.8)
plotDendroAndColors(datExpr_tree, sample_colors, groupLabels = colnames(sample), cex.dendroLabels = 0.8, marAll = c(1, 4, 3, 1), cex.rowText = 0.01, main = "Sample dendrogram and trait heatmap")
dev.off()


#### 4. relating modules to external traits
## a. quantifying module-trait associations
# construct trait matrix for discrete variable
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)
design <- model.matrix(~0 + datTraits$state)                 # design: sample-trait matrix
colnames(design) <- levels(as.factor(datTraits$state))
rownames(design) <- rownames(datTraits)
moduleColors <- labels2colors(net$colors)

# construct module-trait matrix
MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)                                           # MEs: sample-module matrix
moduleTraitCor <- cor(MEs, design , use = "p")                                # moduleTraitCor: module-trait matrix
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

# plot heatmap
pdf(file = "Module-trait-relationships.pdf", width = 9, height = 9)
par(mar = c(4, 7, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor, xLabels = colnames(design), yLabels = names(MEs), ySymbols = names(MEs), colorLabels = FALSE, colors = greenWhiteRed(50),textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 1, zlim = c(-1,1), main = paste("Module-trait relationships"))
dev.off()

# return all the genes in each module
modNames <- substring(names(MEs), 3)[-length(names(MEs))]
for (mymod in modNames){
	interest.gene <- data.frame(gene = names(datExpr)[moduleColors == mymod])
	outputfile <- paste(mymod,"interest.gene.txt",sep = ".")
	write.table(interest.gene, outputfile, row.names = F, quote = F)
}


#### 5. export gene network to external visualization software
# Recalculate topological overlap
TOM <- TOMsimilarityFromExpr(datExpr, power = sft$powerEstimate)      # corelation matrix of all the genes

# save data
save(TOM, file = "exporting-cytoscape-data.RData")
