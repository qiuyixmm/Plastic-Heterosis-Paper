## load packages
library(DESeq2)
library(tidyverse)
library(ggtext)
library(ggalluvial)
library(cowplot)

set.seed(20250621)

###################################################### identify inheritable modes in F1 under ab libitum feeding condition
#### load data
## remove sample: Pek.AL101, Mul.AL.102, Hin.AL.119
liver.count <- read.table("Pekingduck.Muscovyduck.Muleduck.Hinnyduck.AL.counts.txt", header = T, sep = "\t")
liver.count <- liver.count[,c(1,3:20, 22:30, 32:44, 46:54, 56:60)]
rownames(liver.count) <- liver.count$Geneid
liver.count <- liver.count[,-1]


#### generate F1 gene expression matrix
liver.count$Mul.AL.106 <- liver.count$Pek.Mul.AL.106 + liver.count$Mus.Mul.AL.106
liver.count$Mul.AL.110 <- liver.count$Pek.Mul.AL.110 + liver.count$Mus.Mul.AL.110
liver.count$Mul.AL.114 <- liver.count$Pek.Mul.AL.114 + liver.count$Mus.Mul.AL.114
liver.count$Mul.AL.118 <- liver.count$Pek.Mul.AL.118 + liver.count$Mus.Mul.AL.118
liver.count$Mul.AL.82 <- liver.count$Pek.Mul.AL.82 + liver.count$Mus.Mul.AL.82
liver.count$Mul.AL.86 <- liver.count$Pek.Mul.AL.86 + liver.count$Mus.Mul.AL.86
liver.count$Mul.AL.90 <- liver.count$Pek.Mul.AL.90 + liver.count$Mus.Mul.AL.90
liver.count$Mul.AL.94 <- liver.count$Pek.Mul.AL.94 + liver.count$Mus.Mul.AL.94
liver.count$Mul.AL.98 <- liver.count$Pek.Mul.AL.98 + liver.count$Mus.Mul.AL.98

liver.count$Hin.AL.103 <- liver.count$Pek.Hin.AL.103 + liver.count$Mus.Hin.AL.103
liver.count$Hin.AL.107 <- liver.count$Pek.Hin.AL.107 + liver.count$Mus.Hin.AL.107
liver.count$Hin.AL.111 <- liver.count$Pek.Hin.AL.111 + liver.count$Mus.Hin.AL.111
liver.count$Hin.AL.115 <- liver.count$Pek.Hin.AL.115 + liver.count$Mus.Hin.AL.115
liver.count$Hin.AL.83 <- liver.count$Pek.Hin.AL.83 + liver.count$Mus.Hin.AL.83
liver.count$Hin.AL.87 <- liver.count$Pek.Hin.AL.87 + liver.count$Mus.Hin.AL.87
liver.count$Hin.AL.91 <- liver.count$Pek.Hin.AL.91 + liver.count$Mus.Hin.AL.91
liver.count$Hin.AL.95 <- liver.count$Pek.Hin.AL.95 + liver.count$Mus.Hin.AL.95
liver.count$Hin.AL.99 <- liver.count$Pek.Hin.AL.99 + liver.count$Mus.Hin.AL.99

#### generate parental species and F1 gene expression matrix
parents.F1.reads <- liver.count[, c(1:18, 55:72)]      # 13,417 genes
parents.F1.reads.filter <- parents.F1.reads[apply(parents.F1.reads, 1, min) >= 20,]  # remove genes with less than 20 total (8,139 genes)


#### DESeq2 for gene expression comparison among Muscovy, Peking, and Mule ducks
Peking.Muscovy.mule.duck <- parents.F1.reads.filter[, c(1:27)]
condition <- factor(c(rep("Pekingduck", 9), rep("Muscovyduck", 9), rep("Muleduck", 9)))
colData <- data.frame(row.names = colnames(Peking.Muscovy.mule.duck), condition)

Peking.Muscovy.mule.duck.dds <- DESeqDataSetFromMatrix(Peking.Muscovy.mule.duck, colData = colData, design= ~ condition)
Peking.Muscovy.mule.duck.dds <- DESeq(Peking.Muscovy.mule.duck.dds)

## Peking ducks vs. mule ducks
Peking.mule.duck.res <- results(Peking.Muscovy.mule.duck.dds, contrast = c("condition", "Pekingduck", "Muleduck"))
Peking.mule.fc <- data.frame(GeneName = rownames(as.data.frame(Peking.mule.duck.res)), Pek.Mul.fc = Peking.mule.duck.res$log2FoldChange)

# downregulated: 244; upregualted: 331; total DEGs: 575
Peking.mule.duck.diff <- subset(Peking.mule.duck.res, padj < 0.01 & abs(log2FoldChange) > 1)
write.csv(Peking.mule.duck.diff, "Peking.mule.duck.AL.DEGs.csv", row.names = T, quote = T)


## Muscovy ducks vs. mule ducks
Muscovy.mule.duck.res <- results(Peking.Muscovy.mule.duck.dds, contrast = c("condition", "Muscovyduck", "Muleduck"))
Muscovy.mule.fc <- data.frame(GeneName = rownames(as.data.frame(Muscovy.mule.duck.res)), Mus.Mul.fc = Muscovy.mule.duck.res$log2FoldChange)

#  downregulated: 470; upregualted: 936; total DEGs: 1406
Muscovy.mule.duck.diff <- subset(Muscovy.mule.duck.res, padj < 0.01 & abs(log2FoldChange) > 1)
write.csv(Muscovy.mule.duck.diff, "Muscovy.mule.duck.AL.DEGs.csv", row.names = T, quote = T)


## merge results
Peking.Muscovy.muleduck.fc <- inner_join(Peking.mule.fc, Muscovy.mule.fc, by = "GeneName")
  

#### heritable mode identification
Peking.Muscovy.muleduck.overdominant <- subset(Peking.Muscovy.muleduck.fc, Mus.Mul.fc < -log2(1.25) & Pek.Mul.fc < -log2(1.25))
Peking.Muscovy.muleduck.overdominant$mode <- "overdominant"

Peking.Muscovy.muleduck.underdominant <- subset(Peking.Muscovy.muleduck.fc, Mus.Mul.fc > log2(1.25) & Pek.Mul.fc > log2(1.25))
Peking.Muscovy.muleduck.underdominant$mode <- "underdominant"

Peking.Muscovy.muleduck.conserved <- subset(Peking.Muscovy.muleduck.fc, abs(Mus.Mul.fc) < log2(1.25) & abs(Pek.Mul.fc) < log2(1.25))
Peking.Muscovy.muleduck.conserved$mode <- "conserved"

Peking.Muscovy.muleduck.additive <- subset(Peking.Muscovy.muleduck.fc, (Pek.Mul.fc > log2(1.25) & Mus.Mul.fc < -log2(1.25)) | (Pek.Mul.fc < -log2(1.25) & Mus.Mul.fc > log2(1.25)))
Peking.Muscovy.muleduck.additive$mode <- "additive"

Peking.Muscovy.muleduck.M.dominant <- subset(Peking.Muscovy.muleduck.fc, (abs(Pek.Mul.fc) > log2(1.25) & abs(Mus.Mul.fc) < log2(1.25)))
Peking.Muscovy.muleduck.M.dominant$mode <- "Muscovy.dominant"

Peking.Muscovy.muleduck.P.dominant <- subset(Peking.Muscovy.muleduck.fc, (abs(Mus.Mul.fc) > log2(1.25) & abs(Pek.Mul.fc) < log2(1.25)))
Peking.Muscovy.muleduck.P.dominant$mode <- "Peking.dominant"

# Muscovy.dominant: 1,072; Peking.dominant: 2,383
# overdominant: 121; underdominant: 487
# additive: 2,211; conserved: 1,865
Peking.Muscovy.muleduck.H.mode <- rbind(Peking.Muscovy.muleduck.P.dominant, Peking.Muscovy.muleduck.M.dominant,
                                        Peking.Muscovy.muleduck.overdominant, Peking.Muscovy.muleduck.underdominant,
                                        Peking.Muscovy.muleduck.additive, Peking.Muscovy.muleduck.conserved)

write.csv(Peking.Muscovy.muleduck.H.mode, "Peking.Muscovy.muleduck.AL.heritable.mode.csv", row.names = F, quote = F)
write.table(Peking.Muscovy.muleduck.H.mode, "Peking.Muscovy.muleduck.AL.heritable.mode.txt", row.names = F, quote = F, sep = "\t")



#### DESeq2 for gene expression comparison among Muscovy, Peking, and Hinny ducks
Peking.Muscovy.hinny.duck <- parents.F1.reads.filter[, c(1:18, 28:36)]
condition <- factor(c(rep("Pekingduck", 9), rep("Muscovyduck", 9), rep("Hinnyduck", 9)))
colData <- data.frame(row.names = colnames(Peking.Muscovy.hinny.duck), condition)

Peking.Muscovy.hinny.duck.dds <- DESeqDataSetFromMatrix(Peking.Muscovy.hinny.duck, colData = colData, design= ~ condition)
Peking.Muscovy.hinny.duck.dds <- DESeq(Peking.Muscovy.hinny.duck.dds)

## Peking ducks vs. hinny ducks
Peking.hinny.duck.res <- results(Peking.Muscovy.hinny.duck.dds, contrast = c("condition", "Pekingduck", "Hinnyduck"))
Peking.hinny.fc <- data.frame(GeneName = rownames(as.data.frame(Peking.hinny.duck.res)), Pek.Hin.fc = Peking.hinny.duck.res$log2FoldChange)

# downregulated: 238; upregualted: 338; total DEGs: 576
Peking.hinny.duck.diff <- subset(Peking.hinny.duck.res, padj < 0.01 & abs(log2FoldChange) > 1)
write.csv(Peking.hinny.duck.diff, "Peking.hinny.duck.AL.DEGs.csv", row.names = T, quote = T)


## Muscovy ducks vs. hinny ducks
Muscovy.hinny.duck.res <- results(Peking.Muscovy.hinny.duck.dds, contrast = c("condition", "Muscovyduck", "Hinnyduck"))
Muscovy.hinny.fc <- data.frame(GeneName = rownames(as.data.frame(Muscovy.hinny.duck.res)), Mus.Hin.fc = Muscovy.hinny.duck.res$log2FoldChange)

#  downregulated: 432; upregualted: 914; total DEGs: 1346
Muscovy.hinny.duck.diff <- subset(Muscovy.hinny.duck.res, padj < 0.01 & abs(log2FoldChange) > 1)
write.csv(Muscovy.hinny.duck.diff, "Muscovy.hinny.duck.AL.DEGs.csv", row.names = T, quote = T)


## merge results
Peking.Muscovy.hinnyduck.fc <- inner_join(Peking.hinny.fc, Muscovy.hinny.fc, by = "GeneName")


#### heritable mode identification
Peking.Muscovy.hinnyduck.overdominant <- subset(Peking.Muscovy.hinnyduck.fc, Mus.Hin.fc < -log2(1.25) & Pek.Hin.fc < -log2(1.25))
Peking.Muscovy.hinnyduck.overdominant$mode <- "overdominant"

Peking.Muscovy.hinnyduck.underdominant <- subset(Peking.Muscovy.hinnyduck.fc, Mus.Hin.fc > log2(1.25) & Pek.Hin.fc > log2(1.25))
Peking.Muscovy.hinnyduck.underdominant$mode <- "underdominant"

Peking.Muscovy.hinnyduck.conserved <- subset(Peking.Muscovy.hinnyduck.fc, abs(Mus.Hin.fc) < log2(1.25) & abs(Pek.Hin.fc) < log2(1.25))
Peking.Muscovy.hinnyduck.conserved$mode <- "conserved"

Peking.Muscovy.hinnyduck.additive <- subset(Peking.Muscovy.hinnyduck.fc, (Pek.Hin.fc > log2(1.25) & Mus.Hin.fc < -log2(1.25)) | (Pek.Hin.fc < -log2(1.25) & Mus.Hin.fc > log2(1.25)))
Peking.Muscovy.hinnyduck.additive$mode <- "additive"

Peking.Muscovy.hinnyduck.M.dominant <- subset(Peking.Muscovy.hinnyduck.fc, (abs(Pek.Hin.fc) > log2(1.25) & abs(Mus.Hin.fc) < log2(1.25)))
Peking.Muscovy.hinnyduck.M.dominant$mode <- "Muscovy.dominant"

Peking.Muscovy.hinnyduck.P.dominant <- subset(Peking.Muscovy.hinnyduck.fc, (abs(Mus.Hin.fc) > log2(1.25) & abs(Pek.Hin.fc) < log2(1.25)))
Peking.Muscovy.hinnyduck.P.dominant$mode <- "Peking.dominant"


# Muscovy.dominant: 1,115; Peking.dominant: 2,381
# overdominant: 115; underdominant: 481
# additive: 2,188; conserved: 1,859
Peking.Muscovy.hinnyduck.H.mode <- rbind(Peking.Muscovy.hinnyduck.P.dominant, Peking.Muscovy.hinnyduck.M.dominant,
                                         Peking.Muscovy.hinnyduck.overdominant, Peking.Muscovy.hinnyduck.underdominant,
                                         Peking.Muscovy.hinnyduck.additive, Peking.Muscovy.hinnyduck.conserved)

write.csv(Peking.Muscovy.hinnyduck.H.mode, "Peking.Muscovy.hinnyduck.AL.heritable.mode.csv", row.names = F, quote = F)
write.table(Peking.Muscovy.hinnyduck.H.mode, "Peking.Muscovy.hinnyduck.AL.heritable.mode.txt", row.names = F, quote = F, sep = "\t")


#################### GO enrichment
## load packages
library(clusterProfiler)
library(org.Hs.eg.db)


## For Peking.dominant genes in mule ducks
Mule.peking_dominant.gene <- subset(Peking.Muscovy.muleduck.H.mode, mode == "Peking.dominant")
Mule.peking_dominant.gene.go <- enrichGO(Mule.peking_dominant.gene$GeneName, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                                         pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
write.csv(Mule.peking_dominant.gene.go, "Muleduck.AL.Peking_dominant.gene.go.res.csv", row.names = F, quote = F)


## For Peking.dominant genes in mule ducks
Mule.muscovy_dominant.gene <- subset(Peking.Muscovy.muleduck.H.mode, mode == "Muscovy.dominant")
Mule.muscovy_dominant.gene.go <- enrichGO(Mule.muscovy_dominant.gene$GeneName, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                                          pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
write.csv(Mule.muscovy_dominant.gene.go, "Muleduck.AL.Muscovy_dominant.gene.go.res.csv", row.names = F, quote = F)


## For overdominant genes in mule ducks
Mule.overdominant.gene <- subset(Peking.Muscovy.muleduck.H.mode, mode == "overdominant")
Mule.overdominant.gene.go <- enrichGO(Mule.overdominant.gene$GeneName, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                                      pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
write.csv(Mule.overdominant.gene.go, "Muleduck.AL.overdominant.gene.go.res.csv", row.names = F, quote = F)


## For underdominant genes in mule ducks
Mule.underdominant.gene <- subset(Peking.Muscovy.muleduck.H.mode, mode == "underdominant")
Mule.underdominant.gene.go <- enrichGO(Mule.underdominant.gene$GeneName, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                                       pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
write.csv(Mule.underdominant.gene.go, "Muleduck.AL.underdominant.gene.go.res.csv", row.names = F, quote = F)


## For additive genes in mule ducks
Mule.additive.gene <- subset(Peking.Muscovy.muleduck.H.mode, mode == "additive")
Mule.additive.gene.go <- enrichGO(Mule.additive.gene$GeneName, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                                  pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
write.csv(Mule.additive.gene.go, "Muleduck.AL.additive.gene.go.res.csv", row.names = F, quote = F)


## For conserved genes in mule ducks
Mule.conserved.gene <- subset(Peking.Muscovy.muleduck.H.mode, mode == "conserved")
Mule.conserved.gene.go <- enrichGO(Mule.conserved.gene$GeneName, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                                   pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
write.csv(Mule.conserved.gene.go, "Muleduck.AL.conserved.gene.go.res.csv", row.names = F, quote = F)



## For Peking.dominant genes in Hinny ducks
Hinny.peking_dominant.gene <- subset(Peking.Muscovy.hinnyduck.H.mode, mode == "Peking.dominant")
Hinny.peking_dominant.gene.go <- enrichGO(Hinny.peking_dominant.gene$GeneName, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                                          pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
write.csv(Hinny.peking_dominant.gene.go, "Hinnyduck.AL.Peking_dominant.gene.go.res.csv", row.names = F, quote = F)


## For Peking.dominant genes in Hinny ducks
Hinny.muscovy_dominant.gene <- subset(Peking.Muscovy.hinnyduck.H.mode, mode == "Muscovy.dominant")
Hinny.muscovy_dominant.gene.go <- enrichGO(Hinny.muscovy_dominant.gene$GeneName, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                                           pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
write.csv(Hinny.muscovy_dominant.gene.go, "Hinnyduck.AL.Muscovy_dominant.gene.go.res.csv", row.names = F, quote = F)


## For overdominant genes in Hinny ducks
Hinny.overdominant.gene <- subset(Peking.Muscovy.hinnyduck.H.mode, mode == "overdominant")
Hinny.overdominant.gene.go <- enrichGO(Hinny.overdominant.gene$GeneName, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                                       pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
write.csv(Hinny.overdominant.gene.go, "Hinnyduck.AL.overdominant.gene.go.res.csv", row.names = F, quote = F)


## For underdominant genes in Hinny ducks
Hinny.underdominant.gene <- subset(Peking.Muscovy.hinnyduck.H.mode, mode == "underdominant")
Hinny.underdominant.gene.go <- enrichGO(Hinny.underdominant.gene$GeneName, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                                        pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
write.csv(Hinny.underdominant.gene.go, "Hinnyduck.AL.underdominant.gene.go.res.csv", row.names = F, quote = F)


## For additive genes in Hinny ducks
Hinny.additive.gene <- subset(Peking.Muscovy.hinnyduck.H.mode, mode == "additive")
Hinny.additive.gene.go <- enrichGO(Hinny.additive.gene$GeneName, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                                   pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
write.csv(Hinny.additive.gene.go, "Hinnyduck.AL.additive.gene.go.res.csv", row.names = F, quote = F)


## For conserved genes in Hinny ducks
Hinny.conserved.gene <- subset(Peking.Muscovy.hinnyduck.H.mode, mode == "conserved")
Hinny.conserved.gene.go <- enrichGO(Hinny.conserved.gene$GeneName, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                                    pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
write.csv(Hinny.conserved.gene.go, "Hinnyduck.AL.conserved.gene.go.res.csv", row.names = F, quote = F)
