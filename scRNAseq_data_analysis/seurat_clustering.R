library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
set.seed(1234)

# ND5 sample
df1 <- Read10X("resutls/scRNA/cellranger-results/ND5/outs/filtered_feature_bc_matrix")
df1 <- CreateSeuratObject(counts = df1, min.cells = 3, min.features = 200, project = "ND5")
df1
head(df1)
df1[["percent.mt"]] <- PercentageFeatureSet(df1, pattern = "^mt-")
df1[["percent.ribo"]] <- PercentageFeatureSet(df1, pattern = "^Rp[sl]")
saveRDS(df1, file = "resutls/scRNA/seurat-results/clustering-single-sample/ND5.rds")

# TrnE sample
df1 <- Read10X("resutls/scRNA/cellranger-results/TrnE/outs/filtered_feature_bc_matrix")
df1 <- CreateSeuratObject(counts = df1, min.cells = 3, min.features = 200, project = "TrnE")
df1
head(df1)
df1[["percent.mt"]] <- PercentageFeatureSet(df1, pattern = "^mt-")
df1[["percent.ribo"]] <- PercentageFeatureSet(df1, pattern = "^Rp[sl]")
saveRDS(df1, file = "resutls/scRNA/seurat-results/clustering-single-sample/TrnE.rds")

# WT sample
df1 <- Read10X("resutls/scRNA/cellranger-results/wt/outs/filtered_feature_bc_matrix")
df1 <- CreateSeuratObject(counts = df1, min.cells = 3, min.features = 200, project = "WT")
df1
head(df1)
df1[["percent.mt"]] <- PercentageFeatureSet(df1, pattern = "^mt-")
df1[["percent.ribo"]] <- PercentageFeatureSet(df1, pattern = "^Rp[sl]")
saveRDS(df1, file = "resutls/scRNA/seurat-results/clustering-single-sample/WT.rds")

ND5 <- readRDS("resutls/scRNA/seurat-results/clustering-single-sample/ND5.rds")
TrnE <- readRDS("resutls/scRNA/seurat-results/clustering-single-sample/TrnE.rds")
WT <- readRDS("resutls/scRNA/seurat-results/clustering-single-sample/WT.rds")


ND5 <- subset(ND5, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & nCount_RNA < 40000 & percent.mt < 5)
TrnE <- subset(TrnE, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & nCount_RNA < 40000 & percent.mt < 5)
WT <- subset(WT, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & nCount_RNA < 40000 & percent.mt < 5)

all.merged <- merge(ND5, c(TrnE, WT),
                    add.cell.ids = c("ND5",
                                     "TrnE",
                                     "WT"
                    ))
all.merged <- NormalizeData(all.merged)
g2m.genes <- str_to_title(tolower(cc.genes$g2m.genes))
s.genes <- str_to_title(tolower(cc.genes$s.genes))
head(all.merged)
all.merged <- CellCycleScoring(all.merged, s.features = s.genes, g2m.features = g2m.genes)
head(all.merged)
all.merged <- FindVariableFeatures(all.merged, selection.method = "vst", nfeatures = 2000)
all.merged <- ScaleData(all.merged)
all.merged <- RunPCA(all.merged)
all.merged.harmony <- RunHarmony(all.merged, "orig.ident")
all.merged.harmony

all.merged.harmony <- RunUMAP(all.merged.harmony, reduction = "harmony", dims = 1:15)
all.merged.harmony <- RunTSNE(all.merged.harmony, reduction = "harmony", dims = 1:15)
all.merged.harmony <- FindNeighbors(all.merged.harmony, reduction = "harmony", dims = 1:15)
all.merged.harmony <- FindClusters(all.merged.harmony, resolution = 0.8)
all.merged.harmony
all.merged.harmony <- FindClusters(all.merged.harmony, resolution = 1.2)
all.merged.harmony
saveRDS(all.merged.harmony, file = "resutls/scRNA/seurat-results/clustering-sample-merged/three_samples_merged_harmony_dims15.rds")





