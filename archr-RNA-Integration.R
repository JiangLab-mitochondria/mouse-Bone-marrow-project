print("################# Start! #########################")
as.character(R.version["version.string"])
packageVersion("ArchR")
library(ArchR)
addArchRGenome("mm10")
addArchRThreads(threads = 1)
library(Seurat)
set.seed(1234)
projBM <- loadArchRProject(path = "./Bone-Marrow-4-sample-filterDoublets")
seRNA <- readRDS("/results/scRNA/seurat-results/clustering-sample-merged/three_samples_merged_harmony_dims15_celltype.rds")
seRNA
table(seRNA$cell_type)

# Unconstrained Integration
projBM <- addGeneIntegrationMatrix(
  ArchRProj = projBM, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = seRNA,
  addToArrow = FALSE,
  groupRNA = "cell_type",
  nameCell = "predictedCell_Un",
  nameGroup = "predictedGroup_Un",
  nameScore = "predictedScore_Un"
)
df <- as.data.frame(getCellColData(projBM))
df$barcode <- rownames(df)
write.table(df,
            file="Bone-Marrow-4-sample-filterDoublets/Plots/CellColData.txt",
            sep="\t", quote=F, row.names = F, col.names = T)


df_rna <- read.csv("Bone-Marrow-4-sample-filterDoublets/Plots/CellColData.txt", sep = "\t")
df_manual <- read.csv("Save-filterDoublets-Bone-Marrow-5/Plots/CellColData.txt", sep = "\t")
identical(df_rna$barcode, df_manual$barcode)
df_manual$predictedGroup_Un <- df_rna$predictedGroup_Un
table(df_manual$Clusters3)
df_manual <- subset(df_manual, !(Clusters3 %in% c("Unassigned")))
## confusion matrix
cm <- table(df_manual$Clusters3, df_manual$predictedGroup_Un)
prop.table(cm, 1) * 100
library(pheatmap)
p1 <- pheatmap(prop.table(cm, 1) * 100,
               cluster_rows = FALSE,
               border_color = "black",
               cluster_cols = FALSE
)

p1
pdf(file = "Save-filterDoublets-Bone-Marrow-5/Plots/GeneIntegration_confusion_matrix.pdf", 
    width = 7, height = 6)
p1
dev.off()