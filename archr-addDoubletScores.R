as.character(R.version["version.string"])
packageVersion("ArchR")
library(ArchR)
addArchRGenome("mm10")
addArchRThreads(threads = 4)
set.seed(1234)
setwd("/results/mtscATAC/archr-results/arrow-files/")

ArrowFiles <- c(
  "Bone-Marrow-100W-ND5-G12918A-2135-66.arrow",
  "Bone-Marrow-10W-TrnE-G14102A-1739-84.arrow",
  "Bone-Marrow-9W-ND5-G12918A-3683-62.arrow",
  "Bone-Marrow-100W-TrnE-G14102A-4303-79.arrow"
)
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
  LSIMethod = 1
)