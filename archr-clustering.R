as.character(R.version["version.string"])
packageVersion("ArchR")
library(ArchR)
addArchRGenome("mm10")
addArchRThreads(threads = 6)
set.seed(1234)


ArrowFiles <- c(
  "/results/mtscATAC/archr-results/arrow-files/Bone-Marrow-100W-ND5-G12918A-2135-66.arrow",
  "/results/mtscATAC/archr-results/arrow-files/Bone-Marrow-10W-TrnE-G14102A-1739-84.arrow",
  "/results/mtscATAC/archr-results/arrow-files/Bone-Marrow-9W-ND5-G12918A-3683-62.arrow",
  "/results/mtscATAC/archr-results/arrow-files/Bone-Marrow-100W-TrnE-G14102A-4303-79.arrow"
)

projBM <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "Bone-Marrow",
  copyArrows = TRUE
)
projBM <- filterDoublets(projBM)
# projBM <- saveArchRProject(ArchRProj = projBM, outputDirectory = "Bone-Marrow", load = TRUE)
## change iterations to 4
projBM <- addIterativeLSI(
  ArchRProj = projBM,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 4, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(2), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:30
)

projBM <- addClusters(
  input = projBM,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 1.2,
  force = TRUE
)

projBM <- addUMAP(
  ArchRProj = projBM, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 40, 
  minDist = 0.4, 
  metric = "cosine",
  force = TRUE
)

p1 <- plotEmbedding(ArchRProj = projBM, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = projBM, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
plotPDF(p1,p2, name = "Plot-UMAP-Sample-Clusters-filterDoublets-IterativeLSI.pdf", ArchRProj = projBM, addDOC = FALSE, width = 5, height = 5)
projBM <- saveArchRProject(ArchRProj = projBM, outputDirectory = "Save-filterDoublets-Bone-Marrow-1", load = TRUE)

