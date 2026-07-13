print("################# Start! #########################")
as.character(R.version["version.string"])
packageVersion("ArchR")
library(ArchR)
addArchRGenome("mm10")
addArchRThreads(threads = 4)
set.seed(1234)
setwd("/results/mtscATAC/archr-results/")
projBM <- loadArchRProject(path = "./Save-filterDoublets-Bone-Marrow-1")

## Custom Function
cell_type_define <- function(IterativeLSI="IterativeLSI") {
  projBM <- addClusters(
    input = projBM,
    reducedDims = IterativeLSI,
    method = "Seurat",
    name = "Clusters",
    resolution = 1.2,
    force = TRUE
  )
  
  projBM <- addUMAP(
    ArchRProj = projBM, 
    reducedDims = IterativeLSI, 
    name = "UMAP", 
    nNeighbors = 40, 
    minDist = 0.4, 
    metric = "cosine",
    force = TRUE
  )
  p1 <- plotEmbedding(ArchRProj = projBM, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
  p2 <- plotEmbedding(ArchRProj = projBM, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
  features <- list(
    HSCScore = c("Hlf", "Mecom")
  )
  projBM <- addModuleScore(projBM,
                           useMatrix = "GeneScoreMatrix",
                           name = "Module",
                           features = features)
  features <- list(
    MPPScore = c("Apoe", "Cd34",
                 "Wfdc17", "Dntt", "Flt3")
  )
  projBM <- addModuleScore(projBM,
                           useMatrix = "GeneScoreMatrix",
                           name = "Module",
                           features = features)
  
  features <- list(
    MkPScore = c(
      "Pf4", "Vwf", "Gp9"
    )
  )
  projBM <- addModuleScore(projBM,
                           useMatrix = "GeneScoreMatrix",
                           name = "Module",
                           features = features)
  
  features <- list(
    earlyERPScore = c(
      "Klf1", "Gata1", "Car1"
    )
  )
  projBM <- addModuleScore(projBM,
                           useMatrix = "GeneScoreMatrix",
                           name = "Module",
                           features = features)
  
  features <- list(
    cMOpScore = c("Ms4a6c", "Cycs", "Csf1r",  "Ly86")
  )
  projBM <- addModuleScore(projBM,
                           useMatrix = "GeneScoreMatrix",
                           name = "Module",
                           features = features)
  
  features <- list(
    preNeuScore = c(
      "Fcnb", "Prtn3", "Ms4a3"
    )
  )
  projBM <- addModuleScore(projBM,
                           useMatrix = "GeneScoreMatrix",
                           name = "Module",
                           features = features)
  
  features <- list(
    MatureNeutrophilScore = c("Ngp", "S100a9", "S100a8", "Ltf")
  )
  projBM <- addModuleScore(projBM,
                           useMatrix = "GeneScoreMatrix",
                           name = "Module",
                           features = features)
  
  features <- list(
    BasoPScore = c(
      "Prss34", "Mcpt8"
    )
  )
  projBM <- addModuleScore(projBM,
                           useMatrix = "GeneScoreMatrix",
                           name = "Module",
                           features = features)
  
  features <- list(
    BcellPScore = c(
      "Vpreb3", "Cd79a", "Vpreb1", "Ebf1"
    )
  )
  projBM <- addModuleScore(projBM,
                           useMatrix = "GeneScoreMatrix",
                           name = "Module",
                           features = features)
  
  features <- list(
    BcellScore = c(
      "Bank1", "Ikzf3", "Pax5"
    )
  )
  projBM <- addModuleScore(projBM,
                           useMatrix = "GeneScoreMatrix",
                           name = "Module",
                           features = features)
  
  features <- list(
    PlasmaBScore = c("Irf4", "Mzb1", "Jchain")
  )
  projBM <- addModuleScore(projBM,
                           useMatrix = "GeneScoreMatrix",
                           name = "Module",
                           features = features)
  
  features <- list(
    TcellPScore = c(
      "Cd3e", "Bcl11b", "Tcf7"
    )
  )
  projBM <- addModuleScore(projBM,
                           useMatrix = "GeneScoreMatrix",
                           name = "Module",
                           features = features)
  features <- list(
    cDCScore = c(
      "Itgb7", "Naaa", "Cd74"
    )
  )
  projBM <- addModuleScore(projBM,
                           useMatrix = "GeneScoreMatrix",
                           name = "Module",
                           features = features)
  
  features <- list(
    pDCScore = c(
      "Siglech", "Cd7", "Irf8"
    )
  )
  projBM <- addModuleScore(projBM,
                           useMatrix = "GeneScoreMatrix",
                           name = "Module",
                           features = features)
  p3 <- plotEmbedding(projBM,
                      embedding = "UMAP",
                      colorBy = "cellColData",
                      name="Module.HSCScore",
                      imputeWeights = getImputeWeights(projBM))
  
  p4 <- plotEmbedding(projBM,
                      embedding = "UMAP",
                      colorBy = "cellColData",
                      name="Module.MPPScore",
                      imputeWeights = getImputeWeights(projBM))
  p5 <- plotEmbedding(projBM,
                      embedding = "UMAP",
                      colorBy = "cellColData",
                      name="Module.MkPScore",
                      imputeWeights = getImputeWeights(projBM))
  p6 <- plotEmbedding(projBM,
                      embedding = "UMAP",
                      colorBy = "cellColData",
                      name="Module.earlyERPScore",
                      imputeWeights = getImputeWeights(projBM))
  p7 <- plotEmbedding(projBM,
                      embedding = "UMAP",
                      colorBy = "cellColData",
                      name="Module.cMOpScore",
                      imputeWeights = getImputeWeights(projBM))
  p8 <- plotEmbedding(projBM,
                      embedding = "UMAP",
                      colorBy = "cellColData",
                      name="Module.preNeuScore",
                      imputeWeights = getImputeWeights(projBM))
  p9 <- plotEmbedding(projBM,
                      embedding = "UMAP",
                      colorBy = "cellColData",
                      name="Module.MatureNeutrophilScore",
                      imputeWeights = getImputeWeights(projBM))
  p10 <- plotEmbedding(projBM,
                       embedding = "UMAP",
                       colorBy = "cellColData",
                       name="Module.BasoPScore",
                       imputeWeights = getImputeWeights(projBM))
  p11 <- plotEmbedding(projBM,
                       embedding = "UMAP",
                       colorBy = "cellColData",
                       name="Module.BcellPScore",
                       imputeWeights = getImputeWeights(projBM))
  p12 <- plotEmbedding(projBM,
                       embedding = "UMAP",
                       colorBy = "cellColData",
                       name="Module.BcellScore",
                       imputeWeights = getImputeWeights(projBM))
  p13 <- plotEmbedding(projBM,
                       embedding = "UMAP",
                       colorBy = "cellColData",
                       name="Module.PlasmaBScore",
                       imputeWeights = getImputeWeights(projBM))
  p14 <- plotEmbedding(projBM,
                       embedding = "UMAP",
                       colorBy = "cellColData",
                       name="Module.TcellPScore",
                       imputeWeights = getImputeWeights(projBM))
  p15 <- plotEmbedding(projBM,
                       embedding = "UMAP",
                       colorBy = "cellColData",
                       name="Module.cDCScore",
                       imputeWeights = getImputeWeights(projBM))
  p16 <- plotEmbedding(projBM,
                       embedding = "UMAP",
                       colorBy = "cellColData",
                       name="Module.pDCScore",
                       imputeWeights = getImputeWeights(projBM))
  
  plotPDF(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16, 
          name = "Plot-UMAP-Sample-Clusters-filterDoublets-IterativeLSI-celltype.pdf", ArchRProj = projBM, addDOC = FALSE, width = 5, height = 5)
  return() 
}


results <- cell_type_define(IterativeLSI="IterativeLSI")





