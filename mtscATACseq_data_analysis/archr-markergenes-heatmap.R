## marker genes heatmap
print("################# Start! #########################")
as.character(R.version["version.string"])
packageVersion("ArchR")
library(ArchR)
addArchRGenome("mm10")
addArchRThreads(threads = 4)
set.seed(1234)
projBM <- loadArchRProject(path = "./Bone-Marrow-4-sample-filterDoublets")
df <- as.data.frame(getCellColData(projBM))
df$barcode <- rownames(df)
df_manual <- read.csv("Save-filterDoublets-Bone-Marrow-5/Plots/CellColData.txt", sep = "\t")
identical(df$barcode, df_manual$barcode)
projBM$Clusters2 <- df_manual$Clusters3
p1 <- plotEmbedding(projBM, colorBy = "cellColData", name = "Clusters2")
p1

markersGS <- getMarkerFeatures(
  ArchRProj = projBM, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Clusters2",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerGenes  <- c(
  "Hlf", "Mecom",	#HSC
  "Apoe", "Cd34", "Wfdc17", "Dntt", "Flt3",	#MPP
  "Pf4", "Vwf", "Gp9",	#MkP
  "Klf1", "Gata1", "Car1",	#early ERP
  # Hbb-bs, Hba-a2, Hba-a1	#ERP
  "Ms4a6c", "Cycs", "Csf1r", "Ly86",	#cMOp
  "Fcnb", "Prtn3", "Ms4a3",	#preNeu.
  # Prg3, Epx, Prg2	#Eosi. P.
  "Prss34", "Mcpt8",	#Baso. P.
  "Vpreb3", "Cd79a", "Vpreb1", "Ebf1",	#Bcell P.
  "Bank1", "Ikzf3", "Pax5",	#B cell
  "Irf4", "Mzb1", "Jchain",	#Plasma B
  "Cd3e", "Bcl11b", "Tcf7",	#Tcell P.
  # Gata3, Id2, Il7r	#ILC
  # Ncr1, Klrd1, Klrb1c	#NKP
  "Itgb7", "Naaa", "Cd74",	#cDC
  "Siglech", "Cd7", "Irf8",	#pDC
  # C1qa, C1qb	#Mac
  # Cd9, Itga2	#Mature Baso.
  # Emcn, Esam, Kdr	#Endo
  "Ngp", "S100a9", "S100a8", "Ltf"	#Mature Neu.
  
)

heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
  labelMarkers = markerGenes,
  transpose = TRUE
)
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
heatmapGS@row_order <- c(
  12, 1, 3, 4, 11, 14, 10, 7, 5, 2, 6, 9, 8, 13, 15
)
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")

plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap-4-samples", width = 9, height = 6, ArchRProj = projBM, addDOC = FALSE)
