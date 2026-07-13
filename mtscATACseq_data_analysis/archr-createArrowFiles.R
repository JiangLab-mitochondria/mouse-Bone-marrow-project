as.character(R.version["version.string"])
packageVersion("ArchR")
library(ArchR)
addArchRGenome("mm10")
addArchRThreads(threads = 4)
set.seed(1234)

exclude.Chr <- c(
  "chrM", "chrY",
  "JH584299.1","GL456233.1","GL456211.1","GL456350.1","JH584293.1","GL456221.1",
  "JH584297.1","JH584296.1","GL456354.1","JH584294.1","JH584298.1","GL456219.1",
  "GL456210.1","JH584303.1","GL456212.1","JH584304.1","GL456216.1","JH584292.1",
  "JH584295.1"
)
archr_createArrowFiles <- function(sample) {
  sample <- as.character(sample)
  print(sample)
  print("==================start=======================")
  inputfile <- paste0("/results/mtscATAC/cellranger-atac-results/",
                      sample, "/outs/fragments.tsv.gz")
  sample.name <- sample
  ArrowFiles <- createArrowFiles(
    inputFiles = inputfile,
    sampleNames = sample.name,
    force = T,
    excludeChr = exclude.Chr
  )
  
  print(sample)
  print("==================end=======================")
  print("")
  return()
}



samples <- c(
  "Bone-Marrow-100W-ND5-G12918A-2135-66",
  "Bone-Marrow-10W-TrnE-G14102A-1739-84",
  "Bone-Marrow-9W-ND5-G12918A-3683-62",
  "Bone-Marrow-100W-TrnE-G14102A-4303-79"
)
for (i in samples) {
  result <- archr_createArrowFiles(sample = i)
}
print("END")
