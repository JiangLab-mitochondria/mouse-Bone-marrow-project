library(ggplot2)
set.seed(1234)

df <- read.csv("results/mtscATAC/mgatk-results/four_mut_samples_mutation_load_cell_type_final.txt", sep = "\t")
colnames(df)
## filter cells: Total_coverage >= 10
df1 <- subset(df, Total_coverage >= 10)
table(df1$Sample)
df1$Sample <- factor(df1$Sample, levels = c(
  "Bone-Marrow-9W-ND5-G12918A-3683-62", "Bone-Marrow-100W-ND5-G12918A-2135-66",
  "Bone-Marrow-10W-TrnE-G14102A-1739-84", "Bone-Marrow-100W-TrnE-G14102A-4303-79"
))
## UMAP 
p1 <- ggplot(df1, aes(x=UMAP_1, y=UMAP_2)) +
  geom_point(aes(color = Mutation_load), position = "jitter", size=0.3) +
  scale_color_gradient(high = "yellow", low = "darkblue")  +
  theme_classic() +
  facet_wrap(~Sample, ncol = 2) +
  xlab("UMAP_1") +
  ylab("UMAP_2")
p1
pdf(file = "pdf-files/Umap-mutation-load-cov10-reclustering.pdf", width = 8, height = 6)
p1
p1 + theme(legend.position="none")

dev.off()


## Cumulative Distribution
library(ArchR)
addArchRThreads(threads = 3)
library(dplyr)
set.seed(1234)

## remove unassigned cells
table(df1$Clusters3)
pal <- paletteDiscrete(values = df1$Clusters3)
pal
df1 <- subset(df1, !(Clusters3 %in% c("Unassigned")))
## filter cell types wiht cell count < 10
table(df1$Sample, df1$Clusters3)
df_filtered <- df1 %>%
  group_by(Sample, Clusters3) %>%                
  filter(n() >= 10) %>%                          
  ungroup()


p1 <- ggplot(df_filtered, aes(x=Mutation_load, col=Clusters3)) +
  stat_ecdf(linewidth=1, pad = F) +
  theme_classic() +
  xlab("Mutation load") +
  ylab("Cumulative Probability") +
  scale_color_manual(values=pal) +
  facet_wrap(~Sample, ncol = 2) +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    # axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
    legend.title = element_blank()
  ) 
# labs(title = "single cell mutation load") 

p1  

pdf(file = "pdf-files/ecdf-mutation-load-cov10-cell10-reclustering.pdf", width = 8, height = 6)
p1
p1 + theme(legend.position="none")
dev.off()


## Histogram
df_filtered <- df1 %>%
  group_by(Sample, Clusters3) %>%                
  filter(n() >= 10) %>%                          
  ungroup()
pal2 <- c(
  "Bone-Marrow-9W-ND5-G12918A-3683-62"="#D51F26",
  "Bone-Marrow-100W-ND5-G12918A-2135-66"="#272E6A"
)
df3 <- as.data.frame(df_filtered)
table(df3$Sample)
df4 <- subset(df3, Sample %in% c(
  "Bone-Marrow-9W-ND5-G12918A-3683-62",
  "Bone-Marrow-100W-ND5-G12918A-2135-66"
))
p1 <- ggplot(df4, aes(fill = Sample, x = Mutation_load, color = Sample)) +
  geom_histogram(alpha=0.3, binwidth  = 3, aes(y=after_stat(density)), position = "identity") +
  geom_density(alpha=.2) +
  theme_classic() +
  xlab("Mutation load") +
  scale_fill_manual(values=pal2) +
  scale_color_manual(values=pal2) +
  facet_wrap(~Clusters3, ncol = 5, scales = "free_y") +
  theme(plot.title = element_text(hjust=0.5)) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "bottom")
p1


pal2 <- c(
  "Bone-Marrow-10W-TrnE-G14102A-1739-84"="#D51F26",
  "Bone-Marrow-100W-TrnE-G14102A-4303-79"="#272E6A"
)
df3 <- as.data.frame(df_filtered)
table(df3$Sample)
df4 <- subset(df3, Sample %in% c(
  "Bone-Marrow-10W-TrnE-G14102A-1739-84",
  "Bone-Marrow-100W-TrnE-G14102A-4303-79"
))
p2 <- ggplot(df4, aes(fill = Sample, x = Mutation_load, color = Sample)) +
  geom_histogram(alpha=0.3, binwidth  = 3, aes(y=after_stat(density)), position = "identity") +
  geom_density(alpha=.2) +
  theme_classic() +
  xlab("Mutation load") +
  scale_fill_manual(values=pal2) +
  scale_color_manual(values=pal2) +
  facet_wrap(~Clusters3, ncol = 5, scales = "free_y") +
  theme(plot.title = element_text(hjust=0.5)) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = "bottom")
p2
pdf(file = "pdf-files/histogram-mutation-load-cov10-cell10-reclustering.pdf", width = 8, height = 6)
p1
p1 + theme(legend.position="none")
p2
p2 + theme(legend.position="none")

dev.off()
