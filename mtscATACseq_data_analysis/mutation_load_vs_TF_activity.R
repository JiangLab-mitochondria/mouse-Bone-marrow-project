library(dplyr)
library(ggplot2)
set.seed(1234)
setwd("/storage2/jiangminLab/chaiguoshi/projects/lv-qun-yu-projects/results/mtscATAC/archr-results/")
sessionInfo()
# R version 4.2.1 (2022-06-23)

print("HSC start")
df <- read.csv("SubsetHSC_mutation_load_type_path_revised_MotifMatrix.txt", sep = "\t")
table(df$Sample)
sample_names <- unique(df$Sample)
n_samples <- length(sample_names)

for (sample_id in 1:n_samples){
  sample_name <- sample_names[sample_id]
  print(paste0(sample_name, " start"))
  df1 <- subset(df, Sample %in% c(sample_name))
  df2 <- subset(df1, Total_coverage >= 5)
  TFs <- colnames(df2)
  TFs <- TFs[1:884]
  n_tfs <- length(TFs)
  
  final_tf_names <- character(n_tfs)
  final_z_scores   <- numeric(n_tfs)
  final_p_values   <- numeric(n_tfs)
  
  for (i in 1:n_tfs){
    current_tf_name <- TFs[i]
    TF_zscore <- df2[[current_tf_name]]
    Mutation_load <- df2[["Mutation_load"]]
    r_obs <- cor(TF_zscore, Mutation_load, method = "spearman")
    
    # ==========================================
    #  Permutation Test
    # ==========================================
    num_permutations <- 10000 # Number of permutations (recommended >= 1000; larger values yield higher precision but increase computation time)
    r_null <- numeric(num_permutations) # Used to store the correlation coefficients after each permutation
    
    
    message("Permutation testing is in progress...")
    
    
    for (j in 1:num_permutations) {
      # This simulates the null hypothesis (H0): there is no association between the two.
      shuffled_het <- sample(Mutation_load) 
      
      # Compute the correlation coefficient after shuffling
      r_null[j] <- cor(TF_zscore, shuffled_het, method = "spearman")
    }
    
    
    # ==========================================
    #  Calculation of Z-scores and P-values
    # ==========================================
    
    
    # Compute the mean and standard deviation of the null distribution
    mean_null <- mean(r_null)
    sd_null <- sd(r_null)
    
    
    # Calculate the Z-score: (observed value - mean of the null distribution) / standard deviation of the null distribution
    z_score <- (r_obs - mean_null) / sd_null
    
    
    # Calculate the empirical P-value (two-tailed test)
    # Method: Compute the proportion of permutation results with absolute values greater than or equal to |r_obs|
    # Adding 1 prevents the P-value from being zero, serving as a conservative estimation approach
    p_value <- (sum(abs(r_null) >= abs(r_obs)) + 1) / (num_permutations + 1)
    
    final_tf_names[i] <- current_tf_name
    final_z_scores[i] <- z_score
    final_p_values[i] <- p_value
  }
  results_df <- data.frame(
    TF_Name = final_tf_names,
    Z_Score = final_z_scores,
    P_Value = final_p_values,
    stringsAsFactors = FALSE
  )
  write.table(results_df,
              file=paste0("response-to-reviewer1/SubsetHSC_mutation_load_type_", sample_name, "-TFactivity-Zcore-Pvalue.txt"),
              sep="\t", quote=F, row.names = F, col.names = T)
  
  print(paste0(sample_name, " end"))
  
}
print("HSC end")






## added a sample
library(dplyr)
library(ggplot2)
set.seed(1234)
setwd("/storage2/jiangminLab/chaiguoshi/projects/lv-qun-yu-projects/results/mtscATAC/archr-results/")
sessionInfo()
# R version 4.2.1 (2022-06-23)

print("HSC start")
df <- read.csv("response-to-reviewer1/Bone-Marrow-100W-ND5-3673-66-Female-SubsetHSC_mutation_load_type_MotifMatrix.txt", sep = "\t")
table(df$Sample)
sample_names <- unique(df$Sample)
n_samples <- length(sample_names)

for (sample_id in 1:n_samples){
  sample_name <- sample_names[sample_id]
  print(paste0(sample_name, " start"))
  df1 <- subset(df, Sample %in% c(sample_name))
  df2 <- subset(df1, Total_coverage >= 5)
  TFs <- colnames(df2)
  TFs <- TFs[1:884]
  n_tfs <- length(TFs)
  
  final_tf_names <- character(n_tfs)
  final_z_scores   <- numeric(n_tfs)
  final_p_values   <- numeric(n_tfs)
  
  for (i in 1:n_tfs){
    current_tf_name <- TFs[i]
    TF_zscore <- df2[[current_tf_name]]
    Mutation_load <- df2[["Mutation_load"]]
    r_obs <- cor(TF_zscore, Mutation_load, method = "spearman")
    
    # ==========================================
    #  Permutation Test
    # ==========================================
    num_permutations <- 10000 
    r_null <- numeric(num_permutations)
    
    
    message("Permutation testing is in progress...")
    
    
    for (j in 1:num_permutations) {
      
      shuffled_het <- sample(Mutation_load) 
      
      
      r_null[j] <- cor(TF_zscore, shuffled_het, method = "spearman")
    }
    
    
    # ==========================================
    #  Calculation of Z-scores and P-values
    # ==========================================
    
    
    
    mean_null <- mean(r_null)
    sd_null <- sd(r_null)
    
    
    # 
    z_score <- (r_obs - mean_null) / sd_null
    
    
    # 
    #
    p_value <- (sum(abs(r_null) >= abs(r_obs)) + 1) / (num_permutations + 1)
    
    final_tf_names[i] <- current_tf_name
    final_z_scores[i] <- z_score
    final_p_values[i] <- p_value
  }
  results_df <- data.frame(
    TF_Name = final_tf_names,
    Z_Score = final_z_scores,
    P_Value = final_p_values,
    stringsAsFactors = FALSE
  )
  write.table(results_df,
              file=paste0("response-to-reviewer1/SubsetHSC_mutation_load_type_", sample_name, "-TFactivity-Zcore-Pvalue.txt"),
              sep="\t", quote=F, row.names = F, col.names = T)
  
  print(paste0(sample_name, " end"))
  
}
print("HSC end")


df <- read.csv("response-to-reviewer1/SubsetHSC_mutation_load_type_Bone-Marrow-100W-ND5-G12918A-3673-66-Female-TFactivity-Zcore-Pvalue.txt", sep = "\t")
df <- df %>%
  mutate(Type = case_when(
    
    Z_Score >= 1.96 & P_Value <= 0.05 ~ "PositiveCorrelation",
   
    Z_Score <= -1.96 & P_Value <= 0.05 ~ "NegativeCorrelation",
   
    TRUE ~ "NoCorrelation"
  ))

table(df$Type)


highlight_genes <- c("Eomes_768", "Arid3a_7", "Tcf3_31", "Klf10_810", "Tcf12_59", "Tcf4_88",
                     "Fos_104", "Fosl1_107", "Fosb_98", "Fosl2_113", 
                     "Jund_135", "Jun_126", "Junb_127",
                     
                     "Klf10_810", "Tcf3_31", "Nfil3_131", "Zbtb7a_808", "Id3_22", "Tbx21_762",
                     "Fos_104", "Fosl1_107", "Fosb_98", "Fosl2_113", 
                     "Jund_135", "Jun_126", "Junb_127"
)
show_genes <- subset(df, TF_Name %in% highlight_genes)


p1 <- ggplot(df) +
  geom_point(data = subset(df, Type == "NoCorrelation"), aes(x = Z_Score, y = -log10(P_Value)), size = 1, color="grey") +
  geom_point(data = subset(df, Type == "PositiveCorrelation"), aes(x = Z_Score, y = -log10(P_Value)), size = 2, alpha=0.8, color = "#D51F26") +
  geom_point(data = subset(df, Type == "NegativeCorrelation"), aes(x = Z_Score, y = -log10(P_Value)), size = 2, alpha=0.8, color = "#8A9FD1") +
  labs(title="Aged ND5 HSC Female", x="Z Score", y = "-log10(P Value)") +
  theme_classic()  +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") + 
  geom_vline(xintercept = c(1.96, -1.96),
             linetype = "dashed") +
  theme(plot.title = element_text(hjust=0.5),
        legend.position = "top") +
  geom_text_repel(data = show_genes,
                  aes(x=Z_Score, y=-log10(P_Value), label = TF_Name),hjust= 0.30,
                  max.overlaps=100, size=5,
                  min.segment.length = 0, colour="black",
  ) 
p1
## save df
write.table(df,
            file="response-to-reviewer1/SubsetHSC_mutation_load_type_Bone-Marrow-100W-ND5-G12918A-3673-66-Female-TFactivity-Zcore-Pvalue-added.txt",
            sep="\t", quote=F, row.names = F, col.names = T)
pdf(file = "response-to-reviewer1/TF-activity-vs-mutation-load-HSC-Female-highlight-cov5-Zscore-Pvalue.pdf",
    width =5, height = 5)
p1
dev.off()

