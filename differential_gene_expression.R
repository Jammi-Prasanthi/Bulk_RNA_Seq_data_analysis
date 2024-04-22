Once the count matrix is generated we can proceed with differential gene expression, pathway analysis, etc using R studio
```
#set the directories
setwd("./")
getwd()

#get required packages
library(DESeq2)
library(RColorBrewer)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(data.table)
library(ggpubr)
library(tidyr)
library(gprofiler2)
library(pheatmap)
library(svglite)
source("https://raw.githubusercontent.com/ankitasks1/utilities/main/upset_out.R") 
library(UpSetR)
library(DGEobj.utils) 
library(tidyverse)
library(gridExtra)
library(readxl)
library(BiocManager)
```

#DATA LOADING AND DATA CLEANING
#The gene ids of the species that I have been working are not recognized by many tools for GO and pathway analysis, hence I have added corresponding gene ids of closely related species and I am retriving this information in the below code stored under annotation variable.

```
count_data <- read_excel("master_count_matrix.xlsx", sheet = 1)
count_data <- as.data.frame(count_data)
dim(count_data)
colnames(count_data)

# retaining gene id, corresponding id for later analysis
annotation <- count_data[, c("gene_id", "locus_ID", "gene_descript")]
rownames(annotation) <- annotation$gene_id
dim(annotation)

#preparing the count matrix in suitable form
count_data_df <- count_data
count_data_df <- count_data_df[,-c(14:30)] # removed any extra columns that are not counts
rownames(count_data_df) <- count_data_df$gene_id #make geneids as rows of the data frame

#remove the column containing the Geneid values, which have been set as row names
count_data_df <- count_data_df[,-1]
dim(count_data_df)
colnames(count_data_df)

#create a matrix of the data frame
count_data_mat <- as.matrix(count_data_df)

```
#PREPARING/LOADING THE METADATA
If a metadata file is existing, load it else create one

```
#you can create a metadata here or load a pre-existing one using fread()
coldata <- read.csv("coldata.csv")
coldata <- data.frame(coldata) #create a data frame
rownames(coldata) <- coldata$samplename #make the sample names as rows

#converting the column values of conditions, replicates and type into factor datatype
coldata$conditions <- factor(coldata$conditions)
coldata$replicates <- factor(coldata$replicates)
coldata$type <- factor(coldata$disease)
coldata$samplename <- factor(coldata$samplename)
coldata$genotype <- factor(coldata$genotype)
coldata$tolerance_level <- factor(coldata$tolerance_level)
coldata

# Reorder the columns in count_data_mat to match the row order in coldata
count_data_mat <- count_data_mat[, rownames(coldata)]

#check whether the rownames of coldata and column names of count matrix are same using all() function 
all(rownames(coldata) == colnames(count_data_mat)) #should print TRUE
```
#NORMALIZE THE DATA FOR THE SEQUENCE LENGTH AND LIBRARY SIZE
TPM calculation requires count matrix and lengths of genes (avaible in the count matrix geberated with featureCounts)
```
# use tpm.direct() function of DGEobj.utils package
plant_tpm <- tpm.direct(countsMatrix = count_data_mat, geneLength = count_data$gene_length[match(rownames(count_data_mat), count_data$gene_id)], collapse = FALSE)
plant_tpm_df <- as.data.frame(plant_tpm) # Convert the TPM matrix to a dataframe
```
#DIMENSIONALITY REDUCTION: PCA

```
log_tpm_df = log2(plant_tpm_df + 1) #log transformation of TPM 
pcs = prcomp(t(log_tpm_df), center = TRUE) # Transpose the dataframe & Perform PCA
percentVar = round(((pcs$sdev) ^ 2 / sum((pcs$sdev) ^ 2)* 100), 2) # Calculate percent variance
pca_df <- data.frame(pcs$x) # Create a data frame from the PCA results
sample_names <- rownames(pcs$x) # Get the names of the samples from the PCA results

# Match these sample names to the sample names in coldata
matched_conditions <- coldata$conditions[match(sample_names, rownames(coldata))]
pca_df$Condition <- matched_conditions # Add these matched conditions to the PCA results

#PCA plot of TPM values
plant_TPM_PCAplot <- ggplot(pca_df, aes(x = PC1, y = PC2, colour = Condition, shape = Condition, fill = Condition)) + 
                          xlab(paste0("PC1: ", percentVar[1], "% variance")) +
                          ylab(paste0("PC2: ", percentVar[2], "% variance")) +
                          geom_point(size = 5, alpha = 1) +
                          scale_fill_manual(values = c("#00FFFF","#0000FF", "pink", "red"))+
                          scale_color_manual(values = c("black","black","black","black"))+
                          scale_shape_manual(values=c(21,21, 24,24)) +
                          theme(legend.text = element_text(size = 16, face = "bold"),
                                legend.title = element_text(size = 16, colour = "black", face = "bold"),
                                axis.title = element_text(size = 18, face = "bold"),
                                axis.text.x = element_text(size = 16, face = "bold", color = "black"),
                                axis.text.y = element_text(size = 16, face = "bold", color = "black"),
                                plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))
plant_TPM_PCAplot

#PCA plot saved in below pdf
ggsave("plant_log2TPM_12hr_alternaria_PCAplot.svg", width=17*1.25, height=12*1.25, units="cm") 
```
#CALCULATE THE AVERAGES OF TPM

```
# Create a copy of the original DataFrame
plant_TPM_avg <- log_tpm_df

# Replace specific identifiers with generic ones
plant_TPM_avg$Group1_Var1_avg <- rowMeans(log_tpm_df[, grep("Group1_Var1_", colnames(log_tpm_df))])
plant_TPM_avg$Group1_Var2_avg <- rowMeans(log_tpm_df[, grep("Group1_Var2_", colnames(log_tpm_df))])
plant_TPM_avg$Group2_Var1_avg <- rowMeans(log_tpm_df[, grep("Group2_Var1_", colnames(log_tpm_df))])
plant_TPM_avg$Group2_Var2_avg <- rowMeans(log_tpm_df[, grep("Group2_Var2_", colnames(log_tpm_df))])
plant_TPM_avg$Treatment1_Var1_avg <- rowMeans(log_tpm_df[, grep("Treatment1_Var1_", colnames(log_tpm_df))])
plant_TPM_avg$Treatment1_Var2_avg <- rowMeans(log_tpm_df[, grep("Treatment1_Var2_", colnames(log_tpm_df))])
plant_TPM_avg$Treatment2_Var1_avg <- rowMeans(log_tpm_df[, grep("Treatment2_Var1_", colnames(log_tpm_df))])
plant_TPM_avg$Treatment2_Var2_avg <- rowMeans(log_tpm_df[, grep("Treatment2_Var2_", colnames(log_tpm_df))])

# Removing the initial replicate columns
columns_to_remove <- grep("_\\d$", colnames(plant_TPM_avg), value = TRUE)
plant_TPM_avg <- plant_TPM_avg[ , !(colnames(plant_TPM_avg) %in% columns_to_remove)]

```

#DIFFERENTIAL GENE EXPRESSION USING DESeq2
#Its better to write a loop to iterate over all in case of multiple constrasts 
```
dds <- DESeqDataSetFromMatrix(countData = count_data_mat, colData = coldata, design = ~ conditions) 
keep <- rowSums(counts(dds)) >= 10 #filter out low count genes: keep only those rows whose sum is >=10
dds <- dds[keep,]
dds

my_DE_list <- list()
my_DE_list_0.05 <- list()
my_DE_list_0.05_up <- list()
my_DE_list_0.05_down <- list()
plot_MA_list <- list()
for (base1 in unique(as.character(coldata$conditions))){
  for (base2 in unique(as.character(coldata$conditions))){
    if (base1 != base2){
      comparison_key <- paste0(base1, "vs", base2)
      print(paste0(base1,",", base2))
      res <- results(dds, contrast=c("conditions", base1, base2))
     
      my_DE_list[[comparison_key]] <- res
      res_sorted = res[order(rownames(res)),]
      res_sorted$threshold <- as.logical(res_sorted$padj < 0.05)
      
      deseq2_results_res0.05 <- res_sorted[which(res_sorted$threshold == TRUE),]
      my_DE_list_0.05[[comparison_key]] <- deseq2_results_res0.05
      my_DE_list_0.05_up[[comparison_key]] <- deseq2_results_res0.05[which(deseq2_results_res0.05$log2FoldChange > 1),]
      my_DE_list_0.05_down[[comparison_key]] <- deseq2_results_res0.05[which(deseq2_results_res0.05$log2FoldChange < -1),]
      
      plotMA_temp3 <- ggmaplot(res_sorted,
                               fdr = 0.05, fc = 2, size = 0.3,
                               palette = c("#D22B2B", "#1465AC", "darkgray"),
                               genenames = rownames(res_sorted),
                               legend = "top", top = 20,
                               font.label = c("bold", 5),
                               font.legend = "bold",
                               font.main = "bold",
                               ggtheme = ggplot2::theme_minimal())
      plot_MA_list[[comparison_key]]  <- plotMA_temp3
    }
  }
}
```
#GENERATE VOLCANO PLOTS
```
#Contrasts of interest to create MA plots
#plotMA_subset <- c("Condition1_vs_Condition2", "Condition3_vs_Condition4")
plotMA_subset <- c(add your desired contrasts)

plot_MA_new <- list()
for(name1 in names(plot_MA_list)){
  for(name2 in plotMA_subset){
    if(name1 == name2){
      print(name1)
      plot_MA_new[[name1]] <- plot_MA_list[[name2]]
    }
  }
}
plotMA_new <- names(plot_MA_list)[names(plot_MA_list) %in% plotMA_subset]
pdf(file = "var1_var2_maplot_plot_list_rnaseq_noOutliers.pdf", height = 20, width = 20)
ggpubr::ggarrange(plotlist = plot_MA_new, ncol=5, nrow=5, common.legend = F, labels=names(plot_MA_new),
                  vjust = 1,hjust=-0.5,font.label = list(size = 10, color = "black", face = "bold", family = NULL))
dev.off()
```
#GENERATE LISTS OF UP AND DOWN REGULATED FILES
#You may upload these files in DAVID, shinyGO or other web based tools and explore different pathways, GO terms, etc.
```
# Save lists of genes to directory (the script you provided)
if (!dir.exists("gene_lists")) {
  dir.create("gene_lists")
}

# Save upregulated genes
for (contrast in plotMA_subset) {
  upregulated_genes <- annotation$locus_ID[match(rownames(my_DE_list_0.05_up[[contrast]]), annotation$gene_id)]
  write.csv(upregulated_genes, paste0("gene_lists/", contrast, "_up_genes.txt"), quote = FALSE, row.names = FALSE)
}

# Save downregulated genes
for (contrast in plotMA_subset) {
  downregulated_genes <- annotation$locus_ID[match(rownames(my_DE_list_0.05_down[[contrast]]), annotation$gene_id)]
  write.csv(downregulated_genes, paste0("gene_lists/", contrast, "_down_genes.txt"), quote = FALSE, row.names = FALSE)
}
```
#GENE SET ENRICHMENT ANALYSIS WITH GPROFILER
```
gprofiler_list_up <- list()
DE_gene_list_up <- list()
for (i in plotMA_subset){
  i_df_up <- data.frame(my_DE_list_0.05_up[[i]])
  print(i)
  i_df_up["gene_id"] <- rownames(i_df_up) 
  i_df_up_ann <- merge(i_df_up, annotation, by = "gene_id")
  i_df_up_ann_gene <- i_df_up_ann$locus_ID
  DE_gene_list_up[[i]] <- i_df_up_ann_gene
  gprofiler_list_up[[i]] <- gost(query = i_df_up_ann$locus_ID, organism = "organism_code", user_threshold = TRUE) # for example organism_code is "athaliana" for arabidopis
}
```
#GENERATING UPSET PLOTS
```
UpSetR::upset(fromList(DE_gene_list_up), order.by = "freq")
upset_up <- Upsetout(DE_gene_list_up)
```
#GENERATING HEATMAPS
#Here I have a desired set of genes and I would like to generate heatmap representing their expression from two different experiments 

#Get the zscores
```
# Normalizing the data and extracting the counts from the DESeq2 object
ddsNorm_exp1 <- estimateSizeFactors(dds)
ddsNormcounts_exp1 <- counts(ddsNorm_exp1, normalized=TRUE)
ddsNormcounts_exp1 <- data.frame(ddsNormcounts_exp1)

# Calculating the means of normalized counts for generic groups
column_patterns <- c("Group1_Var1", "Group1_Var2", "Group2_Var1", "Group2_Var2", "Treatment1_Var1", "Treatment1_Var2", "Treatment2_Var1", "Treatment2_Var2")
for (i in seq_along(column_patterns)) {
    ddsNormcounts_exp1_avg[[paste0("Avg_", column_patterns[i])]] <- rowMeans(ddsNormcounts_exp1[, grep(column_patterns[i], colnames(ddsNormcounts_exp1))])
}

# Removing the initial replicate columns and only keep calculated average values
columns_to_remove_exp1 <- grep("_\\d$", colnames(ddsNormcounts_exp1_avg), value = TRUE)
ddsNormcounts_exp1_avg <- ddsNormcounts_exp1_avg[ , !(colnames(ddsNormcounts_exp1_avg) %in% columns_to_remove_exp1)]

# Calculate z_score by centering and scaling and convert the data into data frame
z_TddsNormcounts_exp1 = scale(t(ddsNormcounts_exp1_avg), center = TRUE, scale = TRUE)
z_ddsNormcounts_exp1 <- data.frame(t(z_TddsNormcounts_exp1))
z_ddsNormcounts_exp1["gene"] <- rownames(z_ddsNormcounts_exp1)

```
#Generate heatmaps

```
# Create a data frame for gene IDs and symbols of interest
gene_ids_df <- data.frame(gene_id = c("list of gene ids that you want a heatmap for"),
                          gene_symbol = c("corresponding gene symbols"))

# Extracting corresponding identifiers
gene_ids_annotation <- merge(annotation, gene_ids_df, by.y = "gene_id", by.x = "locus_ID")

# Prepare data from Experiment 1
z_normalized_counts_exp1 <- merge(z_ddsNormcounts_exp1, gene_ids_annotation, by.x = "gene", by.y = "gene_id")
rownames(z_normalized_counts_exp1) <- z_normalized_counts_exp1$gene
z_normalized_counts_exp1 <- z_normalized_counts_exp1[, c(2:9,12)]  # Adjust columns as needed

# Import z-scores from Experiment 2
z_normalized_counts_exp2 <- read.csv("z_ddsNormcounts_exp2.csv")
z_normalized_counts_exp2 <- data.frame(z_normalized_counts_exp2)
rownames(z_normalized_counts_exp2) <- z_normalized_counts_exp2$X  # Assume 'X' is the column for gene IDs
z_normalized_counts_exp2 <- z_normalized_counts_exp2[,-1]  # Remove the ID column

# Merge data from Experiments 1 and 2
z_normalized_counts_combined <- merge(z_normalized_counts_exp1, z_normalized_counts_exp2, by = "gene_symbol")
rownames(z_normalized_counts_combined) <- z_normalized_counts_combined$gene_symbol
z_normalized_counts_combined <- z_normalized_counts_combined[,-1]  # Adjust columns as needed

# Specify breaks for color scaling in the heatmap
breaksList = seq(-1, 1, by = 0.1)

# Prepare sample annotations and color coding
my_sample_col1 <- data.frame(treatments = unlist(lapply(strsplit(colnames(z_normalized_counts_combined), "_"), function(x) x[1])))
rownames(my_sample_col1) <- colnames(z_normalized_counts_combined)
my_colour1 = list(treatments = c("Group1" = "green", "Group2" = "#00e6ac", "Group3" ="pink", "Group4" = "darkred", "Group5" = "grey", "Group6" = "black"))

# Define a generalized column order for clarity in presentation
ordered_columns <- c("Group1_Treatment1_avg", "Group2_Treatment1_avg", "Group1_Treatment2_avg", "Group2_Treatment2_avg", 
                     "Group3_Treatment1_avg", "Group4_Treatment1_avg", "Group3_Treatment2_avg", "Group4_Treatment2_avg", 
                     "Group5_Treatment1_avg", "Group6_Treatment1_avg", "Group5_Treatment2_avg", "Group6_Treatment2_avg")

z_normalized_counts_ordered <- z_normalized_counts_combined[, ordered_columns]

# Create the heatmap with ordered data
combined_genes_heatmap <- pheatmap(z_normalized_counts_ordered,
                                   color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
                                   breaks = breaksList,
                                   cluster_cols = FALSE, 
                                   cluster_rows = TRUE,
                                   fontsize = 12,
                                   clustering_distance_cols = "euclidean", 
                                   clustering_method = "ward.D", 
                                   border_color = "black",
                                   annotation_colors = my_colour1,
                                   annotation_col = my_sample_col1)

# Save the heatmap
ggsave(file = "generalized_gene_expression_heatmap.svg", plot = combined_genes_heatmap$gtable, width = 10, height = 8, units = "in")

```
