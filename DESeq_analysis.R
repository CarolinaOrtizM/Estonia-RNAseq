##### Differential expression analysis #########
###paired end reads #####
##

rm(list=ls(all=TRUE)); # graphics.off()

## requiered pacakges
library (DESeq2)
library (ggplot2)
library(VennDiagram)
library(ComplexHeatmap)

setwd("~/brain/Star_Alignment2")

#counts from STAR and subread #### 
Counts <- read.delim("CountTable.csv",header = T, row.names = 1, sep = ",")

Counts <-Counts[which(rowSums(Counts) >50),] 

Countsdf <- as.data.frame(Counts) 
write.csv(Countsdf, file = "Counts_all.csv")

#condition is the factor of my treatments
# C= cold, M= moderate (control), W = warm 

condition <- factor(c("C","C","C","C","M","M","M","M","W","W","W","W"))

coldata <- data.frame(row.names = colnames(Counts),condition)


dds <- DESeqDataSetFromMatrix(countData = Counts, colData = coldata, design = ~condition)
dds <- DESeq(dds)

vsdata <- vst(dds, blind = FALSE)

####some QC of the data 

plotPCA(vsdata, intgroup ="condition" , returnData = TRUE)

#To change the colors...
levels(vsdata$condition)
vsdata$condition <- factor(vsdata$condition, levels = c("W","M","C"))
vsdata$condition <- factor(vsdata$condition, levels = c("C","M","W"))

plotPCA(vsdata, intgroup ="condition")

#Customize PCA with ggplot

PCA <- plotPCA(vsdata, intgroup ="condition")

PCA + 
  scale_color_manual(values=c("C"="#0072b2", "M"="#009e73", "W"="#d55e00")) + 
  stat_ellipse(type="t", linetype=1, linewidth=1, level=0.95) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.title = element_blank(),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid.major = element_blank(),
        axis.ticks.length = unit(.25, "cm"),
        axis.ticks = element_line(color = "black"))
ggsave("PCA_DEGs.png", width=4, height=8)

##### results 

#### results by cold treatment with the "control" treatment that is moderate 

res1 <- results(dds, contrast = c("condition", "M", "C"))

#condition C vs M, 14893 rows

#results by Warm treatment vs moderate 

res2 <- results (dds, contrast = c("condition", "M", "W"))

#results by the two conditions cold and warm... being warm here the "important"

res3 <- results (dds, contrast = c("condition", "C", "W"))

#ommiting the NAs

sigs1 <- na.omit(res1)
sigs2 <- na.omit(res2)
sigs3 <- na.omit(res3)

#only the genes with the padjust <0.05 significance

sigs1 <- sigs1[sigs1$padj < 0.05,]

#186 diff expressed genes for cold compared to moderate

sigs2 <- sigs2[sigs2$padj < 0.05,]

#2.219 diff exp genes for moderate compared with warm 

sigs3 <- sigs3[sigs3$padj < 0.05,]

#1691 for warm, compared with cold 

#Save the data
write.csv(sigs1, file = "deseq_results_MvC.csv") 
write.csv(sigs2, file = "deseq_results_MvW.csv")
write.csv(sigs3, file = "deseq_results_CvW.csv")

##### Venn diagram for the diff expressed genes 
str(sigs1)
class(sigs1)

# Convert DESeqResults object to a dataframe... use these ones for GO analysis
df_sigs1 <- as.data.frame(sigs1)      
df_sigs2 <-as.data.frame(sigs2)
df_sigs3 <-as.data.frame(sigs3)

#getting the numbers of up and down regulated  
head(df_sigs1)
str(df_sigs1)

# Count of up-regulated genes per group of DEGs
num_up_regulated <- sum(df_sigs3$log2FoldChange > 0)

# Count of down-regulated genes
num_down_regulated <- sum(df_sigs3$log2FoldChange < 0)

# Print the results
cat("Number of up-regulated genes:", num_up_regulated, "\n")
cat("Number of down-regulated genes:", num_down_regulated, "\n")

#extract the gene id from the data frames 
genes_cold = row.names(df_sigs1)
genes_moderate = row.names(df_sigs2)
genes_warm = row.names(df_sigs3)

# Create a Venn diagram
genes_cold = as.character(genes_cold)
genes_moderate = as.character(genes_moderate)
genes_warm = as.character(genes_warm)

# Calculate intersections
intersect_cm = base::intersect(genes_cold, genes_moderate)
intersect_cw = base::intersect(genes_cold, genes_warm)
intersect_mw = base::intersect(genes_moderate, genes_warm)

# Intersection of all three
intersect_all = base::intersect(intersect_cm, genes_warm)

# Calculate unique counts for the Venn diagram
only_cold = setdiff(genes_cold, union(intersect_cm, intersect_cw))
only_moderate = setdiff(genes_moderate, union(intersect_cm, intersect_mw))
only_warm = setdiff(genes_warm, union(intersect_cw, intersect_mw))
n12 = length(intersect_cm) - length(intersect_all)  # Shared between Cold and Moderate, not including all three
n23 = length(intersect_mw) - length(intersect_all)  # Shared between Moderate and Warm, not including all three
n13 = length(intersect_cw) - length(intersect_all)  # Shared between Cold and Warm, not including all three
n123 = length(intersect_all)  # Shared among all

# Setup for the Venn diagram with previously calculated intersections
venn.plot <- draw.triple.venn(
  area1 = length(genes_cold),
  area2 = length(genes_moderate),
  area3 = length(genes_warm),
  n12 = n12,
  n23 = n23,
  n13 = n13,
  n123 = n123,
  category = c("Cold", "Moderate", "Warm"),
  fill = c("#0072b2", "#009e73", "#d55e00"),
  cex = 1.8,  # Adjust text size inside Venn
  cat.cex = 1.6,  # Adjust category text size
  fontfamily ="sans serif",
  cat.fontfamily = "sans serif",
  cat.pos = c(-20, 20, 180),  # Adjust these angles to position labels better
  cat.dist = c(0.04, 0.04, 0.04) #Adjust the distance of labels from the center of the circles
)

# Start a new PNG device to save the plot
png("VennDiagram.png")

# Draw the plot on the new device
grid.draw(venn.plot)

# Close the PNG device to finalize the file
dev.off()

######Fisher's Exact Test 

##########moderate vs Cold 

# Intersections
intersect_cm = intersect(genes_cold, genes_moderate)
intersect_cw = intersect(genes_cold, genes_warm)
intersect_mw = intersect(genes_moderate, genes_warm)
intersect_all = intersect(intersect(intersect_cm, genes_warm), genes_warm)

# Unique counts
unique_cold = setdiff(genes_cold, union(intersect_cm, intersect_cw))
unique_moderate = setdiff(genes_moderate, union(intersect_cm, intersect_mw))
unique_warm = setdiff(genes_warm, union(intersect_cw, intersect_mw))

# Total genes (replace with actual total if known)
total_genes = length(union(union(genes_cold, genes_moderate), genes_warm))

# Contingency table for Cold vs. Moderate
cold_moderate_matrix <- matrix(c(
  length(intersect_cm) - length(intersect_all),  # Shared between Cold and Moderate, not in Warm
  length(unique_cold),                          # Only in Cold
  length(unique_moderate),                      # Only in Moderate
  total_genes - (length(intersect_cm) + length(unique_cold) + length(unique_moderate) - length(intersect_all))
), nrow = 2, byrow = TRUE)

# Naming the dimensions for clarity
dimnames(cold_moderate_matrix) <- list(
  "In Cold" = c("In Moderate", "Not in Moderate"),
  "Not in Cold" = c("In Moderate", "Not in Moderate")
)

# Perform Fisher's Exact Test
fishers_test_result = fisher.test(cold_moderate_matrix)
print(fishers_test_result)

######### cold vs warm

# Intersections
intersect_cw = intersect(genes_cold, genes_warm)
intersect_all = intersect(intersect(intersect(genes_cold, genes_moderate), genes_warm), genes_warm)

# Unique counts
unique_cold = setdiff(genes_cold, union(intersect_cw, intersect(genes_moderate, genes_warm)))
unique_warm = setdiff(genes_warm, union(intersect_cw, intersect(genes_cold, genes_moderate)))

# Assuming you know the total number of genes, update the total_genes variable
total_genes = length(union(union(genes_cold, genes_moderate), genes_warm))

# Contingency table for Cold vs. Warm
cold_warm_matrix <- matrix(c(
  length(intersect_cw) - length(intersect_all),
  length(unique_cold),
  length(unique_warm),
  total_genes - (length(intersect_cw) + length(unique_cold) + length(unique_warm) - length(intersect_all))
), nrow = 2, byrow = TRUE)

# Naming the dimensions for clarity
dimnames(cold_warm_matrix) <- list(
  "In Cold" = c("In Warm", "Not in Warm"),
  "Not in Cold" = c("In Warm", "Not in Warm")
)

# Perform Fisher's Exact Test
fishers_test_cw = fisher.test(cold_warm_matrix)
print(fishers_test_cw)
#0.0067 significant! 

###### warm vs moderate 

# Intersections
intersect_mw = intersect(genes_moderate, genes_warm)
intersect_all = intersect(intersect(intersect(genes_cold, genes_moderate), genes_warm), genes_warm)

# Unique counts
unique_moderate = setdiff(genes_moderate, union(intersect_mw, intersect(genes_cold, genes_warm)))
unique_warm = setdiff(genes_warm, union(intersect_mw, intersect(genes_cold, genes_moderate)))

# Contingency table for Moderate vs. Warm
moderate_warm_matrix <- matrix(c(
  length(intersect_mw) - length(intersect_all),
  length(unique_moderate),
  length(unique_warm),
  total_genes - (length(intersect_mw) + length(unique_moderate) + length(unique_warm) - length(intersect_all))
), nrow = 2, byrow = TRUE)

# Naming the dimensions for clarity
dimnames(moderate_warm_matrix) <- list(
  "In Moderate" = c("In Warm", "Not in Warm"),
  "Not in Moderate" = c("In Warm", "Not in Warm")
)

# Perform Fisher's Exact Test
fishers_test_mw = fisher.test(moderate_warm_matrix)
print(fishers_test_mw)
#< 2.2e-16 #highly significant 

#########################################
#simple heatmap 

#get sigs as data frame

sigs1.df <- as.data.frame(sigs1)

sigs2.df <- as.data.frame(sigs2)

sigs3.df <- as.data.frame (sigs3)

#cluster genes and samples, it groups similar patterns together 
#actual heatmap

##### cold vs moderate sigs1.df

sigs1.df <- sigs1.df[(sigs1.df$baseMean > 100) & (abs(sigs1.df$log2FoldChange) > 1.5),]

mat1 <- counts(dds, normalized = T)[rownames(sigs1.df),]
mat1.z <- t(apply(mat1, 1, scale))
colnames(mat1.z) <- rownames(coldata)
mat1.z

#only cold vs moderate
cold_moderate_cols <- colnames(mat1.z)[grepl("Cold|Moderate", colnames(mat1.z))]

# Subset mat1.z for Cold vs Moderate
mat1_cm <- mat1.z[, cold_moderate_cols]
colnames(mat1_cm) <- cold_moderate_cols

# Generate the heatmap
h11 <- Heatmap(mat1_cm, cluster_rows = TRUE, cluster_columns = TRUE, column_labels = colnames(mat1_cm), name = "Z-score")
print(h11)

#### Warm vs moderate sigs2.df 

sigs2.df <- sigs2.df[(sigs2.df$baseMean > 100) & (abs(sigs2.df$log2FoldChange) > 1.5),]

mat2 <- counts(dds, normalized = T)[rownames(sigs2.df),]
mat2.z <- t(apply(mat2, 1, scale))
colnames(mat2.z) <- rownames(coldata)
mat2.z

h2 <- Heatmap(mat2.z, cluster_rows = T, cluster_columns = T, column_labels = colnames(mat2.z), name = "Z-score")

#only moderate and warm samples
warm_moderate_cols <- colnames(mat2.z)[grepl("Warm|Moderate", colnames(mat2.z))]

# Subset mat2.z for Warm vs Moderate
mat2_wm <- mat2.z[, warm_moderate_cols]
colnames(mat2_wm) <- warm_moderate_cols

# Generate the heatmap
h22 <- Heatmap(mat2_wm, cluster_rows = TRUE, cluster_columns = TRUE, column_labels = colnames(mat2_wm), name = "Z-score")
print(h22)

##### cold vs warm 

sigs3.df <- sigs3.df[(sigs3.df$baseMean > 100) & (abs(sigs3.df$log2FoldChange) > 1.5),]

mat3 <- counts(dds, normalized = T)[rownames(sigs3.df),]
mat3.z <- t(apply(mat3, 1, scale))
colnames(mat3.z) <- rownames(coldata)
mat3.z

h3 <- Heatmap(mat3.z, cluster_rows = T, cluster_columns = T, column_labels = colnames(mat3.z), name = "Z-score")

#only the cold and warm samples 

cold_warm_cols <- colnames(mat3.z)[grepl("Cold|Warm", colnames(mat3.z))]

# Subset mat3.z for Cold vs Warm
mat3_cw <- mat3.z[, cold_warm_cols]
colnames(mat3_cw) <- cold_warm_cols

# Generate the heatmap
h33 <- Heatmap(mat3_cw, cluster_rows = TRUE, cluster_columns = TRUE, column_labels = colnames(mat3_cw), name = "Z-score")
print(h33)

#########
###### vulcano plots

library(EnhancedVolcano)
citation("EnhancedVolcano")

### cold vs moderate
res1 <- res1[res1$baseMean > 10,]
res1

res1.df <- as.data.frame(res1)

# Define thresholds for significance
pvalueCutoff <- 1e-4
log2FCCutoff <- 1

# Remove rows with any NA values in the specified columns to avoid errors in computations
res1.df <- res1.df[!is.na(res1.df$log2FoldChange) & !is.na(res1.df$padj), ]

# Create a color mapping based on conditions
res1.df$Color <- ifelse(res1.df$log2FoldChange > log2FCCutoff & res1.df$padj < pvalueCutoff, 'red',
                        ifelse(res1.df$log2FoldChange < -log2FCCutoff & res1.df$padj < pvalueCutoff, 'blue', 'grey50'))

# Calculate the number of up-regulated and down-regulated genes
num_upregulated <- sum(res1.df$log2FoldChange > log2FCCutoff & res1.df$padj < pvalueCutoff, na.rm = TRUE)
num_downregulated <- sum(res1.df$log2FoldChange < -log2FCCutoff & res1.df$padj < pvalueCutoff, na.rm = TRUE)

# Create the plot with lines and text annotations
volcano_plot1 <- ggplot(res1.df, aes(x = log2FoldChange, y = -log10(padj), color = Color)) +
  geom_point(alpha = 0.5) +
  scale_color_identity() +
  geom_vline(xintercept = log2FCCutoff, linetype = "dashed", color = "black") +
  geom_vline(xintercept = -log2FCCutoff, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(pvalueCutoff), linetype = "dashed", color = "black") +
  geom_text(aes(x = 3, y = -log10(pvalueCutoff) + 5, label = paste(num_upregulated, "genes upregulated")), vjust = 1, hjust = 1, color = "red") +
  geom_text(aes(x = -3, y = -log10(pvalueCutoff) + 5, label = paste(num_downregulated, "genes downregulated")), vjust = 1, hjust = 0, color = "blue") +
  theme_classic() +
  labs(x = "Log2 Fold Change", y = "-Log10 adjusted P-value") +  # Removed title
  theme(legend.position = "none")  # Remove the legend
# Print the plot
print(volcano_plot1)


tiff("VolcanoPlotgg_MvC.tiff", res = 250, width = 1300, height = 1500)
print(volcano_plot1)
dev.off()

########## In a volcano plot, the most upregulated genes are towards the right, the most downregulated genes are towards the left,
#and the most statistically significant genes are towards the top.

### warm vs moderate

res2 <- res2[res2$baseMean > 50,]
res2

res2.df <- as.data.frame(res2)
#EnhancedVolcano(res2.df, x = "log2FoldChange", y = "padj", lab = rownames(res2.df))
#V2 <- EnhancedVolcano(res2.df, x = "log2FoldChange", y = "padj", lab = rownames(res2.df), pCutoff = 1e-4, FCcutoff = 1)

# Define thresholds for significance
pvalueCutoff <- 1e-4
log2FCCutoff <- 1

# Remove rows with any NA values in the specified columns to avoid errors in computations
res2.df <- res2.df[!is.na(res2.df$log2FoldChange) & !is.na(res2.df$padj), ]

# Create a color mapping based on conditions
res2.df$Color <- ifelse(res2.df$log2FoldChange > log2FCCutoff & res2.df$padj < pvalueCutoff, 'red',
                        ifelse(res2.df$log2FoldChange < -log2FCCutoff & res2.df$padj < pvalueCutoff, 'blue', 'grey50'))

# Calculate the number of up-regulated and down-regulated genes
num_upregulated <- sum(res2.df$log2FoldChange > log2FCCutoff & res2.df$padj < pvalueCutoff, na.rm = TRUE)
num_downregulated <- sum(res2.df$log2FoldChange < -log2FCCutoff & res2.df$padj < pvalueCutoff, na.rm = TRUE)

# Create the plot with lines and text annotations
volcano_plot2 <- ggplot(res2.df, aes(x = log2FoldChange, y = -log10(padj), color = Color)) +
  geom_point(alpha = 0.5) +
  scale_color_identity() +
  geom_vline(xintercept = log2FCCutoff, linetype = "dashed", color = "black") +
  geom_vline(xintercept = -log2FCCutoff, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(pvalueCutoff), linetype = "dashed", color = "black") +
  geom_text(aes(x = 5, y = -log10(pvalueCutoff) + 5, label = paste(num_upregulated, "genes upregulated")), vjust = 1, hjust = 1, color = "red") +
  geom_text(aes(x = -5, y = -log10(pvalueCutoff) + 5, label = paste(num_downregulated, "genes downregulated")), vjust = 1, hjust = 0, color = "blue") +
  theme_classic() +
  labs(x = "Log2 Fold Change", y = "-Log10 adjusted P-value") +  # Removed title
  theme(legend.position = "none")  # Remove the legend
# Print the plot
print(volcano_plot2)


tiff("VolcanoPlotgg_WvM.tiff", res = 250, width = 1300, height = 1500)
print(volcano_plot2)
dev.off()


#### Cold vs warm

res3 <- res3[res3$baseMean > 50,]

res3.df <- as.data.frame(res3)
#EnhancedVolcano(res3.df, x = "log2FoldChange", y = "padj", lab = rownames(res3.df))
#V3 <- EnhancedVolcano(res3.df, x = "log2FoldChange", y = "padj", lab = rownames(res3.df), pCutoff = 1e-4, FCcutoff = 1)

# Define thresholds for significance
pvalueCutoff <- 1e-4
log2FCCutoff <- 1

# Remove rows with any NA values in the specified columns to avoid errors in computations
res3.df <- res3.df[!is.na(res3.df$log2FoldChange) & !is.na(res3.df$padj), ]

# Create a color mapping based on conditions
res3.df$Color <- ifelse(res3.df$log2FoldChange > log2FCCutoff & res3.df$padj < pvalueCutoff, 'red',
                        ifelse(res3.df$log2FoldChange < -log2FCCutoff & res3.df$padj < pvalueCutoff, 'blue', 'grey50'))

# Calculate the number of up-regulated and down-regulated genes
num_upregulated <- sum(res3.df$log2FoldChange > log2FCCutoff & res3.df$padj < pvalueCutoff, na.rm = TRUE)
num_downregulated <- sum(res3.df$log2FoldChange < -log2FCCutoff & res3.df$padj < pvalueCutoff, na.rm = TRUE)

# Create the plot with lines and text annotations
volcano_plot <- ggplot(res3.df, aes(x = log2FoldChange, y = -log10(padj), color = Color)) +
  geom_point(alpha = 0.5) +
  scale_color_identity() +
  geom_vline(xintercept = log2FCCutoff, linetype = "dashed", color = "black") +
  geom_vline(xintercept = -log2FCCutoff, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(pvalueCutoff), linetype = "dashed", color = "black") +
  geom_text(aes(x = 6, y = -log10(pvalueCutoff) + 10, label = paste(num_upregulated, "genes upregulated")), vjust = 1, hjust = 1, color = "red") +
  geom_text(aes(x = -6, y = -log10(pvalueCutoff) + 10, label = paste(num_downregulated, "genes downregulated")), vjust = 1, hjust = 0, color = "blue") +
  theme_classic() +
  labs(x = "Log2 Fold Change", y = "-Log10 adjusted P-value") +  # Removed title
  theme(legend.position = "none")  # Remove the legend
# Print the plot
print(volcano_plot)

tiff("VolcanoPlotgg_CvW.tiff", res = 250, width = 1300, height = 1500)
print(volcano_plot)
dev.off()

#################
########## 
#####ontology of genes (GO) (see the GO script)
