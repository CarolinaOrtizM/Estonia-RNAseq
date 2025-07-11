##### Differential expression analysis #########
###paired end reads #####
##

rm(list=ls(all=TRUE)); # graphics.off()

## requiered pacakges
library (DESeq2)
library (ggplot2)
library(VennDiagram)
library(ComplexHeatmap)

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

#Customize PCA with ggplot

PCA <- plotPCA(vsdata, intgroup ="condition")

###with box and labels 

PCA +
  scale_color_manual(values = c("C" = "#0072b2", "M" = "#009e73", "W" = "#d55e00")) +
  geom_text_repel(  # Replaces geom_text() with ggrepel's version
    aes(label = name),  # Assuming 'name' is the column with labels
    size = 6,          # Adjust font size
    box.padding = 0.5,  # Controls spacing around labels
    max.overlaps = 20,  # Increase if too many labels are hidden
    segment.color = "grey50",  # Color of label connector lines
    min.segment.length = 0.2   # Shortens connector lines
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16),
    axis.ticks.length = unit(.25, "cm"),
    axis.text.y = element_text(size = 16, color = "black"),
    axis.text.x = element_text(size = 16, color = "black"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)  # Optional: adds full border
  )

##### results #######

#### results by Cold treatment with the "control" treatment that is "Moderate" 

res1 <- results(dds, contrast = c("condition", "C", "M"))

#results by Warm treatment vs moderate 

res2 <- results (dds, contrast = c("condition", "W", "M"))

#results by warm vs cold 

res3 <- results (dds, contrast = c("condition", "W", "C"))

# omitting the NAs

sigs1 <- na.omit(res1)
sigs2 <- na.omit(res2)
sigs3 <- na.omit(res3)

#only the genes with the p-adjust <0.05 significance 

sigs1.1 <- sigs1[sigs1$padj < 0.05 & (sigs1$log2FoldChange > 1 | sigs1$log2FoldChange < -1), ]
sigs1 <- sigs1 [sigs1$padj < 0.05, ]

#186 diff expressed genes for cold compared to moderate (only padj)
#43 diff expressed with log2 >< 1

sigs2.2 <- sigs2[sigs2$padj < 0.05 & (sigs2$log2FoldChange > 1 | sigs2$log2FoldChange < -1),] 
sigs2 <- sigs2 [sigs2$padj < 0.05, ]

#2.219 diff exp genes for warm compared with moderate (only padj)
#453 diff expressed with log2 >< 1

sigs3.3 <- sigs3[sigs3$padj < 0.05 & (sigs3$log2FoldChange > 1 | sigs3$log2FoldChange < -1),] 
sigs3 <- sigs3 [sigs3$padj < 0.05, ]

#1691 for warm, compared with cold (only padj)
#422 diff expressed with log2 >< 1

#Save the data
write.csv(sigs1, file = "deseq_results_CvM.csv") 
write.csv(sigs2, file = "deseq_results_WvM.csv")
write.csv(sigs3, file = "deseq_results_WvC.csv")

##### Venn diagram for the diff expressed genes ########
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
  lwd = c(2, 2, 2),  # Border line widths
  col = c(adjustcolor("black", alpha.f = 0.5),
          adjustcolor("black", alpha.f = 0.5),
          adjustcolor("black", alpha.f = 0.5)),
  category = c("cold/moderate (189)", "warm/moderate (2219)","warm/cold (1691)"),
  fill = c("#0072b2", "#d55e00","#8E44AD"),
  cex = 1.9,  # Adjust text size inside Venn
  cat.cex = 1.8,  # Adjust category text size
  fontfamily ="sans serif",
  cat.fontfamily = "sans serif",
  cat.pos = c(-10, 10, 180),  # Adjust these angles to position labels better
  cat.dist = c(0.04, 0.04, 0.04) #Adjust the distance of labels from the center of the circles
)

grid.newpage()
pushViewport(viewport(width = unit(1, "snpc"), height = unit(1, "snpc")))
grid.draw(venn.plot)
popViewport() 

######Fisher's Exact Test ###########

########## Cold vs moderate  

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

######### warm vs cold 

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

# cluster genes and samples, it groups similar patterns together

#get sigs as data frame

sigs1.df <- as.data.frame(sigs1)

sigs2.df <- as.data.frame(sigs2)

sigs3.df <- as.data.frame (sigs3)


##### Cold vs moderate 

sigs1.df <- sigs1.df[(sigs1.df$baseMean > 100) & (abs(sigs1.df$log2FoldChange) > 1),]

mat1 <- counts(dds, normalized = T)[rownames(sigs1.df),]
mat1.z <- t(apply(mat1, 1, scale))
colnames(mat1.z) <- rownames(coldata)
mat1.z

h1 <- Heatmap(mat1.z, cluster_rows = T, cluster_columns = T, column_labels = colnames(mat1.z), name = "Z-score")

## to save the figure in a nice way 
png("simple_heatmap111.png", res = 250, width = 1000, height = 1500)
print(h1)
dev.off() 

#only cold vs moderate samples 
cold_moderate_cols <- colnames(mat1.z)[grepl("Cold|Moderate", colnames(mat1.z))]

# Subset mat1.z for Cold vs Moderate
mat1_cm <- mat1.z[, cold_moderate_cols]
colnames(mat1_cm) <- cold_moderate_cols

# Generate the heatmap
h11 <- Heatmap(mat1_cm, cluster_rows = TRUE, cluster_columns = TRUE, column_labels = colnames(mat1_cm), name = "Z-score")
h11

png("heatmapColdvsModerate.png", res = 250, width = 1000, height = 1500)
print(h11)
dev.off()

#### Warm vs moderate 

sigs2.df <- sigs2.df[(sigs2.df$baseMean > 100) & (abs(sigs2.df$log2FoldChange) > 2.5),]

mat2 <- counts(dds, normalized = T)[rownames(sigs2.df),]
mat2.z <- t(apply(mat2, 1, scale))
colnames(mat2.z) <- rownames(coldata)
mat2.z

#h2 <- Heatmap(mat2.z, cluster_rows = T, cluster_columns = T, column_labels = colnames(mat2.z), name = "Z-score")

#png("simple_heatmap222.png", res = 250, width = 1000, height = 1500)
#print(h2)
#dev.off() 

#only moderate and warm samples
warm_moderate_cols <- colnames(mat2.z)[grepl("Warm|Moderate", colnames(mat2.z))]

# Subset mat2.z for Warm vs Moderate
mat2_wm <- mat2.z[, warm_moderate_cols]
colnames(mat2_wm) <- warm_moderate_cols

# Generate the heatmap
h22 <- Heatmap(mat2_wm, cluster_rows = TRUE, cluster_columns = TRUE, column_labels = colnames(mat2_wm), name = "Z-score")
h22

png("HeatmapWarm_mode.png", res = 250, width = 1000, height = 1800)
h22
dev.off()

##### warm vs cold 

sigs3.df <- sigs3.df[(sigs3.df$baseMean > 100) & (abs(sigs3.df$log2FoldChange) > 1.5),]

mat3 <- counts(dds, normalized = T)[rownames(sigs3.df),]
mat3.z <- t(apply(mat3, 1, scale))
colnames(mat3.z) <- rownames(coldata)
mat3.z

h3 <- Heatmap(mat3.z, cluster_rows = T, cluster_columns = T, column_labels = colnames(mat3.z), name = "Z-score")
h3


#png("simple_heatmap333.png", res = 250, width = 1000, height = 1500)
#print(h3)
#dev.off()

#only the cold and warm samples 

cold_warm_cols <- colnames(mat3.z)[grepl("Cold|Warm", colnames(mat3.z))]

# Subset mat3.z for Cold vs Warm
mat3_cw <- mat3.z[, cold_warm_cols]
colnames(mat3_cw) <- cold_warm_cols

# Generate the heatmap
h33 <- Heatmap(mat3_cw, cluster_rows = TRUE, cluster_columns = TRUE, column_labels = colnames(mat3_cw), name = "Z-score")
h33

png("heatmapColdvsWarm.png", res = 250, width = 1000, height = 1500)
print(h33)
dev.off()

############### counts per genes ########

head(res3[order(res3$padj),], 10)

countx <- plotCounts(dds, which.min(res1$padj), 
                     intgroup = c("condition"), returnData = TRUE)

######QORs ######
#Quinone oxidoreductases (QORs) are a group of enzymes that play a crucial role in detoxification 
#and energy production in various organisms

# Custom colors for conditions
custom_colors <- c("C" = "#0072b2", "M" = "#009e73", "W" = "#d55e00")

QORs<-ggplot(countx, aes(x = condition, y = count, 
                   color = condition, fill = condition, 
                   group = condition)) + 
  #geom_boxplot(width = 0.1, alpha = 0.2, outlier.shape = NA) +  # Add boxplot with transparency
  geom_quasirandom(size = 2, alpha = 0.7, stroke = 1.2, width = 0.1) +  # Jitter with controlled width
  scale_y_log10(limits = c(2000, 5000), breaks = c(1000, 2000, 3000, 4000, 5000))+
  scale_color_manual(values = custom_colors) +  # Apply colors to shape borders
  scale_fill_manual(values = custom_colors) +  # Apply colors to shape fills
  theme_classic() +  # Remove gray background
  ggtitle("quinone oxidoreductases", "LOC129981218, QORs") +
  theme(
    axis.title.x = element_blank(),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.text.x = element_text(size = 11, color = "black"),
    legend.position = "none"
  )
QORs

kruskal.test(countx$count ~ countx$condition)
#Kruskal-Wallis chi-squared = 7.5385, df = 2, p-value = 0.02307

dunnTest(countx$count ~ countx$condition, method="bh")

####NRF6 ####
#lipid transport 
#Plays a role in the uptake of a range of molecules, including lipids and xenobiotic compounds 
#from the intestine to surrounding tissues. Mediates the  transport of lipids from the intestine to the reproductive tract. 
#Required for efficient yolk transport into oocytes. Vital for embryonic development.

counts1 <- plotCounts(dds, "LOC129962681", 
                      intgroup = c("condition"), returnData = TRUE)

ggplot(counts1, aes(x = condition, y = count, 
                   color = condition, fill = condition, 
                   group = condition)) + 
  geom_boxplot(width = 0.1, alpha = 0.2, outlier.shape = NA) +  # Add boxplot with transparency
  geom_quasirandom(size = 2, alpha = 0.7, stroke = 1.2, width = 0.1) +  # Jitter with controlled width
  scale_color_manual(values = custom_colors) +  # Apply colors to shape borders
  scale_fill_manual(values = custom_colors) +  # Apply colors to shape fills
  theme_classic() +  # Remove gray background
  ggtitle("Nose resistant fluoxetine", "Gene: 129962681, NRF6") +
  theme(
    axis.title.x = element_blank(),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.text.x = element_text(size = 11, color = "black"),
    legend.position = "none"
  )
#######apoLp ####
#apolipophorins-like
#protein components of lipoprotein particles that are essential for lipid transportation in the insect body.
counts2 <- plotCounts(dds, "LOC129959427", 
                      intgroup = c("condition"), returnData = TRUE)

apoLp<-ggplot(counts2, aes(x = condition, y = count, 
                    color = condition, fill = condition, 
                    group = condition)) + 
  #geom_boxplot(width = 0.1, alpha = 0.2, outlier.shape = NA) +  # Add boxplot with transparency
  geom_quasirandom(size = 2, alpha = 0.7, stroke = 1.2, width = 0.1) +  # Jitter with controlled width
  scale_color_manual(values = custom_colors) +  # Apply colors to shape borders
  scale_fill_manual(values = custom_colors) +  # Apply colors to shape fills
  theme_classic() +  # Remove gray background
  ggtitle("apolipophorins-like", "LOC129959427, apoLp") +
  theme( 
    axis.title.x = element_blank(),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.text.x=element_blank(), 
    axis.ticks.x=element_blank(),
    legend.position = "none"
  )
apoLp

head(res3[order(res3$padj),], 10)

######GPR52####
#crucial signaling molecules in arthropods, playing roles in sensing environmental cues, regulating metabolic processes
#essential transmembrane proteins in arthropods that detect various signals and activate cellular responses.
#growth, development, stress responses, and behaviors. 

counts3 <- plotCounts(dds, "LOC129963668", 
                      intgroup = c("condition"), returnData = TRUE)

ggplot(counts3, aes(x = condition, y = count, 
                    color = condition, fill = condition, 
                    group = condition)) + 
  geom_boxplot(width = 0.1, alpha = 0.2, outlier.shape = NA) +  # Add boxplot with transparency
  geom_quasirandom(size = 2, alpha = 0.7, stroke = 1.2, width = 0.1) +  # Jitter with controlled width
  scale_color_manual(values = custom_colors) +  # Apply colors to shape borders
  scale_fill_manual(values = custom_colors) +  # Apply colors to shape fills
  theme_classic() +  # Remove gray background
  ggtitle("G-protein coupled receptor 52-like", "Gene: 129963668, GPR52") +
  theme(
    axis.title.x = element_blank(),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.text.x = element_text(size = 11, color = "black"),
    legend.position = "none"
  )

#####PLA2####
#Phospholipase A2
#superfamily of enzymes that hydrolyze phospholipids, releasing free fatty acids and lysophospholipids
#catalyses the cleavage of fatty acids in position 2 of phospholipids, 
#hydrolyzing the bond between the second fatty acid "tail" and the glycerol molecule
#releasing arachidonic acid and lysophosphatidyl choline, a precursor of lysophosphatidic acid. 
#downstream modification by cyclooxygenases or lipoxygenases, 
#arachidonic acid is modified into active compounds called eicosanoids -->> categorized as anti-inflammatory and inflammatory mediators

#129960889 another one 
#129964053

counts4 <- plotCounts(dds, "LOC129964053", 
                      intgroup = c("condition"), returnData = TRUE)

PLA2<-ggplot(counts4, aes(x = condition, y = count, 
                    color = condition, fill = condition, 
                    group = condition)) + 
  geom_boxplot(width = 0.1, alpha = 0.2, outlier.shape = NA) +  # Add boxplot with transparency
  geom_quasirandom(size = 2, alpha = 0.7, stroke = 1.2, width = 0.1) +  # Jitter with controlled width
  scale_color_manual(values = custom_colors) +  # Apply colors to shape borders
  scale_fill_manual(values = custom_colors) +  # Apply colors to shape fills
  theme_classic() +  
  ggtitle("Phospholipase A2", "Gene: 129964053, PLA2") +
    theme(
      axis.title.x = element_blank(),
      axis.title = element_text(size = 14, face = "bold"),
      axis.text.y = element_text(size = 11, color = "black"),
      axis.text.x = element_text(size = 11, color = "black"),
      legend.position = "none"
    )
PLA2

kruskal.test(counts4$count ~ counts4$condition)
#Kruskal-Wallis chi-squared = 9.2692, df = 2, p-value = 0.00971

dunnTest(counts4$count ~ counts4$condition, method="bh")

######LPAT####
#lysophospholipid acyltransferase
# involved in phospholipid remodeling, specifically incorporating fatty acyl chains into lysophospholipids.
#contributing to the dynamic balance of cellular lipid composition

counts5 <- plotCounts(dds, "LOC129987920", 
                      intgroup = c("condition"), returnData = TRUE)

LPAT<-ggplot(counts5, aes(x = condition, y = count, 
                    color = condition, fill = condition, 
                    group = condition)) + 
  #geom_boxplot(width = 0.1, alpha = 0.2, outlier.shape = NA) +  # Add boxplot with transparency
  geom_quasirandom(size = 2, alpha = 0.7, stroke = 1.2, width = 0.1) +  # Jitter with controlled width
  scale_y_log10(limits = c(800, 2000), breaks = c(800, 1000, 1200, 1500, 2000))+
  scale_color_manual(values = custom_colors) +  # Apply colors to shape borders
  scale_fill_manual(values = custom_colors) +  # Apply colors to shape fills
  theme_classic() +  # Remove gray background
  ggtitle("lysophospholipid acyltransferase", "LOC129987920, LPAT") +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.text.x=element_blank(), 
    axis.ticks.x=element_blank(),
    legend.position = "none"
  )
LPAT

kruskal.test(counts5$count ~ counts5$condition)
#Kruskal-Wallis chi-squared = 9.8462, df = 2, p-value = 0.007277

dunnTest(counts5$count ~ counts5$condition, method="bh")


####AACS####
#Acetoacetyl-CoA synthetase
#an enzyme that converts acetoacetate to acetoacetyl-CoA in the cytosol. This process is crucial for utilizing ketone bodies (KBs) in lipogenesis, 
#the synthesis of lipids. AACS plays a role in the synthesis of cholesterol and fatty acids.
#lipogenic enzyme found in various tissues

counts6 <- plotCounts(dds, "LOC129964270", 
                      intgroup = c("condition"), returnData = TRUE)
AACS<-ggplot(counts6, aes(x = condition, y = count, 
                     color = condition, fill = condition, 
                     group = condition)) + 
  #geom_boxplot(width = 0.1, alpha = 0.2, outlier.shape = NA) +  # Add boxplot with transparency
  geom_quasirandom(size = 2, alpha = 0.7, stroke = 1.2, width = 0.1) +  # Jitter with controlled width
  scale_y_log10(limits = c(30, 2000), breaks = c(30, 50, 100, 500, 1000, 1500, 2000))+
  scale_color_manual(values = custom_colors) +  # Apply colors to shape borders
  scale_fill_manual(values = custom_colors) +  # Apply colors to shape fills
  theme_classic() +  
  ggtitle("acetoacetyl-CoA synthetase", "LOC129964270, AACS") +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    #axis.title = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.text.x = element_text(size = 11, color = "black"),
    legend.position = "none"
  )
AACS
kruskal.test(counts6$count ~ counts6$condition)
#Kruskal-Wallis chi-squared = 8, df = 2, p-value = 0.01832

dunnTest(counts6$count ~ counts6$condition, method="bh")

####FAAH2 ####
#Fatty acid amide hydrolase 2
#is an enzyme primarily known for its role in metabolizing fatty acid amides, 
#including anandamide (AEA), an endocannabinoid, and related amidated signaling lipids
counts7 <- plotCounts(dds, "LOC129976441", 
                      intgroup = c("condition"), returnData = TRUE)

ggplot(counts7, aes(x = condition, y = count, 
                    color = condition, fill = condition, 
                    group = condition)) + 
  geom_boxplot(width = 0.1, alpha = 0.2, outlier.shape = NA) +  # Add boxplot with transparency
  geom_quasirandom(size = 2, alpha = 0.7, stroke = 1.2, width = 0.1) +  # Jitter with controlled width
  scale_color_manual(values = custom_colors) +  # Apply colors to shape borders
  scale_fill_manual(values = custom_colors) +  # Apply colors to shape fills
  theme_classic() +  
  ggtitle("Fatty acid amide hydrolase 2", "Gene: 129976441, FAAH2 ") +
  theme(
    axis.title.x = element_blank(),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.text.x = element_text(size = 11, color = "black"),
    legend.position = "none"
  )
#####Miro1####
##Mitochondrial Rho GTPase 1
#a protein crucial for mitochondrial transport and homeostasis in cells
#It acts as an anchor on the mitochondrial outer membrane, facilitating the movement of mitochondria along microtubules via motor proteins like kinesin. 
#Miro1 plays a key role in mitochondrial transport. Essential for energy production, neuronal function, and overall cell survival. 

counts8 <- plotCounts(dds, "LOC129971588", 
                      intgroup = c("condition"), returnData = TRUE)

Miro1<-ggplot(counts8, aes(x = condition, y = count, 
                    color = condition, fill = condition, 
                    group = condition)) + 
  geom_boxplot(width = 0.1, alpha = 0.2, outlier.shape = NA) +  # Add boxplot with transparency
  geom_quasirandom(size = 2, alpha = 0.7, stroke = 1.2, width = 0.1) +  # Jitter with controlled width
  scale_color_manual(values = custom_colors) +  # Apply colors to shape borders
  scale_fill_manual(values = custom_colors) +  # Apply colors to shape fills
  theme_classic() +  
  ggtitle("Mitochondrial Rho GTPase 1", "Gene: 129971588, Miro1 ") +
  theme(
    axis.title.x = element_blank(),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.text.x = element_text(size = 11, color = "black"),
    legend.position = "none"
  )
Miro1
kruskal.test(counts8$count ~ counts8$condition)
#Kruskal-Wallis chi-squared = 9.8462, df = 2, p-value = 0.007277

dunnTest(counts8$count ~ counts8$condition, method="bh")

####LDLR ####
#Low-density lipoprotein receptor
# is a transmembrane protein that facilitates the uptake of cholesterol-rich low-density lipoproteins (LDL) into cells. This process, receptor-mediated endocytosis, 
#is crucial for regulating cholesterol levels and preventing the buildup of LDL

counts9 <- plotCounts(dds, "LOC129957447", 
                      intgroup = c("condition"), returnData = TRUE)

LDLR<- ggplot(counts9, aes(x = condition, y = count, 
                    color = condition, fill = condition, 
                    group = condition)) + 
  #geom_boxplot(width = 0.1, alpha = 0.2, outlier.shape = NA) +  # Add boxplot with transparency
  geom_quasirandom(size = 2, alpha = 0.7, stroke = 1.2, width = 0.1) +  # Jitter with controlled width
  scale_color_manual(values = custom_colors) +  # Apply colors to shape borders
  scale_fill_manual(values = custom_colors) +  # Apply colors to shape fills
  theme_classic() +  # Remove gray background
  ggtitle("low-density lipoprotein receptor", "LOC129957447, LDLR") +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    #axis.title = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.text.x = element_text(size = 11, color = "black"),
    legend.position = "none"
  )
LDLR

kruskal.test(counts9$count ~ counts9$condition)
#Kruskal-Wallis chi-squared = 8, df = 2, p-value = 0.01832

dunnTest(counts9$count ~ counts9$condition, method="bh")

#####P450 ########
#cytochrome P450 3A21-like
#They are key to metabolizing xenobiotics
#Low temperature exposure has been shown to impact cytochrome P450 expression
#shown to be upregulated in response to low temperature stress, 
#potentially aiding in the detoxification of compounds produced under stress conditions

counts10 <- plotCounts(dds, "LOC129961890", 
                      intgroup = c("condition"), returnData = TRUE)

P450<- ggplot(counts10, aes(x = condition, y = count, 
                     color = condition, fill = condition, 
                     group = condition)) + 
  #geom_boxplot(width = 0.1, alpha = 0.2, outlier.shape = NA) +  # Add boxplot with transparency
  geom_quasirandom(size = 2, alpha = 0.7, stroke = 1.2, width = 0.1) +  # Jitter with controlled width
  scale_color_manual(values = custom_colors) +  # Apply colors to shape borders
  scale_fill_manual(values = custom_colors) +  # Apply colors to shape fills
  theme_classic() +  
  ggtitle("cytochrome P450", "LOC129961890, P450") +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.text.x=element_blank(), 
    axis.ticks.x=element_blank(),
    legend.position = "none"
  )
P450
kruskal.test(counts10$count ~ counts10$condition)
#Kruskal-Wallis chi-squared = 7.7308, df = 2, p-value = 0.02095

dunnTest(counts10$count ~ counts10$condition, method="bh")

####DGK1#####
counts11 <- plotCounts(dds, "LOC129984733", 
                       intgroup = c("condition"), returnData = TRUE)

DGK1<-ggplot(counts11, aes(x = condition, y = count, 
                     color = condition, fill = condition, 
                     group = condition)) + 
  #geom_boxplot(width = 0.1, alpha = 0.2, outlier.shape = NA) +  # Add boxplot with transparency
  geom_quasirandom(size = 2, alpha = 0.7, stroke = 1.2, width = 0.1) +  # Jitter with controlled width
  scale_color_manual(values = custom_colors) +  # Apply colors to shape borders
  scale_fill_manual(values = custom_colors) +  # Apply colors to shape fills
  theme_classic() +  
  ggtitle("diacylglycerol kinase 1", "LOC129984733, DGK1") +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.text.x = element_text(size = 11, color = "black"),
    legend.position = "none"
  )
DGK1
kruskal.test(counts11$count ~ counts11$condition)
#Kruskal-Wallis chi-squared = 7.7308, df = 2, p-value = 0.02095

dunnTest(counts11$count ~ counts11$condition, method="bh")

###Hsp20####
counts12 <- plotCounts(dds, "LOC129988413", 
                       intgroup = c("condition"), returnData = TRUE)
ggplot(counts12, aes(x = condition, y = count, 
                     color = condition, fill = condition, 
                     group = condition)) + 
  geom_boxplot(width = 0.1, alpha = 0.2, outlier.shape = NA) +  # Add boxplot with transparency
  geom_quasirandom(size = 2, alpha = 0.7, stroke = 1.2, width = 0.1) +  # Jitter with controlled width
  scale_color_manual(values = custom_colors) +  # Apply colors to shape borders
  scale_fill_manual(values = custom_colors) +  # Apply colors to shape fills
  theme_classic() +  
  ggtitle("small heat shock protein", "LOC129988413, Hsp20 ") +
  theme(
    axis.title.x = element_blank(),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.text.x = element_text(size = 11, color = "black"),
    legend.position = "none"
  )

######CRYAA#####

counts13 <- plotCounts(dds, "LOC129984551", 
                       intgroup = c("condition"), returnData = TRUE)
ggplot(counts13, aes(x = condition, y = count, 
                     color = condition, fill = condition, 
                     group = condition)) + 
  geom_boxplot(width = 0.1, alpha = 0.2, outlier.shape = NA) +  # Add boxplot with transparency
  geom_quasirandom(size = 2, alpha = 0.7, stroke = 1.2, width = 0.1) +  # Jitter with controlled width
  scale_color_manual(values = custom_colors) +  # Apply colors to shape borders
  scale_fill_manual(values = custom_colors) +  # Apply colors to shape fills
  theme_classic() +  
  ggtitle("alpha-crystallin A", "LOC129984551, CRYAA ") +
  theme(
    axis.title.x = element_blank(),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.text.x = element_text(size = 11, color = "black"),
    legend.position = "none"
  )
#######AdSS ######

counts14 <- plotCounts(dds, "LOC129983880", 
                       intgroup = c("condition"), returnData = TRUE)
ggplot(counts14, aes(x = condition, y = count, 
                     color = condition, fill = condition, 
                     group = condition)) + 
  #geom_boxplot(width = 0.1, alpha = 0.2, outlier.shape = NA) +  # Add boxplot with transparency
  geom_quasirandom(size = 2, alpha = 0.7, stroke = 1.2, width = 0.1) +  # Jitter with controlled width
  scale_color_manual(values = custom_colors) +  # Apply colors to shape borders
  scale_fill_manual(values = custom_colors) +  # Apply colors to shape fills
  theme_classic() +  
  ggtitle("adenylosuccinate synthetase", "LOC129983880, AdSS") +
  theme(
    axis.title.x = element_blank(),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.text.x = element_text(size = 11, color = "black"),
    legend.position = "none"
  )

#for gene info use: e.g. res1["LOC129959427", ]
########
#LOC129963497

counts15 <- plotCounts(dds, "LOC129963497", 
                       intgroup = c("condition"), returnData = TRUE)
ggplot(counts15, aes(x = condition, y = count, 
                     color = condition, fill = condition, 
                     group = condition)) + 
  #geom_boxplot(width = 0.1, alpha = 0.2, outlier.shape = NA) +  # Add boxplot with transparency
  geom_quasirandom(size = 2, alpha = 0.7, stroke = 1.2, width = 0.1) +  # Jitter with controlled width
  scale_color_manual(values = custom_colors) +  # Apply colors to shape borders
  scale_fill_manual(values = custom_colors) +  # Apply colors to shape fills
  theme_classic() +  
  ggtitle("acyl-CoA dehydrogenase", "LOC129963497, ACAD") +
  theme(
    axis.title.x = element_blank(),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.text.x = element_text(size = 11, color = "black"),
    legend.position = "none"
  )

gridCold <- grid.arrange(apoLp, P450, LPAT,  ncol = 3, nrow = 1)
gridWarm <- grid.arrange(QORs, AACS, LDLR,  ncol = 3, nrow = 1)
gridall <-grid.arrange(apoLp, P450, LPAT, QORs, AACS, DGK1,  ncol = 3, nrow = 2)

###### volcano plots ########

library(EnhancedVolcano)
citation("EnhancedVolcano")

### cold vs moderate
#res1 <- res1[res1$baseMean > 10,]
#res1

res1.df <- as.data.frame(res1)

#first attempt
#EnhancedVolcano(res1.df, x = "log2FoldChange", y = "padj", lab = rownames(res1.df))

#second with enhancedvolcano 

#V1 <- EnhancedVolcano(res1.df, x = "log2FoldChange", y = "padj", lab = rownames(res1.df), pCutoff = 1e-4, FCcutoff = 1)
#selected = c("LOC129963668", "LOC129976206")
#EnhancedVolcano(res1.df, x = "log2FoldChange", y = "padj", lab = rownames(res1.df), pCutoff = 1e-4, FCcutoff = 1,
                #selectLab = selected)

# Define thresholds for significance
pvalueCutoff <- 0.05
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
  geom_text(aes(x = 2, y = -log10(pvalueCutoff) + 6, label = paste(num_upregulated, "")), vjust = 1, hjust = 1, color = "red") +
  geom_text(aes(x = -2, y = -log10(pvalueCutoff) + 6, label = paste(num_downregulated, "")), vjust = 1, hjust = 0, color = "blue") +
  theme_classic() +
  labs(x = "Log2 Fold Change", y = "-Log10 adjusted P-value") +  
  theme(axis.title = element_text(size = 14, face = "bold"),  # Bigger axis titles
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.text.x = element_text(size = 14, color = "black"),
        legend.position = "none")  # Remove the legend

volcano_plot1


tiff("VolcanoPlotgg_CvM.tiff", res = 250, width = 1300, height = 1500)
print(volcano_plot1)
dev.off()

EnhancedVolcano(res1.df, x = "log2FoldChange", y = "padj", lab = rownames(res1.df), pCutoff = 0.05, FCcutoff = 1)


########## In a volcano plot, the most upregulated genes are towards the right, the most downregulated genes are towards the left,
#and the most statistically significant genes are towards the top.

### warm vs moderate

#res2 <- res2[res2$baseMean > 50,]
#res2

res2.df <- as.data.frame(res2)
#EnhancedVolcano(res2.df, x = "log2FoldChange", y = "padj", lab = rownames(res2.df))


# Define thresholds for significance
pvalueCutoff <- 0.05  #1e-4
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
  geom_text(aes(x = 5, y = -log10(pvalueCutoff) + 5, label = paste(num_upregulated, "")), vjust = 1, hjust = 1, color = "red") +
  geom_text(aes(x = -5, y = -log10(pvalueCutoff) + 5, label = paste(num_downregulated, "")), vjust = 1, hjust = 0, color = "blue") +
  theme_classic() +
  labs(x = "Log2 Fold Change", y = "") +  # Removed title
  theme(axis.title = element_text(size = 14, face = "bold"),  # Bigger axis titles
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.text.x = element_text(size = 14, color = "black"),
        legend.position = "none")  # Remove the legend
volcano_plot2


tiff("VolcanoPlotgg_WvM.tiff", res = 250, width = 1300, height = 1500)
print(volcano_plot2)
dev.off()

EnhancedVolcano(res2.df, x = "log2FoldChange", y = "padj", lab = rownames(res2.df), pCutoff = 0.05, FCcutoff = 1)

#### Warm vs cold 

#res3 <- res3[res3$baseMean > 50,]

res3.df <- as.data.frame(res3)
#
#V3 <- EnhancedVolcano(res3.df, x = "log2FoldChange", y = "padj", lab = rownames(res3.df), pCutoff = 1e-4, FCcutoff = 1)

# Define thresholds for significance
pvalueCutoff <- 0.05
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
volcano_plot3 <- ggplot(res3.df, aes(x = log2FoldChange, y = -log10(padj), color = Color)) +
  geom_point(alpha = 0.5) +
  scale_color_identity() +
  geom_vline(xintercept = log2FCCutoff, linetype = "dashed", color = "black") +
  geom_vline(xintercept = -log2FCCutoff, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(pvalueCutoff), linetype = "dashed", color = "black") +
  geom_text(aes(x = 6, y = -log10(pvalueCutoff) + 10, label = paste(num_upregulated, "")), vjust = 1, hjust = 1, color = "red") +
  geom_text(aes(x = -6, y = -log10(pvalueCutoff) + 10, label = paste(num_downregulated, "")), vjust = 1, hjust = 0, color = "blue") +
  theme_classic() +
  labs(x = "Log2 Fold Change", y = "") +  # Removed title
  theme(axis.title = element_text(size = 14, face = "bold"),  # Bigger axis titles
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.text.x = element_text(size = 14, color = "black"),
        legend.position = "none")  # Remove the legend  # Remove the legend
volcano_plot3

tiff("VolcanoPlotgg_WvC.tiff", res = 250, width = 1300, height = 1500)
print(volcano_plot3)
dev.off()

EnhancedVolcano(res3.df, x = "log2FoldChange", y = "padj", lab = rownames(res3.df))

Volcano_grid <- grid.arrange(volcano_plot1, volcano_plot2, volcano_plot, 
                           ncol = 3, nrow = 1,
                           widths = c(1, 1, 1),  # Equal column widths
                           #heights = c(1, 1,1),     # Equal row heights
                           padding = unit(0.5, "line"))  # Reduce space between plots
########## 
#####ontology of genes (GO) (see the GO script)
