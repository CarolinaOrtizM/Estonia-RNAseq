#####KEGG enrichment analysis ###########
#needs to be run after the GO one and the DESEq_analysis. 

#libraries  
library(KEGGREST)
library(clusterProfiler)
library(dplyr)
library(tidyr)
library(ggplot2)
library(biomaRt)
library(enrichplot)
library(ComplexHeatmap)
library(circlize) 

####data (DGE)
# Function to strip 'LOC' prefix from IDs
strip_loc_prefix <- function(ids) {
  return(sub("LOC", "", ids))
}

# Strip 'LOC' prefix from gene IDs in the results
cold_results$Gene_id <- strip_loc_prefix(cold_results$Gene_id)
moderate_results$Gene_id <- strip_loc_prefix(moderate_results$Gene_id)
warm_results$Gene_id <- strip_loc_prefix(warm_results$Gene_id)

# Combine all genes into a universe
gene_universe <- unique(c(cold_results, moderate_results, warm_results))

# Verify the gene universe
print(head(gene_universe))
str(gene_universe)

# Extract gene IDs from the specified sublist positions
gene_ids_list <- lapply(gene_universe[c(1, 8, 15)], unlist)

# Flatten the list and make it unique to get rid of duplicates
gene_ids <- unique(unlist(gene_ids_list))
head(gene_ids)

# Perform KEGG enrichment analysis
tryCatch({
  kk <- enrichKEGG(gene = gene_ids,
                   organism = 'abru',
                   keyType = "kegg",
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.05)
}, error = function(e) {
  cat("Error in enrichKEGG: ", e$message, "\n")
})

# Check the structure of kk to see if there are any anomalies
str(kk@result)

# Convert relevant columns to numeric where necessary
kk@result$GeneRatioCount <- as.numeric(sapply(strsplit(kk@result$GeneRatio, "/"), `[`, 1))
kk@result$BgRatioCount <- as.numeric(sapply(strsplit(kk@result$BgRatio, "/"), `[`, 1))

# only kk@result as a data frame 
data <- as.data.frame(kk@result)
write.csv(data, file =  "KKresukts.csv")

##### Option 1: Selecting the top 40 most significant pathways
top_data <- data %>% 
  arrange(p.adjust) %>%  # Sort by p.adjust to get the most significant results
  head(40)  # Select the top 40

##### Option 2: Filtering by a p-adjust threshold
#filtered_data <- data %>%
  #filter(p.adjust < 0.05)  # Adjust this value based on your desired threshold

# Proceed with plotting using either 'top_data' or 'filtered_data'
plot_data <- top_data  # or filtered_data
#plot_data <- filtered_data #doesn't work...


# Remove the repetitive part from the Description
plot_data <- plot_data %>%
  mutate(Description = gsub(" - Argiope bruennichi \\(wasp spider\\)", "", Description))

# Recreate the plot with updated descriptions
ggplot(plot_data, aes(x = reorder(Description, p.adjust), y = -log10(pvalue))) +
  geom_point(aes(size = Count, color = p.adjust), alpha = 0.6) +
  labs(x = "Pathway Description", y = "-log10(P-value)", title = "Top 40 KEGG Pathway Enrichment Analysis",
       subtitle = "Size of points corresponds to number of genes in the pathway") +
  scale_color_gradient(low = "green", high = "red") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12)) +
  guides(size = guide_legend(title = "Gene Count"))
ggsave("KEGG_all_top40.tiff", width=8, height=10) ##for saving it 

########## for the DGE data ### to know which treatment accounts for which pathway

# Example structure, adjust according to your actual data structure and names
cold_genes <- data.frame(geneID = gene_universe[[1]], logFC = gene_universe[[3]], pvalue = gene_universe[[6]], condition = "cold")
moderate_genes <- data.frame(geneID = gene_universe[[8]], logFC = gene_universe[[10]], pvalue = gene_universe[[14]], condition = "moderate")
warm_genes <- data.frame(geneID = gene_universe[[15]], logFC = gene_universe[[17]], pvalue = gene_universe[[20]], condition = "warm")

# Combine into one data frame
combined_genes <- rbind(cold_genes, moderate_genes, warm_genes)

# 'kk' is the KEGG enrichment result and 'kk@result$geneID' contains gene IDs separated by '/'
enriched_genes <- strsplit(as.character(kk@result$geneID), "/")

# Map each list of genes to the corresponding pathway description
pathway_genes <- setNames(enriched_genes, kk@result$Description)

# Function to find conditions associated with each gene in a pathway
pathway_conditions <- lapply(pathway_genes, function(gene_list) {
  gene_data <- combined_genes[combined_genes$geneID %in% gene_list, ]
  return(gene_data)
})

# Optionally, create a summary of conditions for each pathway
pathway_summary <- lapply(pathway_conditions, function(df) {
  if (nrow(df) > 0) {
    return(table(df$condition))
  } else {
    return(NULL)
  }
})

head(pathway_summary)
str(pathway_summary)

write.csv(final_df, file ="Pathways_conditions_numbers.csv") 

# Save each pathway's gene condition data to a separate CSV file
lapply(names(pathway_conditions), function(pathway) {
  if (!is.null(pathway_conditions[[pathway]])) {
    write.csv(pathway_conditions[[pathway]], sprintf("%s.csv", gsub("/", "_", pathway)), row.names = FALSE)
  }
})

# Add a pathway column and combine all into a single dataframe
all_pathways_df <- do.call(rbind, lapply(names(pathway_conditions), function(pathway) {
  if (!is.null(pathway_conditions[[pathway]])) {
    df <- pathway_conditions[[pathway]]
    df$pathway <- pathway  # Add the pathway name as a new column
    return(df)
  }
}))

# Save the combined dataframe to a CSV file
write.csv(all_pathways_df, file = "AllPath.csv")

# Remove 'pathway' and 'category' columns
all_pathways_df$pathway <- NULL

# Add 'category' and "pathway" column from Pathways_cate
all_pathways_df$pathway <- Pathways_cate$pathway
all_pathways_df$category <- Pathways_cate$category

# Clean row names directly in the matrix
rownames(all_pathways_df) <- gsub(" - Argiope bruennichi \\(wasp spider\\)", "", rownames(all_pathways_df))

# Print cleaned row names for verification
cat("Cleaned Row Names:\n")
print(rownames(heatmap_data))

###to plot
# Libraries loading
library(ComplexHeatmap)
library(dplyr)
library(tidyr)
library(tibble)  # For column_to_rownames
library(circlize)  # For colorRamp2

# Aggregate data by pathway and condition
pathway_summary <- all_pathways_df %>%
  group_by(pathway, condition) %>%
  summarise(mean_logFC = mean(logFC, na.rm = TRUE), .groups = 'drop')

# Prepare the heatmap data
heatmap_data <- pathway_summary %>%
  pivot_wider(names_from = condition, values_from = mean_logFC) %>%
  column_to_rownames("pathway") %>%
  as.matrix()

# Calculate variance and select top 40 pathways
variances <- apply(heatmap_data, 1, var)
top_pathways <- names(sort(variances, decreasing = TRUE))[1:40]
heatmap_data_top <- heatmap_data[top_pathways, ]

# Define color map and annotations
color_map <- colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red"))
categories <- sample(c("Metabolism", "Cellular Processes", "Signaling", "Processing"), size = nrow(heatmap_data_top), replace = TRUE)
category_colors <- c("Metabolism" = "blue", "Cellular Processes" = "pink", "Signaling" = "yellow", "Processing" ="purple")
row_annotation <- HeatmapAnnotation(category = factor(categories),
                                    col = list(category = category_colors),
                                    which = "row",
                                    show_legend = TRUE)

# Generate and draw the heatmap
heatmap_top40 <- Heatmap(heatmap_data_top,
                         name = "log2FC",
                         row_title = "",
                         column_title = "",
                         right_annotation = row_annotation,
                         cluster_rows = TRUE,
                         cluster_columns = TRUE,
                         show_row_names = TRUE,
                         show_column_names = TRUE,
                         clustering_distance_rows = "euclidean",
                         clustering_distance_columns = "euclidean",
                         clustering_method_rows = "complete",
                         clustering_method_columns = "complete",
                         row_dend_reorder = TRUE,
                         column_dend_reorder = TRUE,
                         col = color_map,
                         show_heatmap_legend = TRUE,
                         row_names_gp = gpar(fontsize = 10), # Adjust font size as needed
                         heatmap_width = unit(20, "cm"), # Adjust width as necessary
                         heatmap_height = unit(20, "cm"), # Adjust height as necessary
                         ) 

draw(heatmap_top40)

draw(heatmap_top40, heatmap_legend_side = "left", annotation_legend_side = "bottom")

