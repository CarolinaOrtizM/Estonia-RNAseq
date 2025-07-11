#Go analysis  

# Load required packages
library(dplyr)
library(topGO)
library(GOplot)
library(BiocManager)

# Load the GAF file
gaf_file <- "GCF_947563725.1-RS_2023_05_gene_ontology.gaf"
gaf_data <- read.table(gaf_file, header = FALSE, sep = "\t", comment.char = "!", quote = "")

# Select necessary columns (Gene ID and GO term)
gaf_annotations <- dplyr::select(gaf_data, V3, V5)

# Convert the GAF data into a format suitable for topGO
geneID2GO <- by(gaf_annotations$V5, gaf_annotations$V3, function(x) as.character(x))

# Prepare results data frames from your differential expression analysis
cold_moderate_results <- setNames(cbind(rownames(df_sigs1), df_sigs1, row.names = NULL),
                                  c("Gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"))
warm_moderate_results <- setNames(cbind(rownames(df_sigs2), df_sigs2, row.names = NULL),
                                  c("Gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"))
warm_cold_results <- setNames(cbind(rownames(df_sigs3), df_sigs3, row.names = NULL),
                              c("Gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"))

###GO analysis ####
# Combine all genes into a universe
gene_universe <- unique(c(cold_moderate_results$Gene_id, warm_moderate_results$Gene_id, warm_cold_results$Gene_id))

# Function to create gene lists indicating differentially expressed genes
create_geneList <- function(deg_results, universe) {
  genePresence <- universe %in% deg_results$Gene_id
  geneList <- factor(genePresence, levels = c(FALSE, TRUE), labels = c("0", "1"))
  names(geneList) <- universe
  return(geneList)
}

# Create gene lists for each condition
coldMod_genes <- create_geneList(cold_moderate_results, gene_universe)
WarMod_genes <- create_geneList(warm_moderate_results, gene_universe)
warmCold_genes <- create_geneList(warm_cold_results, gene_universe)

# Function to prepare topGOdata for each condition and ontology
prepare_GOdata <- function(geneList, ontology) {
  GOdata <- new("topGOdata", 
                description = "GO analysis",
                ontology = ontology,
                allGenes = geneList,
                geneSel = function(p) p == "1",
                nodeSize = 10,
                annot = annFUN.gene2GO,
                gene2GO = geneID2GO)
  return(GOdata)
}

# Modified run_topGO function to include gene information
run_topGO <- function(geneList, ontology, deg_results) {
  GOdata <- prepare_GOdata(geneList, ontology)
  result <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  
  # Get significant genes for each GO term
  allGO <- usedGO(GOdata)
  gene_GO_mapping <- genesInTerm(GOdata, allGO)
  
  # Create a data frame with GO terms and associated genes
  go_genes <- lapply(names(gene_GO_mapping), function(go) {
    genes_in_go <- gene_GO_mapping[[go]]
    # Get DE information for these genes
    de_genes <- deg_results[deg_results$Gene_id %in% genes_in_go, ]
    up_genes <- de_genes[de_genes$log2FoldChange > 0, "Gene_id"]
    down_genes <- de_genes[de_genes$log2FoldChange < 0, "Gene_id"]
    
    data.frame(
      GO.ID = go,
      Upregulated_Genes = ifelse(length(up_genes) > 0, paste(up_genes, collapse = ","), NA),
      Downregulated_Genes = ifelse(length(down_genes) > 0, paste(down_genes, collapse = ","), NA),
      Num_Upregulated = length(up_genes),
      Num_Downregulated = length(down_genes),
      stringsAsFactors = FALSE
    )
  })
  
  go_genes_df <- do.call(rbind, go_genes)
  
  # Get the standard results table
  allRes <- GenTable(GOdata, classicFisher = result, orderBy = "classicFisher", topNodes = length(allGO))
  
  # Merge with gene information
  merged_results <- merge(allRes, go_genes_df, by = "GO.ID")
  
  return(merged_results)
}

# Modified process_all_ontologies function
process_all_ontologies <- function(geneList, deg_results) {
  categories <- c("BP", "MF", "CC")
  results_list <- list()
  
  for (category in categories) {
    results <- run_topGO(geneList, category, deg_results)
    results_list[[category]] <- results
  }
  
  return(results_list)
}

# Running and storing results for each condition with DEG information
results_coldMod <- process_all_ontologies(coldMod_genes, cold_moderate_results)
results_warmMod <- process_all_ontologies(WarMod_genes, warm_moderate_results)
results_warmCold <- process_all_ontologies(warmCold_genes, warm_cold_results)

# Function to prepare data frame with additional metadata columns
prepare_data_frame <- function(df, category, condition) {
  df$Category <- category
  df$Condition <- condition
  return(df)
}

# Prepare individual data frames
results_coldMod_BP <- prepare_data_frame(results_coldMod$BP, "BP", "Cold")
results_coldMod_MF <- prepare_data_frame(results_coldMod$MF, "MF", "Cold")
results_coldMod_CC <- prepare_data_frame(results_coldMod$CC, "CC", "Cold")

results_warmMod_BP <- prepare_data_frame(results_warmMod$BP, "BP", "Warm")
results_WarmMod_MF <- prepare_data_frame(results_warmMod$MF, "MF", "Warm")
results_WarmMod_CC <- prepare_data_frame(results_warmMod$CC, "CC", "Warm")

results_warmCold_BP <- prepare_data_frame(results_warmCold$BP, "BP", "WarmCold")
results_warmCold_MF <- prepare_data_frame(results_warmCold$MF, "MF", "WarmCold")
results_warmCold_CC <- prepare_data_frame(results_warmCold$CC, "CC", "WarmCold")

# Combine all data frames into a single data frame
all_results <- rbind(
  results_coldMod_BP, results_coldMod_MF, results_coldMod_CC,
  results_warmMod_BP, results_WarmMod_MF, results_WarmMod_CC,
  results_warmCold_BP, results_warmCold_MF, results_warmCold_CC
)

# Optionally, sort or order the combined data frame as needed
all_results_sorted <- all_results[order(all_results$Condition, all_results$Category),]

# View or export the combined data frame
print(head(all_results_sorted))
write.csv(all_results_sorted, "GO_analysis_combined_results_with_genes_all.csv", row.names = FALSE)

#### chord plots ####

##1. warm vs cold 

results_WarmCold <-read.csv("GO_WarmCold.csv", header = T) 
head(results_WarmCold)
results_WarmCold$Genes <-as.character(results_WarmCold$Genes)
str(results_WarmCold)

results_WarmCold <- results_WarmCold[c(-2, -3, -6),]
#results_WarmCold <- results_WarmCold[-3,]

df_sigs3 <- as_tibble(df_sigs3, rownames = "ID")

# We need to keep only the relevant columns and rename them
gene_data <- df_sigs3[, c("ID", "log2FoldChange", "padj")]
colnames(gene_data) <- c("ID", "logFC", "adj.P.Val")

gene_data_df<-as.data.frame(gene_data)

# We need to process the Genes column to be a character vector of gene IDs
enrichment_data <- results_WarmCold
#enrichment_data$Genes <- gsub(" ", "", enrichment_data$Genes) # Remove any spaces

# Create the circle data
circ <- circle_dat(enrichment_data, gene_data_df)

# Select the top processes to display (adjust number as needed)
top_processes <- 6  # You can change this number
selected_terms <- unique(circ$term)[1:top_processes]

# Subset the data for these terms
circ_sub <- subset(circ, term %in% selected_terms)

# Create the chord plot
chord <- chord_dat(data = circ_sub, genes = gene_data_df, process = selected_terms)

# Generate the plot
GOChord(chord, 
        space = 0.02,           # Space between gene segments
        gene.order = 'logFC',   # Order genes by logFC
        gene.space = 0.25,      # Space between genes
        gene.size = 3,          # Gene label size
        border.size = 0.1,      # Border thickness
        process.label = 5)      # Process label size
"LOC129964053" %in% circ_sub$genes

###### Warm vs moderate 

results_WarmMod <-read.csv("GO_WarmMode.csv", header = T) 
head(results_WarmMod)
results_WarmMod$Genes <-as.character(results_WarmMod$Genes)
str(results_WarmMod)

results_WarmMod <- results_WarmMod[-3,]

df_sigs2 <- as_tibble(df_sigs2, rownames = "ID")

# We need to keep only the relevant columns and rename them
gene_data <- df_sigs2[, c("ID", "log2FoldChange", "padj")]
colnames(gene_data) <- c("ID", "logFC", "adj.P.Val")

gene_data_df<-as.data.frame(gene_data)

# We need to process the Genes column to be a character vector of gene IDs
enrichment_data <- results_WarmMod
#enrichment_data$Genes <- gsub(" ", "", enrichment_data$Genes) # Remove any spaces

# Create the circle data
circ <- circle_dat(enrichment_data, gene_data_df)

# Select the top processes to display (adjust number as needed)
top_processes <- 7  # You can change this number
selected_terms <- unique(circ$term)[1:top_processes]

# Subset the data for these terms
circ_sub <- subset(circ, term %in% selected_terms)

# Create the chord plot
chord <- chord_dat(data = circ_sub, genes = gene_data_df, process = selected_terms)

# Generate the plot
GOChord(chord, 
        space = 0.02,           # Space between gene segments
        gene.order = 'logFC',   # Order genes by logFC
        gene.space = 0.25,      # Space between genes
        gene.size = 3,          # Gene label size
        border.size = 0.1,      # Border thickness
        process.label = 5)      # Process label size


# Get all unique genes in the plot
genes_in_plot <- unique(circ_sub$genes)
print(genes_in_plot)
# Check if your gene is in the plot
"LOC129964270" %in% circ_sub$genes

# View all GO term connections for your gene
circ_sub[circ_sub$genes == "LOC129987920", ]

###Cold vs moderate

results_ColdMod <-read.csv("GO_ColdMode.csv", header = T) 
head(results_ColdMod)
results_ColdMod$Genes <-as.character(results_ColdMod$Genes)
str(results_ColdMod)


df_sigs1 <- as_tibble(df_sigs1, rownames = "ID")

# We need to keep only the relevant columns and rename them
gene_data <- df_sigs1[, c("ID", "log2FoldChange", "padj")]
colnames(gene_data) <- c("ID", "logFC", "adj.P.Val")

gene_data_df<-as.data.frame(gene_data)

# We need to process the Genes column to be a character vector of gene IDs
enrichment_data <- results_ColdMod
#enrichment_data$Genes <- gsub(" ", "", enrichment_data$Genes) # Remove any spaces

# Create the circle data
circ <- circle_dat(enrichment_data, gene_data_df)

# Select the top processes to display (adjust number as needed)
top_processes <- 11  # You can change this number
selected_terms <- unique(circ$term)[1:top_processes]

# Subset the data for these terms
circ_sub <- subset(circ, term %in% selected_terms)

# Create the chord plot
chord <- chord_dat(data = circ_sub, genes = gene_data_df, process = selected_terms)

# Generate the plot
GOChord(chord, 
        space = 0.02,           # Space between gene segments
        gene.order = 'logFC',   # Order genes by logFC
        gene.space = 0.25,      # Space between genes
        gene.size = 3,          # Gene label size
        border.size = 0.1,      # Border thickness
        process.label = 5)      # Process label size


"LOC129959427" %in% circ_sub$genes


