#Go analysis  

# Load required packages
library(dplyr)
library(topGO)
library(BiocManager)

# Load the GAF file
gaf_file <- "GCF_947563725.1-RS_2023_05_gene_ontology.gaf"
gaf_data <- read.table(gaf_file, header = FALSE, sep = "\t", comment.char = "!", quote = "")

# Select necessary columns (Gene ID and GO term)
gaf_annotations <- dplyr::select(gaf_data, V3, V5)

# Convert the GAF data into a format suitable for topGO
geneID2GO <- by(gaf_annotations$V5, gaf_annotations$V3, function(x) as.character(x))

# Prepare results data frames from your differential expression analysis (see DESeq_analysis.R)
cold_results <- setNames(cbind(rownames(df_sigs1), df_sigs1, row.names = NULL),
                         c("Gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"))
moderate_results <- setNames(cbind(rownames(df_sigs2), df_sigs2, row.names = NULL),
                             c("Gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"))
warm_results <- setNames(cbind(rownames(df_sigs3), df_sigs3, row.names = NULL),
                         c("Gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"))

# Combine all genes into a universe
gene_universe <- unique(c(cold_results$Gene_id, moderate_results$Gene_id, warm_results$Gene_id))

# Function to create gene lists indicating differentially expressed genes
create_geneList <- function(deg_results, universe) {
  genePresence <- universe %in% deg_results$Gene_id
  geneList <- factor(genePresence, levels = c(FALSE, TRUE), labels = c("0", "1"))
  names(geneList) <- universe
  return(geneList)
}

# Create gene lists for each condition
cold_genes <- create_geneList(cold_results, gene_universe)
moderate_genes <- create_geneList(moderate_results, gene_universe)
warm_genes <- create_geneList(warm_results, gene_universe)

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

# Function to run the topGO analysis for a given ontology
run_topGO <- function(geneList, ontology) {
  GOdata <- prepare_GOdata(geneList, ontology)
  result <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  allRes <- GenTable(GOdata, classicFisher = result, orderBy = "classicFisher", topNodes = 10)
  return(allRes)
}

# Function to process all GO categories for a given gene list
process_all_ontologies <- function(geneList) {
  categories <- c("BP", "MF", "CC")
  results_list <- list()
  
  for (category in categories) {
    results <- run_topGO(geneList, category)
    results_list[[category]] <- results
  }
  
  return(results_list)
}

# Running and storing results for each condition
results_cold <- process_all_ontologies(cold_genes)
results_moderate <- process_all_ontologies(moderate_genes)
results_warm <- process_all_ontologies(warm_genes)

# Function to prepare data frame with additional metadata columns
prepare_data_frame <- function(df, category, condition) {
  df$Category <- category
  df$Condition <- condition
  return(df)
}

# Prepare individual data frames
results_cold_BP <- prepare_data_frame(results_cold$BP, "BP", "Cold")
results_cold_MF <- prepare_data_frame(results_cold$MF, "MF", "Cold")
results_cold_CC <- prepare_data_frame(results_cold$CC, "CC", "Cold")

results_moderate_BP <- prepare_data_frame(results_moderate$BP, "BP", "Moderate")
results_moderate_MF <- prepare_data_frame(results_moderate$MF, "MF", "Moderate")
results_moderate_CC <- prepare_data_frame(results_moderate$CC, "CC", "Moderate")

results_warm_BP <- prepare_data_frame(results_warm$BP, "BP", "Warm")
results_warm_MF <- prepare_data_frame(results_warm$MF, "MF", "Warm")
results_warm_CC <- prepare_data_frame(results_warm$CC, "CC", "Warm")

# Combine all data frames into a single data frame
all_results <- rbind(
  results_cold_BP, results_cold_MF, results_cold_CC,
  results_moderate_BP, results_moderate_MF, results_moderate_CC,
  results_warm_BP, results_warm_MF, results_warm_CC
)

# Optionally, sort or order the combined data frame as needed
all_results_sorted <- all_results[order(all_results$Condition, all_results$Category),]

# View or export the combined data frame
print(head(all_results_sorted))
write.csv(all_results_sorted, "GO_analysis_combined_results.csv", row.names = FALSE)
