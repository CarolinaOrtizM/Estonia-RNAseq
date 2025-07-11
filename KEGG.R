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

#KEGG analysis 
##########Cold vs moderate DEGs (df_sigs1)#############

rm(list=ls(all=TRUE)); # graphics.off()

####data (DGE)
# Function to strip 'LOC' prefix from IDs
strip_loc_prefix <- function(ids) {
  return(sub("LOC", "", ids))
}

df_sigs1 <- as_tibble(df_sigs1, rownames = "ID")

# Strip 'LOC' prefix from gene IDs in the results
df_sigs1$ID <- strip_loc_prefix(df_sigs1$ID )


gene_ids1 <- unique(unlist(df_sigs1))
head(gene_ids1)

# Perform KEGG enrichment analysis

kk1 <- enrichKEGG(gene = gene_ids1,
                   organism = 'abru',
                   keyType = "kegg",
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.05)

# Check the structure of kk again to see if there are any anomalies
str(kk1@result)

# Convert relevant columns to numeric 
kk1@result$GeneRatioCount <- as.numeric(sapply(strsplit(kk1@result$GeneRatio, "/"), `[`, 1))
kk1@result$BgRatioCount <- as.numeric(sapply(strsplit(kk1@result$BgRatio, "/"), `[`, 1))

# only kk@result as a data frame 
data1 <- as.data.frame(kk1@result)
write.csv(data1, file =  "KKresuktsCold_mode.csv")
  
####Plotting
filtered_data1 <- data1 %>%
filter(pvalue < 0.05)  # Adjust this value based on your desired threshold

plot_data1 <- filtered_data1 

# Remove the repetitive part from the Description
plot_data1 <- plot_data1 %>%
  mutate(Description = gsub(" - Argiope bruennichi \\(wasp spider\\)", "", Description))

### Dot plot #########
ggplot(plot_data1, aes(x = reorder(Description, pvalue), y = -log10(pvalue))) +
  geom_point(aes(size = Count, color = pvalue), alpha = 0.6) +
  labs(x = "Pathway Description", y = "-log10(P-value)", title = "KEGG Pathway Enrichment Analysis Cold",
       subtitle = "Size of points corresponds to number of genes in the pathway") +
  scale_color_gradient(low = "green", high = "red") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12)) +
  guides(size = guide_legend(title = "Gene Count"))

###Bar plot #####
ggplot(plot_data1, aes(x = reorder(Description, Count), y = Count, fill = category)) +
  geom_bar(stat = "identity", width = 0.4, color = "black",size = 0.5)+
  coord_flip() +  # Horizontal bars
  labs(x = NULL, 
       y = "Number of DEGs in each pathway", 
       title = "") +
  geom_text(aes(label = Count), hjust = -0.2, size = 5) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 12, color = "black",face = "bold"),
        axis.text.x = element_text(size = 14, color = "black"),
        axis.line.y = element_blank(),  # Remove Y-axis line
        axis.ticks.y = element_blank(),  # Remove Y-axis ticks (optional)
        legend.position = "right",
        legend.title = element_blank(),
        plot.margin = margin(5, 10, 5, 5)) +
  scale_fill_manual(values = c("Metabolism" = "#1A5555", 
                               "Genetic Information Processing" = "green",
                               "Cellular Processes" = "#00FFFF")) +
  scale_y_continuous(
    limits = c(0, 30),        # Hard stop at 0 and 6
    expand = c(0, 0),        # No padding
    breaks = seq(0, 30, by = 5)
  )

########Warm vs moderate DEGs (df_sigs2)############
####data (DGE)
# Function to strip 'LOC' prefix from IDs
strip_loc_prefix <- function(ids) {
  return(sub("LOC", "", ids))
}

df_sigs2 <- as_tibble(df_sigs2, rownames = "ID")

# Strip 'LOC' prefix from gene IDs in the results
df_sigs2$ID <- strip_loc_prefix(df_sigs2$ID )

#gene list 
gene_ids2 <- unique(unlist(df_sigs2))
head(gene_ids2)

# Perform KEGG enrichment analysis

kk2 <- enrichKEGG(gene = gene_ids2,
                  organism = 'abru',
                  keyType = "kegg",
                  pAdjustMethod = "BH",
                  qvalueCutoff = 0.05)

# Check the structure of kk again to see if there are any anomalies
str(kk2@result)

# Convert relevant columns to numeric 
kk2@result$GeneRatioCount <- as.numeric(sapply(strsplit(kk2@result$GeneRatio, "/"), `[`, 1))
kk2@result$BgRatioCount <- as.numeric(sapply(strsplit(kk2@result$BgRatio, "/"), `[`, 1))

# only kk@result as a data frame 
data2 <- as.data.frame(kk2@result)
write.csv(data2, file =  "KKresuktsWarm_mode.csv")


####Plotting
filtered_data2 <- data2 %>%
  filter(pvalue < 0.05)  # Adjust this value based on your desired threshold

plot_data2 <- filtered_data2 

# Remove the repetitive part from the Description
plot_data2 <- plot_data2 %>%
  mutate(Description = gsub(" - Argiope bruennichi \\(wasp spider\\)", "", Description))

### Dot plot #########
ggplot(plot_data2, aes(x = reorder(Description, pvalue), y = -log10(pvalue))) +
  geom_point(aes(size = Count, color = pvalue), alpha = 0.6) +
  labs(x = "Pathway Description", y = "-log10(P-value)", title = "KEGG Pathway Enrichment Analysis Warm",
       subtitle = "Size of points corresponds to number of genes in the pathway") +
  scale_color_gradient(low = "green", high = "red") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12)) +
  guides(size = guide_legend(title = "Gene Count"))


###Bar plot #####
ggplot(plot_data2, aes(x = reorder(Description, Count), y = Count, fill = category)) +
  geom_bar(stat = "identity", width = 0.4, color = "black",size = 0.5)+
  coord_flip() +  # Horizontal bars
  labs(x = NULL, 
       y = "Number of DEGs in each pathway", 
       title = "") +
  geom_text(aes(label = Count), hjust = -0.2, size = 5) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 12, color = "black",face = "bold"),
        axis.text.x = element_text(size = 14, color = "black"),
        axis.line.y = element_blank(),  # Remove Y-axis line
        axis.ticks.y = element_blank(),  # Remove Y-axis ticks (optional)
        legend.position = "right",
        legend.title = element_blank(),
        plot.margin = margin(5, 10, 5, 5)) +
  scale_fill_manual(values = c("Metabolism" = "#1A5555", 
                               "Genetic Information Processing" = "green",
                               "Cellular Processes" = "#00FFFF")) +
  scale_y_continuous(
    limits = c(0, 30),        # Hard stop at 0 and 6
    expand = c(0, 0),        # No padding
    breaks = seq(0, 30, by = 5)            # Explicit breaks at every integer
  )

#######Warm vs Cold DEGs (df_sigs3)############

####data (DGE)
# Function to strip 'LOC' prefix from IDs
strip_loc_prefix <- function(ids) {
  return(sub("LOC", "", ids))
}

df_sigs3 <- as_tibble(df_sigs3, rownames = "ID")

# Strip 'LOC' prefix from gene IDs in the results
df_sigs3$ID <- strip_loc_prefix(df_sigs3$ID )

#gene list 
gene_ids3 <- unique(unlist(df_sigs3))
head(gene_ids3)

# Perform KEGG enrichment analysis

kk3 <- enrichKEGG(gene = gene_ids3,
                  organism = 'abru',
                  keyType = "kegg",
                  pvalueCutoff = 0.05,
                  pAdjustMethod = "BH",
                  qvalueCutoff = 0.1)

# Check the structure of kk again to see if there are any anomalies
str(kk3@result)

# Convert relevant columns to numeric 
kk3@result$GeneRatioCount <- as.numeric(sapply(strsplit(kk3@result$GeneRatio, "/"), `[`, 1))
kk3@result$BgRatioCount <- as.numeric(sapply(strsplit(kk3@result$BgRatio, "/"), `[`, 1))

# only kk@result as a data frame 
data3 <- as.data.frame(kk3@result)

write.csv(data3, file =  "KKresuktsWarm_Cold.csv")

#data3 <- read.csv("KKresuktsWarm_Cold.csv", header = T)

####Plotting
filtered_data3 <- data3 %>%
  filter(pvalue < 0.05)  # Adjust this value based on your desired threshold

plot_data3 <- filtered_data3 

# Remove the repetitive part from the Description
plot_data3 <- plot_data3 %>%
  mutate(Description = gsub(" - Argiope bruennichi \\(wasp spider\\)", "", Description))

#filtered by 
top_data3 <- data3 %>% 
  arrange(pvalue) %>%  # Sort by pvalue to get the most significant results
  head(20)  # Select the top 20

plot_data33 <- top_data3 

# Remove the repetitive part from the Description
plot_data33 <- plot_data33 %>%
  mutate(Description = gsub(" - Argiope bruennichi \\(wasp spider\\)", "", Description))

### Dot plot #########
ggplot(plot_data33, aes(x = reorder(Description, pvalue), y = -log10(pvalue))) +
  geom_point(aes(size = Count, color = pvalue), alpha = 0.6) +
  labs(x = "Pathway Description", y = "-log10(P-value)", title = "KEGG Pathway Enrichment Analysis WarmvsCold",
       subtitle = "Size of points corresponds to number of genes in the pathway") +
  scale_color_gradient(low = "green", high = "red") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12)) +
  guides(size = guide_legend(title = "Gene Count"))


###Bar plot #####
ggplot(plot_data3, aes(x = reorder(Description, Count), y = Count, fill = category)) +
  geom_bar(stat = "identity", width = 0.4, color = "black",size = 0.5)+
  coord_flip() +  # Horizontal bars
  labs(x = NULL, 
       y = "Number of DEGs in each pathway", 
       title = "") +
  geom_text(aes(label = Count), hjust = -0.2, size = 5) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 12, color = "black",face = "bold"),
        axis.text.x = element_text(size = 14, color = "black"),
        axis.line.y = element_blank(),  # Remove Y-axis line
        axis.ticks.y = element_blank(),  # Remove Y-axis ticks (optional)
        legend.position = "right",
        legend.title = element_blank(),
        plot.margin = margin(5, 10, 5, 5)) +
  scale_fill_manual(values = c("Metabolism" = "#1A5555", 
                               "Genetic Information Processing" = "green",
                               "Cellular Processes" = "#00FFFF",
                               "Environmental Information Processing" = "purple",
                               "Signaling" ="darkorange")) +
  scale_y_continuous(
    limits = c(0, 30),        # Hard stop at 0 and 6
    expand = c(0, 0),        # No padding
    breaks = seq(0, 30, by = 5)              # Explicit breaks at every integer
  )


# Update row 5 specifically
plot_data3[5, "category"] <- "Cellular Processes"
plot_data3[5, "subcategory"] <- "Cell motility"

# Verify the change
plot_data33[20,"category"] <- "Cellular Processes"
plot_data33[20, "subcategory"] <- "Transport and catabolism"

