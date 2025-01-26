# RNAseq analysis of cold-adapted _Argiope bruennichi_ spiderlings from the edge of the distribution in Estonia subjected to three overwintering regimes.

## List of files 

1. Bioinformatic_bash.md -- script for sequence data processing from raw reads to count table. 
2. CountTable.csv -- gene counts output from STAR and FeatureCounts (Subread).
3. SampleFile.txt -- information about the samples, including the replicates and samples submitted to NCBI.
5. Deseq_analysis.R -- Script for the differential expression analysis
6. GO_analysis.R -- Gene Ontology analysis
7. GO_analysis_combined_results.csv -- all GO found for the differential expressed genes
8. KEGG.R -- Kyoto Encyclopedia of genes and genomes pathway analysis
9. KEGGresults.csv -- all KEGG found for the differential expressed genes
10. AllEggSacsdata.csv -- data for survival and fat content 
11. Survival_Fat_models.R -- scripts survival proportion and fat content models
12. FA_results.csv -- Fatty acids per sample 
13. FA_analysis.R -- Script for fatty acid analysis 
