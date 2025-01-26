## Environment information

Analysis was performed on the University of Greifswald's "High-Performance Computing" (HPC) cluster.

# 1. Merge lanes of raw Illumina reads 

I copied the file to where the raw data are and ran it from there with the command sh.

```
#!/bin/bash

for i in $(find ./ -type f -name "*.fastq.gz" | while read F; do basename $F | rev | cut -c 22- | rev; done | sort | uniq)

    do echo "Merging R1"

cat "$i"_L00*_R1_001.fastq.gz > "$i"_R1.merged.fq.gz

       echo "Merging R2"

cat "$i"_L00*_R2_001.fastq.gz > "$i"_R2.merged.fq.gz

done;
```
# 2. FASTQC of raw reads 

```
#!/usr/bin/env bash

fastqc \
  -o ./raw_qc \
  -t 2 \
  ~/Ilumina_raw/Alignment_1/fastq1/Fastq/S324_S2_R1.merged.fq.gz \
  ~/Ilumina_raw/Alignment_1/fastq1/Fastq/S324_S2_R2.merged.fq.gz \
  ~/Ilumina_raw/Alignment_1/fastq1/Fastq/S325_S3_R1.merged.fq.gz \
  ~/Ilumina_raw/Alignment_1/fastq1/Fastq/S325_S3_R2.merged.fq.gz \
  ~/Ilumina_raw/Alignment_1/fastq1/Fastq/S326_S4_R1.merged.fq.gz \
  ~/Ilumina_raw/Alignment_1/fastq1/Fastq/S326_S4_R2.merged.fq.gz \
  ~/Ilumina_raw/Alignment_1/fastq1/Fastq/S327_S5_R1.merged.fq.gz \
  ~/Ilumina_raw/Alignment_1/fastq1/Fastq/S327_S5_R2.merged.fq.gz \
  ~/Ilumina_raw/Alignment_1/fastq1/Fastq/S328_S6_R1.merged.fq.gz \
  ~/Ilumina_raw/Alignment_1/fastq1/Fastq/S328_S6_R2.merged.fq.gz \
  ~/Ilumina_raw/Alignment_1/fastq1/Fastq/S329_S7_R1.merged.fq.gz \
  ~/Ilumina_raw/Alignment_1/fastq1/Fastq/S329_S7_R2.merged.fq.gz \
  ~/Ilumina_raw/Alignment_1/fastq1/Fastq/S330_S8_R1.merged.fq.gz \
  ~/Ilumina_raw/Alignment_1/fastq1/Fastq/S330_S8_R2.merged.fq.gz \
  ~/Ilumina_raw/Alignment_1/fastq1/Fastq/S331_S9_R1.merged.fq.gz \
  ~/Ilumina_raw/Alignment_1/fastq1/Fastq/S331_S9_R2.merged.fq.gz \
  ~/Ilumina_raw/Alignment_1/fastq1/Fastq/S332_S10_R1.merged.fq.gz \
  ~/Ilumina_raw/Alignment_1/fastq1/Fastq/S332_S10_R2.merged.fq.gz \
  ~/Ilumina_raw/Alignment_1/fastq1/Fastq/S333_S11_R1.merged.fq.gz \
  ~/Ilumina_raw/Alignment_1/fastq1/Fastq/S333_S11_R2.merged.fq.gz \
  ~/Ilumina_raw/Alignment_1/fastq1/Fastq/S334_S12_R1.merged.fq.gz \
  ~/Ilumina_raw/Alignment_1/fastq1/Fastq/S334_S12_R2.merged.fq.gz
 ```
# 3. Genome Alignment with STAR 2.7.11a 
### Generate genome indices: 
Download the genome fasta file: .fna and the genome annotation file .GTF
1. I went into https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_947563725.1/
2. NCBI RefSeq assembly GCF_947563725.1 --> click on actions --> see more files on FTP
3. right click on  GCF_947563725.1_qqArgBrue1.1_genomic.gtf.gz--> copy link address
4. go to the terminal --> in the folder StarAlignment -->
```
weget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/947/563/725/GCF_947563725.1_qqArgBrue1.1/GCF_947563725.1_qqArgBrue1.1_genomic.gtf.gz 
 ```
5. back to NCBI submitted GenBank assembly GCA_947563725.1 --> click on actions --> see more files on FTP
6. right click on GCF_947563725.1_qqArgBrue1.1_genomic.fna.gz --> click on actions --> copy link address
8. go to the terminal --> in the folder StarAlignment -->
```
weget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/947/563/725/GCA_947563725.1_qqArgBrue1.1/GCF_947563725.1_qqArgBrue1.1_genomic.fna.gz
```
10. Unzip the files by gunzip *.gz
11. change the .fna file to fa by cp GCF_947563725.1_qqArgBrue1.1_genomic.fna GCF_947563725.1_qqArgBrue1.1_genomic.fa
12. make a new directory to keep stuff organized for star: mkdir ref
```
STAR --runMode genomeGenerate --genomeDir ref/ --genomeFastaFiles GCF_947563725.1_qqArgBrue1.1_genomic.fa --sjdbGTFfile GCF_947563725.1_qqArgBrue1.1_genomic.gtf --runThreadN 16
```
### Aligning reads 
make a directory for the RNA fasta reads 
mkdir fastq  

```
#!/bin/bash

for base in S323_S1 S324_S2 S325_S3 S326_S4 S327_S5 S328_S6 S329_S7 S330_S8 S331_S9 S332_10 S333_S11 S334_S12
do
  echo $base

  # define R1 fq filename
  fq1=${base}_R1.merged.fq.gz

 # define R2 fastq filename
  fq2=${base}_R2.merged.fq.gz

 # align with STAR

 STAR --runMode alignReads \
 --genomeDir ../ref/ \
 --outSAMtype BAM SortedByCoordinate \
 --readFilesIn $fq1 $fq2 --runThreadN 12 \
 --readFilesCommand zcat \
 --outFileNamePrefix ../mapped/$base"_" 
 done

echo "done!"
```
# 4. Create the counts table with Subread 2.0.6 using FeatureCounts

```
featureCounts -p --countReadPairs -g gene_id -a GCF_947563725.1_qqArgBrue1.1_genomic.gtf -o count.out -T 8 bams/*.bam
```
# 5. Count table
I downloaded the count.out table and organized it in Excel with my sample names: 
Cold_S1	Cold_S2	Cold_S3	Cold_S4_ Moderate_S5 Moderate_S6 Moderate_S7 Moderate_S8 Warm_S9 Warm_S10 Warm_S11 Warm_S12
and save it as .csv

For the next

