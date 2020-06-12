# OE_TRPM5

Analysis code will be available by time of publication

Data will be public via NCBI SRA by time of publication under project accession: PRJNA632936

BioRXIV link: https://www.biorxiv.org/content/10.1101/2020.05.14.096016v2

# Methods overview  

Demultiplexed paired-end FASTQ were obtained from the University of Colorado Genomics and Microarray Core. Adapter sequences were removed and low quality bases were trimmed using BBDuk (BBTools).

## Differential Gene Expression

Transcript quantification was performed on trimmed FASTQ files with Salmon 1.2.1  using a decoy-away transcriptome index of Ensembl GRCm38 release 99. Prior to indexing, the raw FASTA and GTF files were modified to include GFP and mCherry transcripts.  Salmon output was summarized at the gene level using Tximport in R. The Tximport object was used as the imput for DESeq2. Differential expression was performed using default parameters and pairwise comparisons between each of the three cell-type classifications.  To determine the effect of sex, a separate analysis was performed with the model design of 'celltype + sex + celltype::sex.' Significance was determined using alpha = 0.05.  


## Trpm5 Gene coverage

Read mapping was performed on trimmed FASTQ files with STAR 2.7.3 using a genome index of Ensembl GRCm38 release 99 (with addition of GFP and mCherry fasta). STAR output were filtered for primary mapped reads covering the region of the Trpm5 gene using 'Samtools.' To get per base pair information, these BAMs were passed to Samtools 'depth.' For average exon coverage, a customed BED file was generated from the Trpm5 transcript ENSMUST00000009390 'EXON' entries of the GRCm39 release 99 GTF file. The bedfile was used with Samtools 'bedcov' to summarize read coverage over each exon.  Results were parsed and summarized in R across all samples.  95% bootstrap confidence intervals were calculated using the percentile method.

detailed data:
1) trimming/filtering (trim_filter.sh)
2) Salmon (salmon_quant.sh)
3) STAR (STAR_align.sh)
4) Samtools (depth_trpm5.sh)
5) DESEQ (tximport_deseq2_clean.R)
6) Depth (trpm5_coverage.R)
7) TopGO (topGO_clean.R)
6) plotting (plots_clean.R)

