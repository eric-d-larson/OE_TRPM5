# OE_TRPM5

Analysis code will be available by time of publication
Data will be deposited to NCBI GEO by time of publication. GEOXXXXXXX
BioRXIV link: 

Methods overview  

Demultiplexed paired-end FASTQ were obtained from the University of Colorado Genomics and Microarray Core. Adapter sequences were removed and low quality bases were trimmed using BBDuk (BBTools).

Differential Gene Expression

Transcript quantification was performed on trimmed FASTQ files with Salmon 1.2.1  using a decoy-away transcriptome index of Ensembl GRCm38 release 99. Prior to indexing, the raw FASTA and GTF files were modified to include GFP and mCherry transcripts.  Salmon output was summarized at the gene level using Tximport in R. The 'counts' slot of the Tximport object was used as the imput for DESeq2. Differential expression was performed using default parameters and pairwise comparisons between each of the three cell-type classifications.  To determine the effect of sex, a separate analysis was performed with the model design of 'celltype + sex + celltype::sex.' Significance was determined using alpha = 0.05.  

Trpm5 Gene coverage
Read mapping was performed on trimmed FASTQ files with STAR 2.7.3 using a genome index of Ensembl GRCm38 release 99 (with addition of GFP and mCherry fasta). STAR output were filtered for primary mapped reads covering the region of the Trpm5 gene using 'Samtools.' To get per base pair information, these BAMs were passed to Samtools 'depth.' For average exon coverage, a customed BED file was generated from the Trpm5 transcript ENSMUST00000009390 'EXON' entries of the GRCm39 release 99 GTF file. The bedfile was used with Samtools 'bedcov' to summarize read coverage over each exon.  Results were parsed and summarized in R across all samples.  95% bootstrap confidence intervals were calculated using the percentile method.

detailed data:
1) trimming/filtering
2) Salmon 
3) STAR
4) Samtools
5) DESEQ
6) plotting
  a) heatmaps
  b) volcano plot
  c) coverage plot
  d) bar graph

