# OE_TRPM5

Analysis code will be available by time of publication
Data will be deposited to NCBI GEO by time of publication. GEOXXXXXXX
BioRXIV link: 

Methods overview
Demultiplexed paired-end FASTQ were obtained from the University of Colorado Genomics and Microarray Core. Adapter sequences were removed and low quality bases were trimmed using BBDuk (BBTools).

Differential Gene Expression
Transcript quantification was performed on trimmed FASTQ files with Salmon 1.2.1  using a decoy-away transcriptome index of Ensembl GRCm38 release 99. Prior to indexing, the raw FASTA and GTF files were modified to include GFP and mCherry transcripts.  Salmon output was summarized at the gene level using Tximport in R. The 'counts' slot of the Tximport object was used as the imput for DESeq2. Differential expression was performed using default parameters and pairwise comparisons between each of the three cell-type classifications.  
