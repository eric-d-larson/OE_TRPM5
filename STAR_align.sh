#!/usr/bin/env bash

    STAR --outFilterType BySJout \
	--runThreadN 12 \
	--outSAMtype BAM SortedByCoordinate \
	--genomeDir GRCm38_star_idx \
	--sjdbOverhang 100 \
	--quantMode GeneCounts \
	--sjdbGTFfile GRCm38_gtf \
	--outFileNamePrefix sample. \
	--readFilesCommand zcat \
	--readFilesIn R1trim.fastq.gz R2trim.fastq.gz

