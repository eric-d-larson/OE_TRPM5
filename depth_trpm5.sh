#!/usr/bin/env bash

samtools index STARout.bam

samtools view -F 12 -F 256 STARout.bam "7:143069153-143094642" -o trpm5_primary.bam

samtools index trpm5_primary.bam

samtools depth -a trpm5_primary.bam > trpm5.coverage

samtools bedcov trpm5_exons.bed trpm5_primary.bam > trpm5_exon.coverage


