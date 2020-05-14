#!/usr/bin/env bash

salmon quant -i GRCm38_decoy_idx -l A \
    --validateMappings \
    -1 R1trim.fastq.gz \
    -2 R2trim.fastq.gz \
    -p 12 -o salmon_results/sample
