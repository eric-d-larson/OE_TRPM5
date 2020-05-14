#!/usr/bin/env bash

bbduk -Xmx7g t=4 \
    in1=R1.fastq.gz \
    in2=R2.fastq.gz \
    out1=R1trim.fastq.gz \
    out2=R2trim.fastq.gz \
    ref=adapters.fa \
    forcetrimleft=1 \
    ktrim=r \
    k=23 \
    mink=11 \
    hdist=1 \
    tpe=t \
    tbo=t \
    qtrim=rl \
    trimq=5 \
    minlength=50
