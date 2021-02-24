#!/bin/bash

BAM=$1
THREADS=$2

module load SAMtools/1.5

samtools flagstat -@ $THREADS $BAM | grep total | awk '{print $1}'
