#!/bin/bash                                                                                                      

# Load software modules. Comment these lines if your system does not have an LMOD type module system.
module load BEDTools/2.26.0

WD=$(pwd)

GENOME=$1
SAMPLE=$2
BWA_DIR=$3
RESULTS_DIR=$4

# Map against the reference bed to total up reads that map to each TE consensus
bedtools coverage -a $(echo $GENOME | sed 's/fa$/bed/') -b $BWA_DIR/$SAMPLE.pruned.bam > $RESULTS_DIR/$SAMPLE.pruned.te-counts.txt
