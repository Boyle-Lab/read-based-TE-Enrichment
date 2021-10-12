#!/bin/bash                                                                                                      

# Load software modules. Comment these lines if your system does not have an LMOD type module system.
module load BWA/0.7.16a
module load SAMtools/1.5
module load BEDTools/2.26.0
module load picard/2.18.0

WD=$(pwd)

FASTQ_1=$1
FASTQ_2=$2
THREADS=$3
MAPQ=$4
GENOME=$5
SAMPLE=$6
BWA_DIR=$7
RESULTS_DIR=$8

# Run bwa mem and sort output
bwa mem -t $(( $THREADS / 2 )) $GENOME $FASTQ_1 $FASTQ_2 | samtools sort -@ $(( $THREADS / 2 )) -O bam -T $SAMPLE.sort.tmp -o $BWA_DIR/$SAMPLE.bam

# Mark duplicates with Picard
java -jar $EBROOTPICARD/picard.jar MarkDuplicates I=$BWA_DIR/$SAMPLE.bam O=$BWA_DIR/$SAMPLE.md.bam ASSUME_SORTED=true METRICS_FILE=$BWA_DIR/$SAMPLE.dup.metrics VALIDATION_STRINGENCY=LENIENT TMP_DIR=.
samtools index $BWA_DIR/$SAMPLE.md.bam

# Filter reads. 4 = unmapped, 256 = not primary, 1024 = PCR/Optical duplicate, suppplementary alignment
samtools view -b -h -F 4 -F 256 -F 1024 -F 2048 -q $MAPQ $BWA_DIR/$SAMPLE.md.bam > $BWA_DIR/$SAMPLE.pruned.bam
samtools index $BWA_DIR/$SAMPLE.pruned.bam

# Map against the reference bed to total up reads that map to each TE consensus
#bedtools coverage -a $(echo $GENOME | sed 's/fa$/bed/') -b $BWA_DIR/$SAMPLE.pruned.bam > $RESULTS_DIR/$SAMPLE.pruned.te-counts.txt
