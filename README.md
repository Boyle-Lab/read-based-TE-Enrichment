# read-based-TE-Enrichment
Repeat enrichment testing for short-read sequencing data using an assembly-free alignment approach.


# Prerequisites

1. BWA (<http://bio-bwa.sourceforge.net/>)
2. BEDTools (<https://bedtools.readthedocs.io/en/latest/index.html>)
3. SAMtools (<http://www.htslib.org/>)
4. Picard Tools (<https://broadinstitute.github.io/picard/>)
5. Python >= 3.8
6. R >= 4.0

# Preparing to Run the Pipeline

The pipeline works by aligning reads to an artificial genome constructed from a set of repeat consensus sequences and spacer sequences concatenated into pseudochromosomes. This genome is then indexed and used as an alignment reference for the BWA mem algorithm. Reads that map entirely or partially to any coordinates within the artificial genome can be assumed to originate from a copy of a known repeat sequence within the host genome. Matches to each type of repeat are then labeled and counted. When applied to immunoprecipitated (IP) and Input samples, enrichment tests can then be performed comparing the observed (from the IP sample) and expected (from the input sample) count densities. Fisher's exact tests are used by default. 

The first step is to obtain/generate an artificial genome. Several examples are given in the "genomes" folder, or you can create your own by obtaining a repeat library from RepeatMasker, DFAM, etc., and using the "extractAncestral.py" and "makePseudoChromosomes.py" scripts, available as part of the repeatMaskerUtils repository, available [here](https://github.com/adadiehl/repeatMaskerUtils.git). This must then be indexed with bwa index before use.

# Usage

The enrichment test workflow is invoked with the run_enrichment_pipeline.py script. This script performs all necessary alignment, filtering, and analysis steps for a set of immunoprecipitated and control (input) reads. The workflow can handle either single-ended or paired-ended samples, and is capable of integrating data from multiple experimental replicates.

```bash
run_enrichment_pipeline.py [-h] -f ip_sample.fastq
                                  [ip_sample.fastq ...]
                                  [-g ip_sample.fastq [ip_sample.fastq ...]]
                                  -b input_sample.fastq
                                  [input_sample.fastq ...]
                                  [-c input_sample.fastq [input_sample.fastq ...]]
                                  -n genome.fa [-a /path/to/write/alignments]
                                  [-r /path/to/write/results] [-p]
                                  [-t THREADS] [-q N] [-o OUTROOT]
                                  [-s SCRIPTSDIR] [-l CMDLOG]
```

## Required Arguments
Short | Long | Arg(s) | Description
------|------|--------|-------------
-f | --fgSamples | ip_sample.1.fastq [ip_sample.2.fastq ...] | Immunoprecipitated sample fastq file(s) for single-ended or read pair 1. Multiple files will be processed as replicates.
-b | --bgSamples | input_sample.1.fastq [input_sample.2.fastq ...] | Input sample fastq file(s) for single-ended or read pair 1. Multiple files will be processed as replicates.
-n | --genome | genome.fa | Fasta file containing the pseudogenome to which samples will be aligned.

## For Paired-Ended Experiments
Short | Long | Arg(s) | Default | Description
------|------|--------|---------|-------------
-g | --fgSamples2 | ip_sample.1.fastq [ip_sample.2.fastq ...] | Immunoprecipitated sample fastq file(s) for read pair 2. Multiple files will be processed as replicates.
-c | --bgSamples2 | input_sample.1.fastq [input_sample.2.fastq ...] | Input sample fastq file(s) for read pair 2. Multiple files will be processed as replicates.

## Optional Arguments
Short | Long | Arg(s) | Default | Description 
------|------|--------|---------|-------------
-a | --alignmentPath | /path/to/write/alignments | . | Location to write alignment data.
-r | --resultsPath | /path/to/write/results | . | Location to write results.
-p | --paired | | False | Toggle paired-end mode.
-t | --threads | N | 1 | Maximum number of threads to use for subprocesses.
-q | --minQual | N | 30 | Minimum quality score for retaining mapped reads.
-o | --outRoot | OUTROOT | basename of the first IP fastq | Output root to prepend to results.
-s | --scriptsDir | SCRIPTSDIR | ./scripts | Path to shell scripts.
-l | --cmdLog | CMDLOG | run_enrichment_pipeline.cmd.log | File to which command output will be logged.


