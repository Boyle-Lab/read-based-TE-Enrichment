# read-based-TE-Enrichment
Repeat enrichment testing for short-read sequencing data using an assembly-free alignment approach.


# Prerequsisites

BWA
BEDTools
SAMtools
Python >= 3.7

# Running the Pipeline

The pipeline works by aligning reads to an artificial genome constructed from a set of repeat consensus sequences and spacer sequences concatenated into pseudochromosomes. This genome is then indexed and used as an alignment reference for the BWA mem algorithm. Reads that map entirely or partially to any coordinates within the artificial genome can be assumed to originate from a copy of a known repeat sequence within the host genome. Matches to each type of repeat are then labeled and counted. When applied to immunoprecipitated (IP) and Input samples, enrichment tests can then be performed comparing the observed (from the IP sample) and expected (from the input sample) count densities. Fisher's exact tests are used by default. 

The first step is to obtain/generate an artificial genome. Several examples are given in the "genomes" folder, or you can create your own by obtaining a repeat library from RepeatMasker, DFAM, etc., and using the "extractAncestral.py" and "makePseudoChromosomes.py" scripts, available as part of the repeatMaskerUtils repository, available [here](https://github.com/adadiehl/repeatMaskerUtils.git). This must then be indexed with bwa index before use.
