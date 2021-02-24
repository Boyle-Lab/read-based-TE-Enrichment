#!/usr/bin/env python3

import sys, os, re, subprocess
from argparse import ArgumentParser

def main():
    parser = ArgumentParser(description='Run enrichment tests for a set of IP and input fastq files against a library of sequences in the given artificial genome.')
    parser.add_argument('-f', '--fgSamples', metavar='ip_sample.fastq', type=str, nargs='+', required=True,
                        help='Immunoprecipitated sample fastq file(s).')
    parser.add_argument('-b', '--bgSamples', metavar='input_sample.fastq', type=str, nargs='+', required=True,
                        help='Input sample fastq file(s).')
    parser.add_argument('-g', '--genome', metavar="genome.fa", type=str, required=True,
                        help='Path to artificial genome fasta.')
    parser.add_argument('-a', '--alignmentPath', metavar="/path/to/write/alignments", type=str, default='.',
                        help='Location to write alignment data.')
    parser.add_argument('-r', '--resultsPath', metavar="/path/to/write/results", type=str, default='.',
                        help='Location to write results.')
    parser.add_argument('-p', '--paired', action='store_true',
                        help='Sequencing data are paired-ended.')
    parser.add_argument('-t', '--threads', type=int, default=8,
                        help='Max number of threads to launch')
    parser.add_argument('-q', '--minQual', metavar="N", type=str, default='30',
                        help='Minimum quality score for retaining mapped reads.')
    parser.add_argument('-o', '--outRoot', type=str,
                        help='Output root to prepend to results. By default, the basename of the first immunoprecipitated fastq file will be used.')
    parser.add_argument('-s', '--scriptsDir', type=str, default='./scripts',
                        help='Path to shell scripts.')
    parser.add_argument('-l', '--cmdLog', type=str, default='run_enrichment_pipeline.cmd.log',
                        help="File to write command output to.")

    args = parser.parse_args()

    # Sanity checks
    if args.paired and len(args.fgSamples) < 2:
        parser.exit(status=2, message="ERROR: --paired requires two immunoprecipitated files!\n")
    if len(args.fgSamples) == 2 and len(args.bgSamples) == 2 and not args.paired:
        args.paired = True
        sys.stderr.write("WARNING: Two sets of input and immunoprecipitated samples provided. Assuming paired-ended data!\n")
    if args.paired and len(args.bgSamples) < 2:
        parser.exit(status=2, message="ERROR: --paired requires two input files!\n")
    if len(args.fgSamples) != len(args.bgSamples):
        parser.exit(status=2, message="ERROR: Immunoprecipitated and input file lists must be of the same length!\n")        

    # See if the outRoot arg is set. If not, use the basename of the first fastq file.
    if args.outRoot is None:
        args.outRoot = os.path.basename(os.path.splitext(args.fgSamples[0])[0])
        
    # We will utilize external tools to run each step of the pipeline through system calls.

    # Alignment/coverage pipeline for IP and Input samples
    if args.paired:
        cmd_args = [os.path.join(args.scriptsDir, "paired-ended-job.sh"), args.fgSamples[0], args.fgSamples[1], str(args.threads), str(args.minQual), args.genome, args.outRoot + ".fg", args.alignmentPath, args.resultsPath, "2>>" + args.cmdLog, "1>&2"]
        sys.stderr.write("\nExecuting alingment, indexing, and filtering steps for IP samples: {}\n\n".format(" ".join(cmd_args)))
        subprocess.run(cmd_args)
        
        cmd_args = [os.path.join(args.scriptsDir, "paired-ended-job.sh"), args.bgSamples[0], args.bgSamples[1], str(args.threads), str(args.minQual), args.genome, args.outRoot + ".bg", args.alignmentPath, args.resultsPath, "2>>" + args.cmdLog, "1>&2"]
        sys.stderr.write("\nExecuting alingment, indexing, and filtering steps for Input samples: {}\n\n".format(" ".join(cmd_args)))
        subprocess.run(cmd_args)
    else:
        cmd_args = [os.path.join(args.scriptsDir, "single-ended-job.sh"), args.fgSamples[0], str(args.threads), str(args.minQual), args.genome, args.outRoot + ".fg", args.alignmentPath, args.resultsPath, "2>>" + args.cmdLog, "1>&2"]
        sys.stderr.write("\nExecuting alingment, indexing, and filtering steps for IP sample: {}\n\n".format(" ".join(cmd_args)))
        subprocess.run(cmd_args)

        cmd_args = [os.path.join(args.scriptsDir, "single-ended-job.sh"), args.bgSamples[0], str(args.threads), str(args.minQual), args.genome, args.outRoot + ".bg", args.alignmentPath, args.resultsPath, "2>>" + args.cmdLog, "1>&2"]
        sys.stderr.write("\nExecuting alingment, indexing, and filtering steps for Input sample: {}\n\n".format(" ".join(cmd_args)))
        subprocess.run(cmd_args)
    
    # Get read counts from IP and input alignments
    sys.stderr.write("\nCalculating total reads for IP and Input samples\n\n")
    fgBam = os.path.join(args.alignmentPath, args.outRoot) + ".fg.pruned.bam"
    fgReads = subprocess.run([os.path.join(args.scriptsDir, "get_readcount.sh"), fgBam, str(args.threads), "2>>" + args.cmdLog, "1>&2"])
    bgBam = os.path.join(args.alignmentPath, args.outRoot) + ".bg.pruned.bam"
    fgReads = subprocess.run([os.path.join(args.scriptsDir, "get_readcount.sh"), bgBam, str(args.threads), "2>>" + args.cmdLog, "1>&2"])

    # Run enrichment tests on the output. Results are written from R.
    sys.stderr.write("\nRunning enrichment tests.\n\n")
    fgCounts = os.path.join(args.resultsPath, args.outRoot) + ".fg.pruned.te-counts.txt"
    bgCounts = os.path.join(args.resultsPath, args.outRoot) + ".bg.pruned.te-counts.txt"
    outFile = os.path.join(args.resultsPath, args.outRoot) + ".enrichment-tests.txt"
    subprocess.run([os.path.join(args.scriptsDir, "calc_enrichments.Rscript"), fgCounts, bgCounts, fgReads, bgReads, outFile, "2>>" + args.cmdLog])
    
    
if __name__ == "__main__":
    main()
