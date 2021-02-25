#!/usr/bin/env python3

import sys, os, re, subprocess
from argparse import ArgumentParser

def constructReadCountsStr(samplesList, resultsPath, outRoot, stype):
    """ Construct a comma-delimited string of counts files. """
    files = []
    for i in range(len(samplesList)):
        files.append(os.path.join(resultsPath, outRoot) + "." + stype + "." + str(i) + ".pruned.te-counts.txt")
    return ",".join(files)

def runAlignments(samplesList, samplesList2, scriptsDir, threads, minQual, genome, alingmentPath, outRoot, stype, isPaired):
    """ Run the alignment/filtering pipeline for a set of samples. """
    for i in range(len(samplesList)):
        if isPaired:
            cmd_args = [os.path.join(scriptsDir, "paired-ended-job.sh"), samplesList[i], samplesList2[i], str(threads), str(minQual), genome, outRoot + "." + stype + "." + str(i), alignmentPath, resultsPath]
        else:
            cmd_args = [os.path.join(scriptsDir, "single-ended-job.sh"), samplesList[i], str(threads), str(minQual), genome, outRoot + "." + stype + "." + str(i), alignmentPath, resultsPath]
        sys.stderr.write("\nExecuting command: {}\n\n".format(" ".join(cmd_args)))
        subprocess.run(cmd_args)


def gatherReadCounts(samplesList, scriptsDir, threads, alignmentPath, outRoot, stype):
    """ Collect total read counts for filtered reads. """
    reads = 0
    for i in range(len(samplesList)):
        bam = os.path.join(alignmentPath, outRoot) + "." + stype + "." + str(i) + ".pruned.bam"
        reads += subprocess.run([os.path.join(args.scriptsDir, "get_readcount.sh"), bam, str(threads)], capture_output=True, text=True).stdout.strip("\n")
    return reads


def main():
    parser = ArgumentParser(description='Run enrichment tests for a set of IP and input fastq files against a library of sequences in the given artificial genome.')
    parser.add_argument('-f', '--fgSamples', metavar='ip_sample.fastq', type=str, nargs='+', required=True,
                        help='Immunoprecipitated sample fastq file(s) for single-ended or read pair 1.')
    parser.add_argument('-g', '--fgSamples2', metavar='ip_sample.fastq', type=str, nargs='+',
                        help='Immunoprecipitated sample fastq file(s) for read pair 2.')
    parser.add_argument('-b', '--bgSamples', metavar='input_sample.fastq', type=str, nargs='+', required=True,
                        help='Input sample fastq file(s) for single-ended or read pair 1.')
    parser.add_argument('-c', '--bgSamples2', metavar='input_sample.fastq', type=str, nargs='+',
                        help='Input sample fastq file(s) for read pair 2.')
    parser.add_argument('-n', '--genome', metavar="genome.fa", type=str, required=True,
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
    if (args.paired and args.fgSamples2 == None) or (len(args.fgSamples) != len(args.fgSamples2)):
        parser.exit(status=2, message="ERROR: --paired mode requires an equal number of files for pair1 and pair2 immunoprecipitated sample(s)!\n")
    if (args.paired and args.bgSamples2 == None) or (len(args.bgSamples) != len(args.bgSamples2)):
        parser.exit(status=2, message="ERROR: --paired mode requires an equal number of files for pair1 and pair2 input sample(s)!\n")
    if (args.fgSamples2 != None and args.fgSamples != None) and (args.bgSamples2 != None and args.bgSamples != None) and not args.paired:
        args.paired = True
        sys.stderr.write("WARNING: Paired input and immunoprecipitated samples provided. Forcing paired-ended mode!\n")
    
    # See if the outRoot arg is set. If not, use the basename of the first fastq file.
    if args.outRoot is None:
        args.outRoot = os.path.basename(os.path.splitext(args.fgSamples[0])[0])
        
    # We will utilize external tools to run each step of the pipeline through system calls.
    # Alignment/coverage pipeline for IP and Input samples
    runAlignments(args.fgSamples, args.fgSamples2, args.scriptsDir, args.threads, args.minQual, args.genome, args.alignmentPath, args.outRoot, "fg", args.paired)
    runAlignments(args.bgSamples, args.bgSamples2, args.scriptsDir, args.threads, args.minQual, args.genome, args.alignmentPath, args.outRoot, "bg", args.paired)
    
    """
    if args.paired:
        for i in range(len(args.fgSamples)):
            cmd_args = [os.path.join(args.scriptsDir, "paired-ended-job.sh"), args.fgSamples[i], args.fgSamples2[i], str(args.threads), str(args.minQual), args.genome, args.outRoot + ".fg." + str(i), args.alignmentPath, args.resultsPath, "2>>" + args.cmdLog, "1>&2"]
            sys.stderr.write("\nExecuting alingment, indexing, and filtering steps for IP samples: {}\n\n".format(" ".join(cmd_args)))
            subprocess.run(cmd_args)

        for i in range(len(args.bgSamples)):
            cmd_args = [os.path.join(args.scriptsDir, "paired-ended-job.sh"), args.bgSamples[i], args.bgSamples2[i], str(args.threads), str(args.minQual), args.genome, args.outRoot + ".bg." + str(i), args.alignmentPath, args.resultsPath, "2>>" + args.cmdLog, "1>&2"]
            sys.stderr.write("\nExecuting alingment, indexing, and filtering steps for Input samples: {}\n\n".format(" ".join(cmd_args)))
            subprocess.run(cmd_args)

    else:
        for i in range(len(args.fgSamples)):
            cmd_args = [os.path.join(args.scriptsDir, "single-ended-job.sh"), args.fgSamples[i], str(args.threads), str(args.minQual), args.genome, args.outRoot + ".fg." + str(i), args.alignmentPath, args.resultsPath, "2>>" + args.cmdLog, "1>&2"]
            sys.stderr.write("\nExecuting alingment, indexing, and filtering steps for IP sample: {}\n\n".format(" ".join(cmd_args)))
            subprocess.run(cmd_args)
            
        for i in range(len(args.bgSamples)):
            cmd_args = [os.path.join(args.scriptsDir, "single-ended-job.sh"), args.bgSamples[i], str(args.threads), str(args.minQual), args.genome, args.outRoot + ".bg." + str(i), args.alignmentPath, args.resultsPath, "2>>" + args.cmdLog, "1>&2"]
            sys.stderr.write("\nExecuting alingment, indexing, and filtering steps for Input sample: {}\n\n".format(" ".join(cmd_args)))
            subprocess.run(cmd_args)
    """
    
    # Get read counts from IP and input alignments
    sys.stderr.write("\nCalculating total reads for IP and Input samples\n\n")
    fgReads = gatherReadCounts(args.fgSamples, args.scriptsDir, args.threads, args.alignmentPath, args.outRoot, "fg")
    bgReads = gatherReadCounts(args.bgSamples, args.scriptsDir, args.threads, args.alignmentPath, args.outRoot, "bg")
    """
    for i in range(len(args.fgSamples)):
        fgBam = os.path.join(args.alignmentPath, args.outRoot) + ".fg." + str(i) + ".pruned.bam"
        fgReads += subprocess.run([os.path.join(args.scriptsDir, "get_readcount.sh"), fgBam, str(args.threads), "2>>" + args.cmdLog], capture_output=True, text=True).stdout.strip("\n")
    for i in range(len(args.bgSamples)):
        bgBam = os.path.join(args.alignmentPath, args.outRoot) + ".bg." + str(i) + ".pruned.bam"
        bgReads += subprocess.run([os.path.join(args.scriptsDir, "get_readcount.sh"), bgBam, str(args.threads), "2>>" + args.cmdLog], capture_output=True, text=True).stdout.strip("\n")
    """
    sys.stderr.write("FG Reads: {}\nBG Reads: {}\n".format(fgReads, bgReads))

    # Run enrichment tests on the output. Results are written from R.
    sys.stderr.write("\nRunning enrichment tests.\n\n")
    fgCounts = constructReadCountsStr(args.fgSamples, args.resultsPath, args.outRoot, "fg")
    bgCounts = constructReadCountsStr(args.bgSamples, args.resultsPath, args.outRoot, "bg")
    outFile = os.path.join(args.resultsPath, args.outRoot) + ".enrichment-tests.txt"
    cmd_args = ["Rscript", "--vanilla", os.path.join(args.scriptsDir, "calc_enrichments.Rscript"), fgCounts, bgCounts, fgReads, bgReads, outFile, "2>>" + args.cmdLog]
    sys.stderr.write("Invoking Rscript with the following command: {}\n".format(" ".join(cmd_args)))
    subprocess.run(cmd_args)
    
if __name__ == "__main__":
    main()
