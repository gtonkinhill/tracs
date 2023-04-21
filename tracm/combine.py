import os
import sys
import re
import argparse
import glob
import gzip
import math
from collections import defaultdict
from collections import ChainMap
import pyfastx
from joblib import Parallel, delayed


def combine_parser(parser):
    parser.description = "Combine runs of Tracm'm align ready for distance estimation"

    io_opts = parser.add_argument_group("Input/output")

    io_opts.add_argument(
        "-i",
        "--input",
        dest="directories",
        required=True,
        help="Paths to each directory containing the output of the Trac'm align function",
        type=os.path.abspath,
        nargs="+",
    )

    io_opts.add_argument(
        "-o",
        "--output",
        dest="output_dir",
        required=True,
        help="name of the output driectory to store the combined alignments.",
        type=str,
    )

    # Other options
    parser.add_argument(
        "-t",
        "--threads",
        dest="n_cpu",
        help="number of threads to use (default=1)",
        type=int,
        default=1,
    )

    parser.add_argument(
        "--quiet",
        dest="quiet",
        help="turns off some console output",
        action="store_true",
        default=False,
    )

    parser.set_defaults(func=combine)

    return parser


def find_ref(filename):
    pattern = r"posterior_counts_ref_(.+)\.fasta"
    result = re.search(pattern, filename)

    if result:
        ref = result.group(1)
    else:
        print("ERROR: {} is not the expected output of Trac'm align".format(filename))
        sys.exit(1)

    return ref

def sum_after_semicolon(line):
    last_column = line.strip().split()[-1]
    numbers_str = last_column.replace(':',',')
    numbers = map(int, numbers_str.split(',')[2:])
    return sum(numbers)

def calculate_coverage(pileup):
    sample = os.path.dirname(pileup).split(os.sep)[-1]
    ref = re.search(r"ref_(.+)_pileup", os.path.basename(pileup)).group(1)

    try:
        with gzip.open(pileup, "rt") as infile:
            cov = 0
            depth = 0
            
            for line in infile:
                c = sum_after_semicolon(line)
                if c>0 : cov += 1
                depth += c
    except EOFError as e:
        print(str(e))
        print(f"Error: An EOFError occurred reading {pileup}")
        return (sample, ref, math.nan, math.nan)

    return (sample, ref, cov, depth)

def combine(args):
    # check that all directories exist
    if len(args.directories) == 1:
        with open(args.directories[0], "r") as infile:
            args.directories = [line.strip() for line in infile.readlines()]
        
    for directory in args.directories:
        if not os.path.isdir(directory):
            print("ERROR: {} is not a directory".format(directory))
            sys.exit(1)

    # create directory if it isn't present already
    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)
    # make sure trailing forward slash is present
    args.output_dir = os.path.join(args.output_dir, "")

    # collect alignments by reference genome
    alignments = defaultdict(list)
    for directory in args.directories:
        sample = os.path.basename(os.path.normpath(directory))
        for aln in glob.iglob(os.path.join(directory, "*posterior_counts_ref_*.fasta")):
            ref = find_ref(aln)
            alignments[ref].append((sample, aln))

    # write out as gzipped multifasta files and calculate fraction of N's
    ncovs = Parallel(n_jobs=args.n_cpu)(
        delayed(write_alignment)(ref, alns, args.output_dir, args.quiet)
            for ref, alns in alignments.items()
    )
    ncovs = ChainMap(*ncovs)

    # calculate coverage
    coverage = Parallel(n_jobs=args.n_cpu)(
        delayed(calculate_coverage)(pileup)
            for directory in args.directories
            for pileup in glob.iglob(os.path.join(directory, "*_pileup.txt.gz"))
    )

    # merge and divide by length
    coverage_dict = {}
    for sample, ref, cov, depth in coverage:
        if (sample, ref) in ncovs:
            coverage_dict[(sample, ref)] = (cov/ncovs[(sample, ref)][1], depth/cov, depth/ncovs[(sample, ref)][1])
        else:
            coverage_dict[(sample, ref)] = ("NA", depth/cov, "NA")


    # process sourmash results and write to output file
    with open(args.output_dir + "combined_metadata.csv", "w") as outfile:
        outfile.write("sample,accession,intersect_bp,f_orig_query,f_match,f_unique_to_query,coverage,mean_depth,mean_nonzero_depth,frac_N,species\n")
        for directory in args.directories:
            sample = os.path.basename(os.path.normpath(directory))
            for sourmash in glob.iglob(os.path.join(directory, "*_sourmash_hits.csv")):
                with open(sourmash, "r") as infile:
                    next(infile) # skip header
                    for line in infile:
                        line = line.strip().split(",")
                        accession = line[9].split()[0].strip('"')
                        
                        species = line[9].replace(accession, '').replace('"',).strip()
                        
                        if (sample, accession) in coverage_dict:
                            cov = coverage_dict[(sample, accession)]
                            if (sample, accession) in ncovs:
                                ncov = str(ncovs[(sample, accession)][0])
                            else:
                                ncov = "NA"
                            outfile.write(','.join([sample, accession] + line[:4] + [str(cov[0]), str(cov[1]), str(cov[2]), ncov, species]) + '\n')
                        else:
                            outfile.write(','.join([sample, accession] + line[:4] + ["NA", "NA", "NA", "NA",species]) + '\n')
                        

    return


def write_alignment(ref, alns, output_dir, quiet):
    output_file = output_dir + ref + "_combined.fasta.gz"
    ncov = {}

    if not quiet:
        print("Writing combined alignment for {} to {}".format(ref, output_file))

    with gzip.open(output_file, "wt") as fasta_file:
        for aln in alns:
            count = 0
            for name, seq in pyfastx.Fasta(aln[1], build_index=False):
                fasta_file.write(f">{aln[0]}\n{seq}\n")
                count += 1
                if count > 1:
                    print("ERROR: {} contains more than one sequence".format(aln[1]))
                    sys.exit(1)
                ncov[(aln[0], ref)] = (seq.count("N") / len(seq), len(seq))

    return ncov


def main():
    # set up and parse arguments
    parser = argparse.ArgumentParser()
    parser = combine_parser(parser)
    args = parser.parse_args()

    # run distance command
    args.func(args)

    return


if __name__ == "__main__":
    main()
