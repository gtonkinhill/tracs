import os
import sys
import re
import argparse
import glob
import gzip
from collections import defaultdict
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


def combine(args):
    # check that all directories exist
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

    # write out as gzipped multifasta files
    Parallel(n_jobs=args.n_cpu)(
        delayed(write_alignment)(ref, alns, args.output_dir, args.quiet)
            for ref, alns in alignments.items()
    )

    return


def write_alignment(ref, alns, output_dir, quiet):
    output_file = output_dir + ref + "_combined.fasta.gz"

    with gzip.open(output_file, "wt") as fasta_file:
        for aln in alns:
            count = 0
            for name, seq in pyfastx.Fasta(aln[1], build_index=False):
                fasta_file.write(f">{aln[0]}\n{seq}\n")
                count += 1
                if count > 1:
                    print("ERROR: {} contains more than one sequence".format(aln[1]))
                    sys.exit(1)

    if not quiet:
        print("Writing combined alignment for {} to {}".format(ref, output_file))


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
