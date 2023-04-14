import os
import sys
import re
import argparse
import glob
import gzip
from collections import defaultdict

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
        for aln in glob.iglob(os.path.join(directory, "*posterior_counts_ref_*.fasta")):
            ref = find_ref(aln)
            sample = os.path.basename(os.path.normpath(aln))
            alignments[ref].append((sample, aln))

    # write out as gzipped multifasta files
    for ref, alns in alignments.items():
        if not args.quiet:
            print("Writing combined alignment for {} to {}".format(ref, args.output_dir + ref + '_combined.fasta.gz'))
        with gzip.open(args.output_dir + ref + '_combined.fasta.gz', "wt") as fasta_file:
            for seq in alns:
                fasta_file.write(f'>{seq[0]}\n{seq[1]}\n')

    return


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
