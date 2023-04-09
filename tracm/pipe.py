import os
import sys
import argparse
import time
import shutil
import tempfile
import glob
import re
from collections import defaultdict

from .utils import run_gather, check_positive_int, check_positive_float
from .align import align
from .distance import distance
from .cluster import cluster

def pipe_parser(parser):

    parser.description = "A script to run the full Tracm pipeline."

    io_opts = parser.add_argument_group("Input/output")

    io_opts.add_argument(
        "-i",
        "--input",
        dest="input_file",
        required=True,
        help="path to text file containing input file paths",
        type=os.path.abspath
    )

    io_opts.add_argument(
        "--database",
        dest="database",
        help="path to database signatures",
        type=os.path.abspath,
        required=True
    )

    io_opts.add_argument(
        "-o",
        "--output",
        dest="output_dir",
        required=True,
        help="location of an output directory",
        type=os.path.abspath,
    )

    io_opts.add_argument(
        "--meta",
        dest="metadata",
        default=None,
        help="""Location of metadata in csv format. The first column must include the 
        sequence names and the second column must include sampling dates.""",
        type=os.path.abspath,
    )

    alignment = parser.add_argument_group("Alignment options")

    alignment.add_argument(
        "--minimap_preset",
        dest="minimap_preset",
        help="minimap preset to use - one of 'sr' (default), 'map-ont' or 'map-pb'",
        default="sr",
        type=str,
    )

    pileup = parser.add_argument_group("Pileup options")

    pileup.add_argument(
        "-Q",
        "--min_base_qual",
        dest="min_base_qual",
        help="minimum base quality (default=0)",
        type=int,
        default=0,
    )

    pileup.add_argument(
        "-q",
        "--min_map_qual",
        dest="min_map_qual",
        help="minimum mapping quality (default=0)",
        type=int,
        default=0,
    )

    pileup.add_argument(
        "-l",
        "--min_query_len",
        dest="min_query_len",
        help="minimum query length (default=0)",
        type=int,
        default=0,
    )

    pileup.add_argument(
        "-V",
        "--max_div",
        dest="max_div",
        help="ignore queries with per-base divergence > max_div (default=1)",
        type=float,
        default=1,
    )

    pileup.add_argument(
        "--trim",
        dest="trim",
        help="ignore bases within TRIM-bp from either end of a read (default=0)",
        type=int,
        default=0,
    )

    posterior = parser.add_argument_group("Posterior count estimates")

    posterior.add_argument(
        "--min-cov",
        dest="min_cov",
        default=5,
        help=(
            "Minimum read coverage (default=5)."
        ),
        type=int,
    )

    posterior.add_argument(
        "--error-perc",
        dest="error_threshold",
        default=0.01,
        help=(
            "Threshold to exclude likely erroneous variants prior to"
            + " fitting Dirichlet multinomial model"
        ),
        type=float,
    )

    posterior.add_argument(
        "--either-strand",
        dest="require_both_strands",
        help="turns off the requirement that a variant is supported by both strands",
        action="store_true",
        default=False,
    )

    posterior.add_argument(
        "--keep-all",
        dest="keep_all",
        help="turns on filtering of variants with support below the posterior frequency threshold",
        action="store_true",
        default=False,
    )

    # Distance options
    snpdist = parser.add_argument_group("SNP distance options")

    snpdist.add_argument(
        "-D",
        "--snp_threshold",
        dest="snp_threshold",
        help="Only output those transmission pairs with a SNP distance <= D",
        type=check_positive_int,
        default=2147483647,
    )

    snpdist.add_argument(
        "--filter",
        dest="recomb_filter",
        help="Filter out regions with unusually high SNP distances often caused by HGT",
        action="store_true",
        default=False
    )

    transdist = parser.add_argument_group("Transmission distance options")

    transdist.add_argument(
        "--clock_rate",
        dest="clock_rate",
        help="clock rate as defined in the transcluster paper (SNPs/genome/year) default=1e-3 * 29903",
        type=check_positive_float,
        default=1e-3 * 29903,
    )

    transdist.add_argument(
        "--trans_rate",
        dest="trans_rate",
        help="transmission rate as defined in the transcluster paper (transmissions/year) default=73",
        type=check_positive_float,
        default=73.0,
    )

    transdist.add_argument(
        "-K",
        "--trans_threshold",
        dest="trans_threshold",
        help=(
            "Only outputs those pairs where the most likely number of intermediate hosts <= K"
        ),
        type=check_positive_int,
        default=None,
    )

    transdist.add_argument(
        "--precision",
        dest="precision",
        help=(
            "The precision used to calculate E(K) (default=0.01)."
        ),
        type=check_positive_float,
        default=0.01,
    )

    # Cluster options
    cluster_opts = parser.add_argument_group("Cluster options")

    cluster_opts.add_argument(
        "-c",
        "--cluster_threshold",
        dest="threshold",
        help="Distance threshold. Samples will be grouped together if the distance between them is below this threshold. (default=10)",
        type=float,
        default=10,
    )

    cluster_opts.add_argument(
        "--cluster_distance",
        dest="distance",
        help="The type of transmission distance to use. Can be one of 'SNP' (default), 'direct', 'expectedK'",
        choices=["SNP", "direct", "expectedK"],
        type=str,
        default="SNP"
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

    parser.set_defaults(func=pipe)

    return parser


def pipe(args):

    # get working directory and create temp directory
    # create directory if it isn't present already
    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)
    # make sure trailing forward slash is present
    args.output_dir = os.path.join(args.output_dir, "")
    outputdir = args.output_dir

    # check input file
    prefixes = set()
    with open(args.input_file, 'r') as infile:
        next(infile)
        for line in infile:
            line = line.strip().split()
            if line[0] in prefixes:
                raise ValueError("Repeated file name! " + line[0])
            else:
                prefixes.add(line[0])
            if not os.path.isfile(line[1]):
                raise ValueError("Path does not exist or is not a file! " + line[1])
            if (len(line)>2) and not os.path.isfile(line[2]):
                raise ValueError("Path does not exist or is not a file! " + line[2])

    # generate alignments
    with open(args.input_file, 'r') as infile:
        next(infile)
        for line in infile:
            line = line.strip().split()
            args.input_files = line[1:]
            args.prefix = line[0]
            args.output_dir = outputdir + line[0]
            align(args)

    args.output_dir = outputdir

    # concatenate alignments
    references = defaultdict(list)
    for prefix in prefixes:
        for aln in glob.glob(outputdir + prefix + "/*.fasta"):
            print(aln)
            ref = re.search(r"posterior_counts_ref_(.+?)\.fasta", aln).group(1)
            print(ref)
            references[ref].append(aln)
    
    alignments = []
    for ref in references:
        if len(references[ref])<=1: continue
        combined_aln = outputdir + "combined" + ref
        with open(combined_aln, 'w') as outfile:
            for aln in references[ref]:
                outfile.write(open(aln, "r").read())
        alignments.append(combined_aln)

    # find pairwise distances
    args.output_file = outputdir + "transmission_distances.csv"
    args.msa_files = alignments
    args.msa_db = None
    distance(args)

    # cluster using single linkage
    args.distance_file = outputdir + "transmission_distances.csv"
    args.output_file = outputdir + "transmission_clusters.csv"
    cluster(args)


    return


def main():
    # set up and parse arguments
    parser = argparse.ArgumentParser()
    parser = pipe_parser(parser)
    args = parser.parse_args()

    # run align command
    args.func(args)

    return


if __name__ == "__main__":
    main()
