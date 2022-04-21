import os
import sys
import argparse
import time
import shutil
import gzip
from zipfile import ZipFile
import tempfile
import numpy as np

from .utils import sketch_files, gather
from .pileup import align_and_pileup
from .dirichlet_multinomial import find_dirichlet_priors
from MTRAN import calculate_posteriors


def align_parser(parser):

    parser.description = "Uses sourmash to identify reference genomes within a read set and then aligns reads to each reference using minimap2"

    io_opts = parser.add_argument_group("Input/output")

    io_opts.add_argument(
        "-i",
        "--input",
        dest="input_files",
        required=True,
        help="path to query signature",
        type=os.path.abspath,
        nargs="+",
    )

    io_opts.add_argument(
        "--database",
        dest="database",
        help="path to database signatures",
        type=os.path.abspath,
        default=None,
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
        "-p",
        "--prefix",
        dest="prefix",
        default=None,
        help="prefix to describe the input sample read files",
        type=str,
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
        "--threshold",
        dest="expected_freq_threshold",
        default=None,
        help=(
            "Minimum posterior read frequency threshold."
            + " The default is set that a variant at a "
            + "location is discounted if it is not found "
            + "with a coverage of ~100x"
        ),
        type=float,
    )

    posterior.add_argument(
        "--both-strands",
        dest="require_both_strands",
        help="turns on the requirement that a variant is supported by both strands",
        action="store_true",
        default=False,
    )

    posterior.add_argument(
        "--filter-all",
        dest="keep_all",
        help="turns on filtering of variants with support below the posterior frequency threshold",
        action="store_false",
        default=True,
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

    parser.set_defaults(func=align)

    return parser


def align(args):

    alleles = np.array(["A", "C", "G", "T"])
    iupac_codes = {
        "A": "A",
        "C": "C",
        "G": "G",
        "T": "T",
        "AC": "M",
        "AG": "R",
        "AT": "W",
        "CG": "S",
        "CT": "Y",
        "GT": "K",
        "CGT": "B",
        "AGT": "D",
        "ACT": "H",
        "ACG": "V",
        "ACGT": "N",
    }

    # get working directory and create temp directory
    # create directory if it isn't present already
    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)
    # make sure trailing forward slash is present
    args.output_dir = os.path.join(args.output_dir, "")
    # Create temporary directory
    temp_dir = os.path.join(tempfile.mkdtemp(dir=args.output_dir), "")

    # set prefix to file name if not provided
    if args.prefix is None:
        args.prefix = os.path.splitext(os.path.basename(args.input_files[0]))[0]

    # retrieve sourmash database from zipfile
    with ZipFile(args.database, "r") as archive:
        archive.extract("sourmashDB.sbt.zip", temp_dir)

    # run soursmash 'gather' method
    references = gather(
        input_files=args.input_files,
        databasefile=temp_dir + "sourmashDB.sbt.zip",
        output=args.output_dir + args.prefix + "_sourmash_hits",
        temp_dir=temp_dir,
    )

    # retrieve references and perform alignment
    if len(args.input_files) == 1:
        r1 = args.input_files[0]
        r2 = None
    elif len(args.input_files) == 2:
        r1 = args.input_files[0]
        r2 = args.input_files[1]

    with ZipFile(args.database, "r") as archive:
        print(references)
        for ref in references:
            archive.extract(ref + ".fasta.gz", temp_dir)
            align_and_pileup(
                temp_dir + ref + ".fasta.gz",
                temp_dir,
                args.output_dir + args.prefix + "_ref_" + str(ref),
                r1,
                r2=r2,
                aligner="minimap2",
                minimap_preset=args.minimap_preset,
                minimap_params=None,
                Q = args.min_base_qual, #minimum base quality
                q = args.min_map_qual, #minimum mapping quality
                l = args.min_query_len, #minimum query length
                V = args.max_div, #ignore queries with per-base divergence >FLOAT [1]
                T = args.trim, #ignore bases within INT-bp from either end of a read [0]
                n_cpu=args.n_cpu,
                quiet=args.quiet,
            )

    # add empirical Bayes pseudocounts
    npos = {"A": 0, "C": 1, "G": 2, "T": 3}
    for ref in references:
        all_counts = []
        with open(
            args.output_dir + args.prefix + "_ref_" + str(ref) + "_pileup.txt", "r"
        ) as infile:
            for i, line in enumerate(infile):
                line = line.strip().split()
                nucs = line[-2].split(",")
                ncounts = line[-1].split(":")[1:3]
                counts = np.zeros(4, dtype=float)
                for nuc, c1, c2 in zip(
                    nucs, ncounts[0].split(","), ncounts[1].split(",")
                ):
                    if nuc not in npos:
                        continue
                    if args.require_both_strands:
                        if (c1 == 0) or (c2 == 0):
                            c1 = c2 = 0
                    counts[npos[nuc]] = c1 + c2
                all_counts.append(counts)
        all_counts = np.asarray(all_counts)
        alphas = find_dirichlet_priors(all_counts)
        a0 = np.sum(alphas)

        if args.expected_freq_threshold is None:
            args.expected_freq_threshold = alphas[1] / (np.sum(alphas) + 50)

        if not args.quiet:
            print("Calculating posterior frequency estimates...")
            print(
                "Filtering sites with posterior estimates below frequency threshold:",
                args.expected_freq_threshold,
            )
            if args.keep_all:
                print("Keeping all observed alleles")

        # Calculate posterior frequency estimates and filter out those below the threshold
        all_counts = calculate_posteriors(
            all_counts, alphas, args.keep_all, args.expected_freq_threshold
        )

        # save allele counts to file
        if not args.quiet:
            print("saving to file...")
        with gzip.open(
            args.output_dir
            + args.prefix
            + "_posterior_counts_ref_"
            + str(ref)
            + ".csv.gz",
            "wb",
        ) as outfile:
            np.savetxt(
                outfile,
                all_counts.reshape((1, all_counts.size)),
                delimiter=",",
                newline="\n",
                fmt="%0.5f",
            )
            outfile.write(b"\n")

        # generate fasta output
        with open(
            args.output_dir
            + args.prefix
            + "_posterior_counts_ref_"
            + str(ref)
            + ".fasta",
            "w",
        ) as outfile:
            outfile.write(">" + args.prefix + "_" + str(ref) + "\n")
            for i in range(all_counts.shape[0]):
                outfile.write(iupac_codes["".join(alleles[all_counts[i, :] > 0])])

    shutil.rmtree(temp_dir)

    return


def main():
    # set up and parse arguments
    parser = argparse.ArgumentParser()
    parser = align_parser(parser)
    args = parser.parse_args()

    # run align command
    args.func(args)

    return


if __name__ == "__main__":
    main()
