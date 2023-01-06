import os
import sys
import argparse
from datetime import date
import numpy as np

from TRACM import pairsnp
from .transcluster import calculate_trans_prob
from .utils import check_positive_int, check_positive_float


def distance_parser(parser):

    parser.description = "Estimates pairwise SNP and transmission distances between each pair of samples aligned to the same reference genome."

    io_opts = parser.add_argument_group("Input/output")

    io_opts.add_argument(
        "--msa",
        dest="msa_files",
        required=True,
        help="Input fasta files formatted by the align and merge functions",
        type=os.path.abspath,
        nargs="+",
    )

    io_opts.add_argument(
        "--msa-db",
        dest="msa_db",
        help="A database MSA used to compare each sequence to. By default this is not uses and all pairwise comparisons within each MSA are considered.",
        type=os.path.abspath,
        default=None,
    )

    io_opts.add_argument(
        "--meta",
        dest="metadata",
        default=None,
        help="""Location of metadata in csv format. The first column must include the 
        sequence names and the second column must include sampling dates.""",
        type=os.path.abspath,
    )

    io_opts.add_argument(
        "-o",
        "--output",
        dest="output_file",
        required=True,
        help="name of the output file to store the pairwise distance estimates.",
        type=str,
    )

    snpdist = parser.add_argument_group("SNP distance options")

    snpdist.add_argument(
        "-D",
        "--snp_threshold",
        dest="snp_threshold",
        help="Only output those transmission pairs with a SNP distance <= D",
        type=check_positive_int,
        default=2147483647,
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

    # Other options
    parser.add_argument(
        "-t",
        "--threads",
        dest="n_cpu",
        help="number of threads to use (default=1)",
        type=check_positive_int,
        default=1,
    )

    parser.add_argument(
        "--quiet",
        dest="quiet",
        help="turns off some console output",
        action="store_true",
        default=False,
    )

    parser.set_defaults(func=distance)

    return parser


def distance(args):

    # Load dates
    if not args.quiet:
        print("Loading metadata...")

    if args.metadata is not None:
        dates = {}
        with open(args.metadata, "r") as infile:
            next(infile)
            for line in infile:
                line = line.strip().split(",")
                dates[line[0]] = (line[1], date.fromisoformat(line[1]))

    if not args.quiet:
        print("Estimating transmission distances...")
    
    with open(args.output_file, "w") as outfile:
        outfile.write(
            "sampleA,sampleB,date difference,SNP distance,transmission distance,expected K\n"
        )
        for msa in args.msa_files:
            # Estimate SNP distances
            if not args.quiet:
                print("Calculating pairwise snp distances for %s" % msa)
            # I, J, dist, names
            if args.msa_db is not None:
                msas = [msa, args.msa_db]
            else:
                msas = [msa]
            
            snp_dists = pairsnp(
                fasta=msas, n_threads=args.n_cpu, dist=args.snp_threshold
            )
            names = snp_dists[3]

            # Estimate transmission distances
            if args.metadata is not None:
                if not args.quiet:
                    print("Inferring transmission probabilities for %s" % msa)
                transmission_dists, expectedk, datediff = calculate_trans_prob(
                    snp_dists[:3],
                    sample_dates=dates,
                    K=100,
                    lamb=args.clock_rate,
                    beta=args.trans_rate,
                    samplenames=snp_dists[3],
                    log=False,
                    precision=args.precision
                )

            # Write output
            if not args.quiet:
                print("Saving distances for %s" % msa)
            if args.metadata is not None:
                for i, j, dateD, snpD, expK, tranD in zip(
                    snp_dists[0],
                    snp_dists[1],
                    datediff,
                    snp_dists[2],
                    expectedk,
                    transmission_dists,
                ):
                    if (args.trans_threshold is None) or (args.trans_threshold >= expK):
                        outfile.write(
                            ",".join(
                                [
                                    names[i],
                                    names[j],
                                    str(dateD),
                                    str(int(snpD)),
                                    str(tranD),
                                    str(expK),
                                ]
                            )
                            + "\n"
                        )
            else:
                for i, j, snpD in zip(
                    snp_dists[0],
                    snp_dists[1],
                    snp_dists[2]
                ):
                    outfile.write(
                        ",".join(
                            [
                                names[i],
                                names[j],
                                "NA",
                                str(int(snpD)),
                                "NA",
                                "NA",
                            ]
                        )
                        + "\n"
                    )
    return


def main():
    # set up and parse arguments
    parser = argparse.ArgumentParser()
    parser = distance_parser(parser)
    args = parser.parse_args()

    # run distance command
    args.func(args)

    return


if __name__ == "__main__":
    main()
