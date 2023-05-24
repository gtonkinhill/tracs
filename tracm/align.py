import os
import sys
import argparse
import logging
import time
import shutil
import gzip
from zipfile import ZipFile
import tempfile
import numpy as np
from scipy.stats import norm, poisson
import pyfastx as fx
import ncbi_genome_download as ngd
import glob

from .utils import run_gather, generate_reads
from .pileup import align_and_pileup, align_and_pileup_composite
from .dirichlet_multinomial import find_dirichlet_priors
from TRACM import calculate_posteriors

from collections import Counter


os.environ["OPENBLAS_NUM_THREADS"] = "1"


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
        "--refseqs",
        dest="refseqs",
        help="path to reference fasta files",
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
        "--min-cov",
        dest="min_cov",
        default=5,
        help=("Minimum read coverage (default=5)."),
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
        action="store_false",
        default=True,
    )

    posterior.add_argument(
        "--keep-all",
        dest="keep_all",
        help="turns on filtering of variants with support below the posterior frequency threshold",
        action="store_true",
        default=False,
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
        "--loglevel",
        type=str.upper,
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        default="INFO",
        help="Set the logging threshold.",
    )

    parser.set_defaults(func=align)

    return parser


def download_ref(ref, outputdir):
    r = ngd.download(
        groups="bacteria",
        section="genbank",
        file_formats="fasta",
        flat_output=True,
        output=outputdir,
        assembly_accessions=ref,
    )
    if r != 0:
        # try refseq
        r = ngd.download(
            groups="bacteria",
            section="refseq",
            file_formats="fasta",
            flat_output=True,
            output=outputdir,
            assembly_accessions=ref,
        )

    if r != 0:
        raise ValueError("Could not download reference for: ", ref)

    refpath = glob.glob(outputdir + "*fna.gz")[0]

    return refpath


def align(args):

    # set logging up
    logging.basicConfig(
        level=args.loglevel,
        format="%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    # check inputs
    if args.database is None and args.refseqs is None:
        logging.error("Must provide either a database or reference sequences!")
        sys.exit(1)

    if args.database is not None:
        if ".zip" not in args.database:
            logging.error("Database must be a zip file!")
            sys.exit(1)
    
    single_ref = False
    if args.refseqs is not None and args.database is None:
        if ".fna" not in args.refseqs and ".fasta" not in args.refseqs:
            logging.error("Reference sequences must be a fasta file if not using a database!")
            sys.exit(1)
        else:
            single_ref = True
            references = [os.path.splitext(os.path.basename(args.refseqs))[0]]
            ref_locs = {references[0]: args.refseqs}



    b = np.array(
        [
            [0, 0, 0, 0],  # X
            [1, 0, 0, 0],  # A
            [0, 1, 0, 0],  # C
            [0, 0, 1, 0],  # G
            [0, 0, 0, 1],  # T
            [1, 1, 0, 0],  # AC
            [1, 0, 1, 0],  # AG
            [1, 0, 0, 1],  # AT
            [0, 1, 1, 0],  # CG
            [0, 1, 0, 1],  # CT
            [0, 0, 1, 1],  # GT
            [0, 1, 1, 1],  # CGT
            [1, 0, 1, 1],  # AGT
            [1, 1, 0, 1],  # ACT
            [1, 1, 1, 0],  # ACG
            [1, 1, 1, 1],  # ACGT
        ]
    )
    iupac_codes = np.chararray(b.shape[0])
    iupac_codes[np.packbits(b, axis=1, bitorder="little").flatten()] = [
        b"X",
        b"A",
        b"C",
        b"G",
        b"T",
        b"M",
        b"R",
        b"W",
        b"S",
        b"Y",
        b"K",
        b"B",
        b"D",
        b"H",
        b"V",
        b"N",
    ]

    # get working directory and create temp directory
    # create directory if it isn't present already
    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)
    # make sure trailing forward slash is present
    args.output_dir = os.path.join(args.output_dir, "")
    # Create temporary directory
    temp_dir = os.path.join(tempfile.mkdtemp(dir=args.output_dir), "")
    # temp_dir = args.output_dir + "temp/"

    # set prefix to file name if not provided
    if args.prefix is None:
        args.prefix = os.path.splitext(os.path.basename(args.input_files[0]))[0]

    if not single_ref:
        # retrieve sourmash database from zipfile
        if ".sbt.zip" in args.database:
            smdb = args.database
        else:
            with ZipFile(args.database, "r") as archive:
                archive.extract("sourmashDB.sbt.zip", temp_dir)
                smdb = temp_dir + "sourmashDB.sbt.zip"

        # run soursmash 'gather' method
        references = run_gather(
            input_files=args.input_files,
            databasefile=smdb,
            output=args.output_dir + args.prefix + "_sourmash_hits",
            temp_dir=temp_dir,
        )

        ref_locs = {}
        if ".sbt.zip" in args.database:
            logging.warning(
                "No references provided. Tracm will attempt to download references from Genbank"
            )
            if not os.path.exists(args.output_dir + "genbank_references"):
                os.mkdir(args.output_dir + "genbank_references")

            # attempt to download references
            references = [r.split()[0].strip('"') for r in references]

            logging.debug('%s', references)

            for ref in references:
                temprefdir = args.output_dir + "genbank_references/" + ref + "/"
                if not os.path.exists(temprefdir):
                    os.mkdir(temprefdir)
                    ref_locs[ref] = download_ref(ref, temprefdir)
                else:
                    ref_locs[ref] = glob.glob(temprefdir + "*.fna.gz")[0]
        else:
            with ZipFile(args.database, "r") as archive:
                for ref in references:
                    archive.extract(ref + ".fasta.gz", temp_dir)
                    ref_locs[ref] = temp_dir + ref + ".fasta.gz"

    # retrieve references and perform alignment
    if len(args.input_files) == 1:

        logging.debug('%s', os.path.splitext(args.input_files[0])[1])

        if os.path.splitext(args.input_files[0])[1] in [".fasta", ".fa"]:
            # shred fasta to enable alignment step
            r1 = temp_dir + "simulated_" + os.path.basename(args.input_files[0]) + ".gz"
            generate_reads(args.input_files[0], r1)
        else:
            r1 = args.input_files[0]
        r2 = None
    elif len(args.input_files) == 2:
        r1 = args.input_files[0]
        r2 = args.input_files[1]

    composite = False
    if composite:
        align_and_pileup_composite(
            ref_locs,
            temp_dir,
            args.output_dir + args.prefix,
            r1,
            r2=r2,
            aligner="minimap2",
            minimap_preset=args.minimap_preset,
            minimap_params=None,
            Q=args.min_base_qual,  # minimum base quality
            q=args.min_map_qual,  # minimum mapping quality
            l=args.min_query_len,  # minimum query length
            V=args.max_div,  # ignore queries with per-base divergence >FLOAT [1]
            T=args.trim,  # ignore bases within INT-bp from either end of a read [0]
            n_cpu=args.n_cpu
        )
    else:
        for ref in references:
            align_and_pileup(
                ref_locs[ref],
                temp_dir,
                args.output_dir + args.prefix + "_ref_" + str(ref),
                r1,
                r2=r2,
                aligner="minimap2",
                minimap_preset=args.minimap_preset,
                minimap_params=None,
                Q=args.min_base_qual,  # minimum base quality
                q=args.min_map_qual,  # minimum mapping quality
                l=args.min_query_len,  # minimum query length
                V=1,  # ignore queries with per-base divergence >FLOAT [1]
                T=args.trim,  # ignore bases within INT-bp from either end of a read [0]
                max_div=args.max_div,
                n_cpu=args.n_cpu
            )

    # add empirical Bayes pseudocounts
    npos = {"A": 0, "C": 1, "G": 2, "T": 3}
    for ref in references:
        logging.info(f"Analysing reference: {ref}")

        all_counts = {}
        for name, seq in fx.Fasta(ref_locs[ref], build_index=False):
            all_counts[name] = np.zeros((len(seq), 4), dtype=float)

        with gzip.open(
            args.output_dir + args.prefix + "_ref_" + str(ref) + "_pileup.txt.gz", "rt"
        ) as infile:
            for i, line in enumerate(infile):
                line = line.strip().split()
                contig = line[0]
                pos = int(line[1]) - 1
                nucs = line[-2].split(",")
                ncounts = line[-1].split(":")[1:]
                counts = np.zeros(4, dtype=float)
                for nuc, c1, c2 in zip(
                    nucs, ncounts[0].split(","), ncounts[1].split(",")
                ):
                    c1 = int(c1)
                    c2 = int(c2)
                    if (nuc not in npos) or (line[2] not in npos):
                        continue
                    if args.require_both_strands:
                        if (c1 == 0) or (c2 == 0):
                            c1 = c2 = 0
                    counts[npos[nuc]] = c1 + c2
                all_counts[contig][pos, :] = counts
        all_counts = np.concatenate(list(all_counts.values()))

        # minimum coverage
        rs = np.sum(all_counts, 1)
        nz_cov = np.sum(all_counts[rs > 0,], 1)
        total_cov = np.sum(rs > 0) / all_counts.shape[0]
        median_cov = np.median(nz_cov)

        # calculate minimum frequency threshold
        expected_freq_threshold = max(args.min_cov / median_cov, args.error_threshold)
        total_cov_min_threshold = np.sum(rs >= args.min_cov) / all_counts.shape[0]

        logging.info(f"Fraction of genome with read coverage: {total_cov}")
        logging.info(
            f"Fraction of genome with read coverage >= {args.min_cov}: {total_cov_min_threshold}"
        )
        logging.info(f"Median non-zero coverage: {median_cov}")

        if total_cov_min_threshold < 0.25:
            logging.info(
                f"Skipping reference: {ref} as less than 25% of the genome has sufficient read coverage."
            )
            continue

        alphas = find_dirichlet_priors(
            all_counts, method="FPI", error_filt_threshold=args.error_threshold
        )

        if expected_freq_threshold <= alphas[1] / (median_cov + np.sum(alphas)):
            expected_freq_threshold = alphas[1] / (median_cov + np.sum(alphas)) + 0.01
            logging.warning(
                "WARNING: Frequency threshold is set too low! The majority of the genome will be called as ambiguous."
            )
            logging.warning(
                "WARNING: The threshold has been automatically increased to: %s",
                expected_freq_threshold
            )

        # calculate coverage threshold to handle differences in gene presence and absence
        cov_filter_threshold = 50
        if (median_cov > cov_filter_threshold) and (
            alphas[1] / np.sum(alphas) > expected_freq_threshold
        ):
            bad_cov_lower_bound = alphas[1] / expected_freq_threshold - np.sum(alphas)

            lq = np.quantile(nz_cov, [0.25, 0.5])
            bad_cov_upper_bound = lq[0] - 1.5 * (lq[1] - lq[0])

            if bad_cov_lower_bound < bad_cov_upper_bound:
                logging.info(f"Lower coverage bound: {bad_cov_lower_bound}")
                logging.info(f"Upper coverage bound: {bad_cov_upper_bound}")

        logging.info(f"Using frequency threshold: {expected_freq_threshold}")

        logging.info("Calculating posterior frequency estimates...")
        logging.info(
            f"Filtering sites with posterior estimates below frequency threshold: {expected_freq_threshold}"
        )
        if args.keep_all:
            logging.info("Keeping all observed alleles")

        # Calculate posterior frequency estimates and filter out those below the threshold
        all_counts = calculate_posteriors(
            all_counts, alphas, args.keep_all, expected_freq_threshold
        )

        # save allele counts to file
        logging.info("saving to file...")
        
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
                all_counts,
                delimiter=",",
                newline="\n",
                fmt="%0.5f",
            )
            outfile.write(b"\n")

        # apply coverage filter
        if (median_cov > cov_filter_threshold) and (
            alphas[1] / np.sum(alphas) > expected_freq_threshold
        ):
            logging.info(
                "Fraction of genome filtered by coverage: %s",
                np.sum((rs < bad_cov_upper_bound) & (rs > bad_cov_lower_bound))
                / len(rs)
            )
            if bad_cov_upper_bound > bad_cov_lower_bound:
                all_counts[
                    (rs <= bad_cov_upper_bound) & (rs >= bad_cov_lower_bound),
                ] = 1
        all_counts[rs < args.min_cov,] = 1

        # generate fasta outputs
        sequence = (
            iupac_codes[
                np.packbits(all_counts > 0, axis=1, bitorder="little").flatten()
            ]
            .tobytes()
            .decode("utf-8")
        )
        allelecount = Counter(sequence)
        logging.info(f"allelecount: {allelecount}")

        if sequence.count("N") / (float(len(sequence))) > 0.25:
            logging.info(
                f"Skipping reference: {ref} as greater than 25% of the genome has completely ambiguous (N) base calls!"
            )
            continue

        with open(
            args.output_dir
            + args.prefix
            + "_posterior_counts_ref_"
            + str(ref)
            + ".fasta",
            "w",
        ) as outfile:
            outfile.write(">" + args.prefix + "_" + str(ref) + "\n")
            outfile.write(sequence + "\n")

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
