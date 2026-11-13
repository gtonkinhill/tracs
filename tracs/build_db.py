import os
import sys
import argparse
import logging
import subprocess
import tempfile
import shutil
import gzip
from zipfile import ZipFile
from tqdm import tqdm
from joblib import Parallel, delayed
import pyfastx as fx
import numpy as np

from .utils import run_sketch


def build_db_parser(parser):
    parser.description = "Builds a database for tracs"

    parser.add_argument(
        "-i",
        "--input",
        dest="input_files",
        required=True,
        help="path to genome fasta files (one per reference genome).",
        type=os.path.abspath,
        nargs="+",
    )

    parser.add_argument(
        "-o",
        "--output",
        dest="dbname",
        required=True,
        help="name of the database file",
        type=os.path.abspath,
    )

    # parser.add_argument(
    #     "--core-thresh",
    #     dest="core_thresh",
    #     help="the conservation required to call a region of a reference genome as 'core' (default=0.95)",
    #     default=0.95,
    #     type=float,
    # )

    parser.add_argument(
        "--ksize",
        dest="ksize",
        help="the kmer length used in sourmash (default=51)",
        default=51,
        type=int,
    )

    parser.add_argument(
        "--scale",
        dest="scale",
        help="the scale used in sourmash (default=100)",
        default=1000,
        type=int,
    )

    parser.add_argument(
        "--screen_method",
        dest="screen_method",
        help="screening database to include in the bundle (default=sourmash)",
        choices=["sourmash", "sylph"],
        default="sourmash",
        type=str,
    )

    parser.add_argument(
        "--sylph_c",
        dest="sylph_c",
        help="sylph subsampling rate; larger values give a smaller database (default=200)",
        default=200,
        type=int,
    )

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

    parser.set_defaults(func=build_db)

    return parser


def build_sourmash_db(inputs, outputdir, ksize=51, scale=1000, n_cpu=1):
    temp_dir = os.path.join(tempfile.mkdtemp(dir=outputdir), "")

    # sketch_files(input_files, prefix, outputfile, sourmash_params)
    Parallel(n_jobs=n_cpu)(
        delayed(run_sketch)([f], prefix, temp_dir + prefix + ".sig", ksize, scale)
        for f, prefix in tqdm(inputs)
    )

    # build database from signatures
    cmd = "sourmash index "
    cmd += outputdir + "sourmashDB.sbt.zip "
    cmd += temp_dir + "*.sig"

    logging.info(f"running cmd: {cmd}")

    subprocess.run(cmd, shell=True, check=True)

    # clean up
    shutil.rmtree(temp_dir)

    return outputdir + "sourmashDB.sbt.zip"


def build_sylph_db(genome_files, outputdir, c=200, n_cpu=1):
    # genome files must already be named <prefix>.fasta.gz so that the
    # Genome_file reported by `sylph profile` maps back to the bundled fasta
    out_prefix = os.path.join(outputdir, "sylphDB")
    cmd = ["sylph", "sketch", "-c", str(c), "-t", str(n_cpu), "-o", out_prefix, "-g"]
    cmd += list(genome_files)

    logging.info(f"running cmd: {subprocess.list2cmdline(cmd)}")
    subprocess.run(cmd, check=True)

    return out_prefix + ".syldb"


def genome_prefix(f):
    # strip .gz before the fasta extension so X.fasta.gz -> X
    name = os.path.basename(f)
    if name.endswith(".gz"):
        name = name[:-3]
    return os.path.splitext(name)[0]


def build_db(args):
    # set logging up
    logging.basicConfig(
        level=args.loglevel,
        format="%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    # create a temporary directory
    wd = os.path.dirname(os.path.realpath(args.dbname))
    temp_dir = os.path.join(tempfile.mkdtemp(dir=wd), "")

    # read in input files/list
    if len(args.input_files) == 1:
        with open(args.input_files[0], "r") as infile:
            inputs = []
            for line in infile:
                line = line.strip().split(",")
                inputs.append((line[1], line[0]))
    else:
        inputs = [(f, genome_prefix(f)) for f in args.input_files]

    # find core if necessary
    # filtinputs = find_core(inputs, temp_dir, args.core_thresh, args.n_cpu)
    filtinputs = inputs

    # stage genomes as <prefix>.fasta.gz, gzipping if necessary
    genome_dir = os.path.join(temp_dir, "genomes", "")
    os.mkdir(genome_dir)
    staged = []
    for f, prefix in filtinputs:
        newloc = genome_dir + prefix + ".fasta.gz"
        if f.split(".")[-1] == "gz":
            shutil.copyfile(f, newloc)
        else:
            with open(f, "rb") as f_in:
                with gzip.open(newloc, "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
        staged.append((newloc, prefix))

    # build zip file to hold database
    with ZipFile(args.dbname + ".zip", "w") as archive:
        if args.screen_method == "sourmash":
            path_to_sourmashdb = build_sourmash_db(
                inputs, temp_dir, ksize=args.ksize, scale=args.scale, n_cpu=1
            )
            archive.write(path_to_sourmashdb, "sourmashDB.sbt.zip")

        elif args.screen_method == "sylph":
            path_to_sylphdb = build_sylph_db(
                [f for f, _ in staged], temp_dir, c=args.sylph_c, n_cpu=args.n_cpu
            )
            archive.write(path_to_sylphdb, "sylphDB.syldb")

        # add genomes to the bundle
        for f, prefix in staged:
            archive.write(f, prefix + ".fasta.gz")

        # add summary file for navigation
        with open(temp_dir + "summary.tsv", "w") as outfile:
            for f, prefix in inputs:
                outfile.write(prefix + "," + prefix + ".fasta.gz\n")
        archive.write(temp_dir + "summary.tsv", "summary.tsv")

    # clean up
    shutil.rmtree(temp_dir)

    return


def main():
    # set up and parse arguments
    parser = argparse.ArgumentParser()
    parser = build_db_parser(parser)
    args = parser.parse_args()

    # run build_db command
    args.func(args)

    return


if __name__ == "__main__":
    main()
