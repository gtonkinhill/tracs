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
        inputs = [
            (f, os.path.splitext(os.path.basename(f))[0]) for f in args.input_files
        ]

    # find core if necessary
    # filtinputs = find_core(inputs, temp_dir, args.core_thresh, args.n_cpu)
    filtinputs = inputs

    # build zip file to hold database
    with ZipFile(args.dbname + ".zip", "w") as archive:
        # generate sourmash database
        path_to_sourmashdb = build_sourmash_db(
            inputs, temp_dir, ksize=args.ksize, scale=args.scale, n_cpu=1
        )
        archive.write(path_to_sourmashdb, "sourmashDB.sbt.zip")

        # copy genomes into subdirectory and gzip if necessary
        for f, prefix in filtinputs:
            if f.split(".")[-1] == "gz":
                archive.write(f, prefix + ".fasta.gz")
            else:
                newloc = temp_dir + prefix + ".fasta.gz"
                with open(f, "rb") as f_in:
                    with gzip.open(newloc, "wb") as f_out:
                        shutil.copyfileobj(f_in, f_out)
                archive.write(newloc, prefix + ".fasta.gz")
                os.remove(newloc)

        # add summary file for navigation
        with open(temp_dir + "summary.tsv", "w") as outfile:
            for f, prefix in inputs:
                outfile.write(prefix + "," + prefix + ".fasta.gz")
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
