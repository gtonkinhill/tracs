import os
import sys
import argparse
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

    parser.description = "Builds a database for tracm"

    parser.add_argument(
        "-i",
        "--input",
        dest="input_files",
        required=True,
        help="path to query signature",
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

    parser.add_argument(
        "--core-thresh",
        dest="core_thresh",
        help="the conservation required to call a region of a reference genome as 'core' (default=0.95)",
        default=0.95,
        type=float,
    )

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
        "--quiet",
        dest="quiet",
        help="turns off some console output",
        action="store_true",
        default=False,
    )

    parser.set_defaults(func=build_db)

    return parser


def find_core(genomes, outputdir, core_thresh, ncpu, quiet=False):

    coverage = {}
    ngenomes = len(genomes)
    ninputs = []
    for name, seq in fx.Fasta(genomes[0][0], build_index=False):
        coverage[name] = np.zeros(len(seq))

    for genome in genomes[1:]:
        temp_file = tempfile.NamedTemporaryFile(delete=False, dir=outputdir)
        temp_file.close()

        cmd = "minimap2 --secondary=no -cx asm10 "
        cmd += " -t " + str(ncpu)
        cmd += " " + genomes[0][0]
        cmd += " " + genome[0]
        cmd += ' > ' + temp_file.name + ' 2> /dev/null '

        print("initial alignment: " + genome[1])
        subprocess.run(cmd, shell=True, check=True)


        with open(temp_file.name, 'r') as infile:
            for line in infile:
                line = line.strip().split()
                # qname = line[0]
                # qstart = int(line[2])
                # qend = int(line[3])
                tname = line[5]
                tstart = int(line[7])
                tend = int(line[8])
                if tend < tstart:
                    raise ValueError("Alignmet is reversed!")

                coverage[tname][tstart:tend] += 1
        
        os.remove(temp_file.name)
    
    filt_ref = outputdir + genomes[0][1] + "_core.fasta"
    core_size = 0
    with open(filt_ref, 'w') as outfile:
        for name, seq in fx.Fasta(genomes[0][0], build_index=False):
            seq = np.array(list(seq))
            seq[coverage[name]/float(ngenomes-1) < core_thresh] = 'N'
            outfile.write(">" + name + '\n' + ''.join(seq) + '\n')
            core_size += np.sum(coverage[name]/float(ngenomes-1) >= core_thresh)

    for genome in genomes[1:]:
        temp_file = tempfile.NamedTemporaryFile(delete=False, dir=outputdir)
        temp_file.close()

        cmd = "minimap2 --secondary=no -cx asm10 "
        cmd += " -t " + str(ncpu)
        cmd += " " + genome[0]
        cmd += " " + filt_ref
        cmd += ' > ' + temp_file.name + ' 2> /dev/null '

        print("secondary alignment: " + genome[1])
        subprocess.run(cmd, shell=True, check=True)

        nref = {}
        for name, seq in fx.Fasta(genome[0], build_index=False):
            nref[name] = np.zeros(len(seq))

        with open(temp_file.name, 'r') as infile:
            for line in infile:
                line = line.strip().split()
                # qname = line[0]
                # qstart = int(line[2])
                # qend = int(line[3])
                tname = line[5]
                tstart = int(line[7])
                tend = int(line[8])
                if tend < tstart:
                    raise ValueError("Alignmet is reversed!")

                nref[tname][tstart:tend] = 1

        with open(outputdir + genome[1] + "_core.fasta", 'w') as outfile:
            for name, seq in fx.Fasta(genome[0], build_index=False):
                seq = np.array(list(seq))
                seq[nref[name] < 1] = 'N'
                outfile.write(">" + name + '\n' + ''.join(seq) + '\n')

    ninputs = [(outputdir + g[1] + "_core.fasta", g[1]) for g in genomes]

    print("core genome size: " + str(core_size) + "nt")

    return ninputs


def build_sourmash_db(inputs, outputdir, ksize=51, scale=1000, n_cpu=1, quiet=False):
    temp_dir = os.path.join(tempfile.mkdtemp(dir=outputdir), "")


    # sketch_files(input_files, prefix, outputfile, sourmash_params)
    Parallel(n_jobs=n_cpu)(
        delayed(run_sketch)([f], prefix, temp_dir + prefix + ".sig", ksize, scale)
        for f, prefix in tqdm(inputs, disable=quiet)
    )

    # build database from signatures
    cmd = "sourmash index "
    cmd += outputdir + "sourmashDB.sbt.zip "
    cmd += temp_dir + "*.sig"

    if not quiet:
        print("running cmd: " + cmd)

    subprocess.run(cmd, shell=True, check=True)

    # clean up
    shutil.rmtree(temp_dir)

    return outputdir + "sourmashDB.sbt.zip"


def build_db(args):

    # create a temporary directory
    wd = os.path.dirname(os.path.realpath(args.dbname))
    temp_dir = os.path.join(tempfile.mkdtemp(dir=wd), "")

    # read in input files/list
    if len(args.input_files) == 1:
        with open(args.input_files[0], "r") as infile:
            inputs = []
            for line in infile:
                line = line.strip().split(",")
                print(line)
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
            inputs, temp_dir, ksize=args.ksize, scale=args.scale, n_cpu=1, quiet=False
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
