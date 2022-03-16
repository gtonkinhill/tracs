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

from .utils import sketch_files
from .__init__ import __version__



def build_db(inputs, outputdir, ksize=51, scale=1000, n_cpu=1, quiet=False):
    temp_dir = os.path.join(tempfile.mkdtemp(dir=outputdir), "")


    # sourmash index database
    sourmash_params = 'k=' + str(ksize) + ','
    sourmash_params += 'scaled=' + str(scale) + ','
    sourmash_params += 'noabund'

    # sketch_files(input_files, prefix, outputfile, sourmash_params)
    Parallel(n_jobs=n_cpu)(
        delayed(sketch_files)([f],
                            prefix,
                            temp_dir + prefix + '.sig',
                            sourmash_params)
    for f, prefix in tqdm(inputs,disable=quiet))


    # build database from signatures
    cmd = "sourmash index "
    cmd += outputdir + "sourmashDB.sbt.zip "
    cmd += temp_dir + "*.sig"
    
    if not quiet:
        print("running cmd: " + cmd)

    subprocess.run(cmd, shell=True, check=True)

    # clean up
    shutil.rmtree(temp_dir)

    return (outputdir + "sourmashDB.sbt.zip")



def main():

    parser = argparse.ArgumentParser(
        description="Builds a database for Mtrans",
        prog="buildb",
    )

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

    parser.add_argument(
        "--version", action="version", version="%(prog)s " + __version__
    )

    args = parser.parse_args()


    # create a temporary directory
    wd = os.path.dirname(os.path.realpath(args.dbname))
    temp_dir = os.path.join(tempfile.mkdtemp(dir=wd), "")

    # read in input files/list
    if len(args.input_files)==1:
        with open(args.input_files[0], 'r') as infile:
            inputs = []
            for line in infile:
                line = line.strip().split(',')
                print(line)
                inputs.append((line[1], line[0]))
    else:
        inputs = [(f, os.path.splitext(os.path.basename(f))[0]) for f in args.input_files]

    # build zip file to hold database
    with ZipFile(args.dbname + '.zip', 'w') as archive:
        # generate sourmash database
        path_to_sourmashdb = build_db(inputs, temp_dir, 
                ksize=args.ksize, scale=args.scale, n_cpu=1, quiet=False)
        archive.write(path_to_sourmashdb, "sourmashDB.sbt.zip")

        # copy genomes into subdirectory and gzip if necessary
        for f, prefix in inputs:
            if f.split('.')[-1]=='gz':
                archive.write(f, prefix + '.fasta.gz')
            else:
                newloc = temp_dir + prefix + '.fasta.gz'
                with open(f, 'rb') as f_in:
                    with gzip.open(newloc, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                archive.write(newloc, prefix + '.fasta.gz')
                os.remove(newloc)

        # add summary file for navigation
        with open(temp_dir + 'summary.tsv', 'w') as outfile:
            for f, prefix in inputs:
                outfile.write(prefix + ',' + prefix + '.fasta.gz')
        archive.write(temp_dir + 'summary.tsv', 'summary.tsv')

    # clean up
    shutil.rmtree(temp_dir)

    return


if __name__ == "__main__":
    main()