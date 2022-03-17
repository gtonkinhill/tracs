import os
import sys
import argparse
import time
import shutil
import gzip
from zipfile import ZipFile
import tempfile

from .__init__ import __version__
from .utils import sketch_files, gather
from .pileup import align_and_pileup



def main():

    parser = argparse.ArgumentParser(
        description="Uses sourmash to identify reference genomes within a read set and then aligns reads to each reference using minimap2",
        prog="align",
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
        "--database",
        dest="database",
        help="path to database signatures",
        type=os.path.abspath,
        default=None,
    )

    parser.add_argument(
        "-o",
        "--output",
        dest="output_dir",
        required=True,
        help="location of an output directory",
        type=os.path.abspath,
    )

    parser.add_argument(
        "-p",
        "--prefix",
        dest="prefix",
        default=None,
        help="prefix to describe the input sample read files",
        type=str,
    )

    parser.add_argument("--minimap_preset",
                         dest="minimap_preset",
                         help="minimap preset to use - one of 'sr' (default), 'map-ont' or 'map-pb'",
                         default='sr',
                         type=str)

    parser.add_argument("-t",
                        "--threads",
                        dest="n_cpu",
                        help="number of threads to use (default=1)",
                        type=int,
                        default=1)

    parser.add_argument("--quiet",
                        dest="quiet",
                        help="turns off some console output",
                        action='store_true',
                        default=False)

    parser.add_argument('--version',
                        action='version',
                        version='%(prog)s ' + __version__)

    args = parser.parse_args()

    # get working directory and create temp directory
    # create directory if it isn't present already
    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)
    # make sure trailing forward slash is present
    args.output_dir = os.path.join(args.output_dir, "")
    # Create temporary directory
    temp_dir = os.path.join(tempfile.mkdtemp(dir=args.output_dir), "")

    # retrieve sourmash database from zipfile
    with ZipFile(args.database, 'r') as archive:
        archive.extract("sourmashDB.sbt.zip", temp_dir)

    # run soursmash 'gather' method
    references = gather(
        input_files=args.input_files,
        databasefile=temp_dir + "sourmashDB.sbt.zip",
        output=args.output_dir + "sourmash_hits",
        temp_dir=temp_dir
    )

    # retrieve references and perform alignment
    if len(args.input_files)==1:
        r1 = args.input_files[0]
        r2 = None
    elif len(args.input_files)==2:
        r1 = args.input_files[0]
        r2 = args.input_files[1]

    with ZipFile(args.database, 'r') as archive:
        print(references)
        for ref in references:
            archive.extract(ref + '.fasta.gz', temp_dir)
            align_and_pileup(temp_dir + ref + '.fasta.gz',
                temp_dir,
                args.output_dir + "ref_" + str(ref),
                r1, 
                r2=r2,
                aligner='minimap2',
                minimap_preset=args.minimap_preset,
                minimap_params=None,
                n_cpu=args.n_cpu,
                quiet=args.quiet)


    shutil.rmtree(temp_dir)

    return


if __name__ == "__main__":
    main()