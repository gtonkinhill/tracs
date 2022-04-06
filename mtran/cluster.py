import os
import sys
import argparse

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components


from .__init__ import __version__



def index_count(name):
    if "dict" not in index_count.__dict__: index_count.dict = {}
    if "curr" not in index_count.__dict__: index_count.curr = 0

    if name not in index_count.dict:
        index_count.dict[name] = index_count.curr
        index_count.curr += 1
    
    return(index_count.dict[name])


def main():

    parser = argparse.ArgumentParser(
        description="Groups samples into putative transmission clusters using single linkage clustering",
        prog="cluster",
    )

    io_opts = parser.add_argument_group("Input/output")

    io_opts.add_argument(
        "-d",
        "--distances",
        dest="distance_file",
        required=True,
        help="Pairwise distance estimates obtained from running the 'distance' function",
        type=os.path.abspath,
    )

    io_opts.add_argument(
        "-o",
        "--output",
        dest="output_file",
        required=True,
        help="name of the output file to store the resulting cluster assignments",
        type=str,
    )

    cluster = parser.add_argument_group("Cluster options")

    cluster.add_argument(
        "-c",
        "--threshold",
        dest="threshold",
        help="Distance threshold. Samples will be grouped together if the distance between them is below this threshold.",
        type=float,
        required=True,
    )

    cluster.add_argument(
        "-D",
        "--distance",
        dest="distance",
        help="The type of transmission distance to use. Can be one of 'SNP', 'direct', 'expectedK'",
        choices=['SNP', 'direct', 'expectedK'],
        type=str,
        required=True,
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

    parser.add_argument(
        "--version", action="version", version="%(prog)s " + __version__
    )

    args = parser.parse_args()

    match args.distance:
        case 'SNP':
            col_index = 3
        case 'direct':
            col_index = 4
        case 'expectedK':
            col_index = 5

    # Load distances
    I = []
    J = []
    indices = {}
    with open(args.distance_file, "r") as infile:
        next(infile)
        for line in infile:
            line = line.strip().split(",")
            i = index_count(line[0])
            j = index_count(line[1])
            if float(line[col_index]) <= args.threshold:
                I.append(i)
                J.append(j)

    # pull out names in order
    names = list(index_count.dict.keys())
    nsamples = len(names)

    if not args.quiet:
        print("Clustering %d samples..." % nsamples)

    # Build sparse graph and find connected components
    G = csr_matrix((np.ones_like(I), (I, J)), shape=(nsamples, nsamples))
    n_components, labels = connected_components(csgraph=G, directed=False, return_labels=True)

    if not args.quiet:
        print(n_components, " putative transmission clusters found!")

    # write clusters to file
    with open(args.output_file, 'w') as outfile:
        outfile.write('sample,cluster\n')
        for i, lab in enumerate(labels):
            outfile.write(names[i] + ',' + str(lab) + '\n')

    return


if __name__ == "__main__":
    main()
