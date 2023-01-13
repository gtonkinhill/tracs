import os
import sys
import argparse

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components


def index_count(name):
    if "dict" not in index_count.__dict__:
        index_count.dict = {}
    if "curr" not in index_count.__dict__:
        index_count.curr = 0

    if name not in index_count.dict:
        index_count.dict[name] = index_count.curr
        index_count.curr += 1

    return index_count.dict[name]


def cluster_parser(parser):

    parser.description = "Groups samples into putative transmission clusters using single linkage clustering"

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

    cluster_opts = parser.add_argument_group("Cluster options")

    cluster_opts.add_argument(
        "-c",
        "--threshold",
        dest="threshold",
        help="Distance threshold. Samples will be grouped together if the distance between them is below this threshold.",
        type=float,
        required=True,
    )

    cluster_opts.add_argument(
        "-D",
        "--distance",
        dest="distance",
        help="The type of transmission distance to use. Can be one of 'SNP', 'direct', 'expectedK'",
        choices=["SNP", "direct", "expectedK"],
        type=str,
        default='SNP',
        required=True
    )

    # Other options
    parser.add_argument(
        "--quiet",
        dest="quiet",
        help="turns off some console output",
        action="store_true",
        default=False,
    )

    parser.set_defaults(func=cluster)

    return parser


def cluster(args):

    if args.distance == "SNP":
        col_index = 3
    elif args.distance == "direct":
        col_index = 4
    elif args.distance == "expectedK":
        col_index = 5

    # Load distances
    I = []
    J = []
    indices = {}
    count = 0
    with open(args.distance_file, "r") as infile:
        next(infile)
        for line in infile:
            line = line.strip().split(",")
            i = index_count(line[0])
            j = index_count(line[1])
            if float(line[col_index]) <= args.threshold:
                I.append(i)
                J.append(j)
            count += 1
    
    if count <= 0:
        print("No distances available! Abandoning clustering.")
        return

    # pull out names in order
    names = list(index_count.dict.keys())
    nsamples = len(names)

    if not args.quiet:
        print("Clustering %d samples..." % nsamples)

    # Build sparse graph and find connected components
    G = csr_matrix((np.ones_like(I), (I, J)), shape=(nsamples, nsamples))
    n_components, labels = connected_components(
        csgraph=G, directed=False, return_labels=True
    )

    if not args.quiet:
        print(n_components, " putative transmission clusters found!")

    # write clusters to file
    with open(args.output_file, "w") as outfile:
        outfile.write("sample,cluster\n")
        for i, lab in enumerate(labels):
            outfile.write(names[i] + "," + str(lab) + "\n")

    return


def main():
    # set up and parse arguments
    parser = argparse.ArgumentParser()
    parser = cluster_parser(parser)
    args = parser.parse_args()

    # run cluster command
    args.func(args)

    return


if __name__ == "__main__":
    main()
