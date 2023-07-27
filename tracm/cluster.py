import os
import sys
import argparse
import logging

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
        help="The type of transmission distance to use. Can be one of 'snp', 'filter', 'direct', 'expectedK'",
        choices=["snp", "filter", "direct", "expectedK"],
        type=str,
        required=True,
    )

    # Other options
    parser.add_argument(
        "--loglevel",
        type=str.upper,
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        default="INFO",
        help="Set the logging threshold.",
    )

    parser.set_defaults(func=cluster)

    return parser


def cluster(args):
    # set logging up
    logging.basicConfig(
        level=args.loglevel,
        format="%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    if args.distance == "snp":
        col_index = 3
    elif args.distance == "filter":
        col_index = 6
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
        logging.warning("No distances available! Abandoning clustering.")
        return

    # pull out names in order
    names = list(index_count.dict.keys())
    nsamples = len(names)

    logging.info(f"Clustering {nsamples} samples...")

    # Build sparse graph and find connected components
    G = csr_matrix((np.ones_like(I), (I, J)), shape=(nsamples, nsamples))
    n_components, labels = connected_components(
        csgraph=G, directed=False, return_labels=True
    )

    logging.info(f"{n_components} putative transmission clusters found!")

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
