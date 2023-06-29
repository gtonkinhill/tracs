import os, sys
import argparse
from .__init__ import __version__

from .build_db import build_db_parser
from .distance import distance_parser
from .align import align_parser
from .cluster import cluster_parser
from .combine import combine_parser
from .pipe import pipe_parser


def main():

    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help="select a subcommand", dest="command")

    parser.add_argument(
        "--version", action="version", version="%(prog)s " + __version__
    )

    # add subcommands
    align_subparser = subparsers.add_parser("align")
    align_subparser = align_parser(align_subparser)

    combine_subparser = subparsers.add_parser("combine")
    combine_subparser = combine_parser(combine_subparser)

    distance_subparser = subparsers.add_parser("distance")
    distance_subparser = distance_parser(distance_subparser)

    cluster_subparser = subparsers.add_parser("cluster")
    cluster_subparser = cluster_parser(cluster_subparser)

    build_db_subparser = subparsers.add_parser("build-db")
    build_db_subparser = build_db_parser(build_db_subparser)

    pipe_subparser = subparsers.add_parser("pipe")
    pipe_subparser = pipe_parser(pipe_subparser)

    # parse arguments and run function
    args = parser.parse_args()
    args.func(args)

    return


if __name__ == "__main__":
    main()
