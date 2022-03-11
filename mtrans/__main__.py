import os, sys
import argparse
from datetime import datetime
import pyfastx
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

from .pairsnp import run_pairsnp
from .transcluster import calculate_trans_prob, lprob_k_given_N
from .plots import plot_heatmap
from .iqtree import run_iqtree_index_cases, run_iqtree_mrca_cases
from .pileup import pileup_dist

from .__init__ import __version__

SECONDS_IN_YEAR = 31556952


def get_options():
    description = 'Runs the pairsnp and transcluster algorithms.'
    parser = argparse.ArgumentParser(description=description,
                                     prog='mtrans')

    io_opts = parser.add_argument_group('Input/output')
    io_opts.add_argument(
        "--msa",
        dest="msa",
        help="Location of fasta formatted multiple sequence alignment")

    io_opts.add_argument("--pileup",
                         dest="pileup",
                         help=("input pileup files"),
                         type=str,
                         nargs='+')

    io_opts.add_argument(
        "--dates",
        dest="metadata",
        required=True,
        help=
        """Location of metadata in csv format. The first column must include the 
        sequence names and the second column must include sampling dates.""")

    io_opts.add_argument("-o",
                         "--out_dir",
                         dest="output_dir",
                         required=True,
                         help="location of an output directory",
                         type=str)

    io_opts.add_argument(
        "--tree",
        dest="tree",
        type=str,
        default=None,
        choices=['index', 'mrca'],
        help=
        ("Toggles the pipeline to build a phylogeny of the initial sequences for"
         +
         " each transmission cluster. Can be based on either the first sequnece in each"
         + " cluster 'index' or the MRCA of each cluster 'mrca"))

    io_opts.add_argument(
        "--heatmap",
        dest="heatmap",
        default=False,
        action='store_true',
        help=("Toggles the pipeline to generate an interactive heatmap"))

    # transcluster options
    transcluster = parser.add_argument_group('transcluster options')
    transcluster.add_argument(
        "-K",
        "--trans_threshold",
        dest="trans_threshold",
        help=
        ("transmission distance threshold - samples are clustered together" +
         " where the implied number of transmissions k is less than or equal to K "
         + " with a probability of P"),
        type=int,
        default=5)

    transcluster.add_argument(
        "-P",
        "--prob_threshold",
        dest="prob_threshold",
        help=
        ("probability threshold - samples are clustered together" +
         " where the implied number of transmissions k is less than or equal to K "
         + " with a probability of at least P"),
        type=float,
        default=0.8)

    transcluster.add_argument(
        "--clock_rate",
        dest="clock_rate",
        help=
        "clock rate as defined in the transcluster paper (SNPs/genome/year) default=1e-3 * 29903",
        type=float,
        default=1e-3 * 29903)

    transcluster.add_argument(
        "--trans_rate",
        dest="trans_rate",
        help=
        "transmission rate as defined in the transcluster paper (transmissions/year) default=73",
        type=float,
        default=73)

    transcluster.add_argument(
        "--save_probs",
        dest="save_probs",
        help="write out transmission probabilites (can be a large file)",
        action='store_true',
        default=False)

    transcluster.add_argument(
        "--log",
        dest="log",
        help="saves probability in log space to improve precision",
        action='store_true',
        default=False)

    # pairsnp options
    pairsnp = parser.add_argument_group('Pairsnp options')
    pairsnp.add_argument(
        "--snp_threshold",
        dest="snp_threshold",
        help="SNP threshold used to sparsify initial SNP distance matrix",
        type=int,
        default=10)

    # Other options
    parser.add_argument("-t",
                        "--threads",
                        dest="n_cpu",
                        help="number of threads to use (default=1)",
                        type=int,
                        default=1)
    parser.add_argument("--quiet",
                        dest="quiet",
                        help="suppress additional output",
                        action='store_true',
                        default=False)
    parser.add_argument('--version',
                        action='version',
                        version='%(prog)s ' + __version__)

    args = parser.parse_args()

    if (args.msa is None) and (args.pileup is None):
        raise ValueError(
            "One option of '--pileup'or '--msa' must be provided!")
    c = 0
    for param in [args.msa, args.pileup]:
        if param is not None: c += 1
    if c > 1:
        raise ValueError(
            "Only one option of '--pileup' or '--msa' can be provided!"
        )

    return (args)


def main():
    args = get_options()

    # create directory and make sure trailing forward slash is present
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    args.output_dir = os.path.join(os.path.abspath(args.output_dir), "")

    samples = []
    sample_to_index = {}
    # get sample names by index from fasta
    if args.msa is not None:
        for i, seq in enumerate(pyfastx.Fasta(args.msa, build_index=False)):
            samples.append(seq[0])
            sample_to_index[seq[0]] = i
    else:
        # get pairwise distances from pileups
        samples, sample_to_index, sparse_dist = pileup_dist(args.pileup,
            max_dist=args.snp_threshold,
            quiet=args.quiet)

    nsamples = len(samples)

    # load metadata
    sample_dates = {}
    with open(args.metadata, 'r') as infile:
        next(infile)
        for line in infile:
            line = line.strip().split(",")
            if line[0] not in sample_to_index:
                print("Missing date in msa!", line)
            else:
                sample_dates[sample_to_index[line[0]]] = (
                    line[1], datetime.fromisoformat(line[1]).timestamp() /
                    SECONDS_IN_YEAR)

    if args.msa is not None:
        # run pairsnp
        snp_dist_file = args.output_dir + "pairsnp_sparse_dist.csv"
        sparse_dist = run_pairsnp(msa=args.msa,
                                  snp_threshold=args.snp_threshold,
                                  outputfile=snp_dist_file,
                                  ncpu=args.n_cpu)


    # run transcluster algorithm
    if args.save_probs:
        prob_out = args.output_dir + "transcluster_probabilities.csv"
    else:
        prob_out = None
    row_ind, col_ind, data = calculate_trans_prob(
        sparse_snp_dist=sparse_dist,
        sample_dates=sample_dates,
        K=args.trans_threshold,
        lamb=args.clock_rate,
        beta=args.trans_rate,
        threshold=args.prob_threshold,
        samplenames=samples,
        outputfile=prob_out,
        log=args.log)

    # generate clusters using single linkage algorithm
    sparse_dist_matrix = csr_matrix((data, (row_ind, col_ind)),
                                    shape=(nsamples, nsamples))
    n_components, labels = connected_components(csgraph=sparse_dist_matrix,
                                                directed=False,
                                                return_labels=True)
    index_to_cluster = {}
    for index, cluster in enumerate(labels):
        index_to_cluster[index] = cluster

    if not args.quiet:
        print("Number of inferred transmission clusters: ", n_components)

    # write results to file
    cluster_output = args.output_dir + "transclusters.csv"
    with open(cluster_output, 'w') as outfile:
        outfile.write("sample,calendar_date,cluster\n")
        for i, sample in enumerate(samples):
            outfile.write(",".join(
                [sample, sample_dates[i][0],
                 str(index_to_cluster[i] + 1)]) + "\n")

    # if requested build a phylogeny of index cases
    if args.tree == 'index':
        run_iqtree_index_cases(msa=args.msa,
                               clusters=index_to_cluster,
                               dates=sample_dates,
                               outdir=args.output_dir,
                               ncpu=args.n_cpu)
    elif args.tree == 'mrca':
        run_iqtree_mrca_cases(msa=args.msa,
                              clusters=index_to_cluster,
                              dates=sample_dates,
                              sparse_snp_dist=sparse_dist,
                              outdir=args.output_dir,
                              ncpu=args.n_cpu)

    # plot results
    if args.heatmap:
        heatmap_output = args.output_dir + "transmission_cluster.html"
        plot_heatmap(sample_dates, index_to_cluster, heatmap_output)

    return


if __name__ == '__main__':
    main()
