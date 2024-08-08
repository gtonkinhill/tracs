import os
import sys
import argparse
import logging
import pandas as pd
import numpy as np
from scipy import stats
import scipy.optimize as optimize
from scipy.special import logsumexp


def threshold_parser(parser):

    parser.description = "Estimates transmission thresholds."

    io_opts = parser.add_argument_group("Input/output")

    io_opts.add_argument(
        "--close",
        dest="close_file",
        required=True,
        help="path to csv file with distances between isolates mostly linked by recent transmission",
        type=os.path.abspath
    )

    io_opts.add_argument(
        "--distant",
        dest="distant_file",
        required=True,
        help="path to csv file with distances between isolates not related by recent transmission",
        type=os.path.abspath
    )

    io_opts.add_argument(
        "-o",
        "--output",
        dest="output_file",
        required=True,
        help="location of an output file",
        type=os.path.abspath,
    )

    io_opts.add_argument(
        "--column",
        dest="column",
        default=1,
        help="index of column containing SNP distances (default=1)",
        type=int
    )

    parser.set_defaults(func=threshold)

    return parser

# Define the negative binomial log-likelihood function
def negbinom_ll(params, data):
    r, p = params
    if r <= 0 or p <= 0 or p >= 1:
        return np.inf
    return -np.sum(stats.nbinom.logpmf(data, r, p))

# Define the mixture log-likelihood function
def mixture(params, data, r, p):
    q, lambd = params
    logpmf_poisson = np.log(q) + stats.poisson.logpmf(data, mu=lambd)
    logpmf_nbinom = np.log(1 - q) + stats.nbinom.logpmf(data, r, p)
    return(sum(logsumexp([logpmf_poisson, logpmf_nbinom], axis=0)))

# Define the optimizer function using Nelder-Mead
def optimizer_NM(func, x0):
    result = optimize.minimize(func, x0, method="nelder-mead")
    return result.x

def estimate_thresholds(close_file, distant_file, outfile, column):
    
    logging.info("Loading distances...")
    # read close distances
    df = pd.read_csv(close_file)
    close_distances = df.iloc[:, column].astype(float).values

    # read far distances
    df = pd.read_csv(distant_file)
    far_distances = df.iloc[:, column].astype(float).values

    logging.info("Fitting distribution...")

    # Initial guesses for r and p
    initial_params = np.array([100, 0.5])

    # Fit the data to the negative binomial distribution
    far_fitted_params = optimizer_NM(lambda params: negbinom_ll(params, far_distances), initial_params)
    r, p = far_fitted_params

    # Initial guesses for r and p
    initial_params = np.array([0.5, 1])
    
    # Fit the data to the mixture distribution
    mix_fitted_params = optimizer_NM(lambda params: mixture(params, close_distances, r, p), initial_params)
    q, lambd = mix_fitted_params

    logging.info(f"Fitted parameters - r:{r}, p:{p}, q:{q}, lambda:{lambd}")

    snp_threshold = stats.poisson.ppf(0.95, mu=lambd)*3

    logging.info(f"SNP threshold: {snp_threshold}")
    
    return


def threshold(args):

    logging.basicConfig(format='%(asctime)s - %(message)s', 
                        datefmt='%d-%b-%y %H:%M:%S',
                        level=logging.INFO)

    estimate_thresholds(args.close_file, 
                             args.distant_file, 
                             args.output_file,
                             args.column)

    return

def main():
    # set up and parse arguments
    parser = argparse.ArgumentParser()
    parser = threshold_parser(parser)
    args = parser.parse_args()

    # run plot command
    args.func(args)

    return


if __name__ == "__main__":
    main()
