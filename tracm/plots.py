import os
import sys
import argparse
import logging
import plotly.graph_objects as go
import pandas as pd
import plotly.express as px
import plotly.offline as offline
import numpy as np
from collections import Counter
from datetime import datetime
import gzip


def plots_parser(parser):

    parser.description = "Generates plots from a pileup file."

    io_opts = parser.add_argument_group("Input/output")

    io_opts.add_argument(
        "-i",
        "--input",
        dest="input_files",
        required=True,
        help="path to query signature",
        type=os.path.abspath,
        nargs="+",
    )

    io_opts.add_argument(
        "-o",
        "--output",
        dest="output_file",
        required=True,
        help="location of an output file",
        type=os.path.abspath,
    )

    pileup = parser.add_argument_group("Pileup options")
    
    pileup.add_argument(
        "--min-freq",
        dest="min_freq",
        help="minimum frequency to include a variant (default=0.0)",
        type=float,
        default=0.0,
    )

    pileup.add_argument(
        "--either-strand",
        dest="require_both_strands",
        help="turns off the requirement that a variant is supported by both strands",
        action="store_false",
        default=True,
    )

    pileup.add_argument(
        "--contigs",
        dest="contigs",
        default=["All"],
        help="contigs for plotting (default=All)",
        type=str,
        nargs="+",
    )

    parser.set_defaults(func=plots)

    return parser


def plot_heatmap(dates, clusters, outfile):

    cluster_date_count = Counter()
    for i in dates:
        cluster_date_count[(dates[i][0], clusters[i])] += 1

    cluster_counts = []
    dates = []
    cluster = []
    for dc in cluster_date_count:
        cluster_counts.append(cluster_date_count[dc])
        dates.append(datetime.fromisoformat(dc[0]))
        cluster.append(str(dc[1] + 1))

    fig = go.Figure(data=go.Heatmap(
        z=cluster_counts, x=dates, y=cluster, colorscale='Viridis'))

    fig.update_layout(title='Genomes per transmission cluster per day',
                      yaxis_nticks=len(set(cluster)) + 1)

    fig.update_xaxes(title_text='Time')
    fig.update_yaxes(title_text='Transmission Cluster')

    offline.plot(fig, filename=outfile, auto_open=True)

    return



def read_pileup(inputfile, contig_length, require_both_strands=True, keep_contigs='All'):
    npos = {"A": 0, "C": 1, "G": 2, "T": 3}

    # generate count matrix
    all_counts = {}
    for contig in contig_length:
        all_counts[contig] = np.zeros((contig_length[contig], 4), dtype=float)

    with gzip.open(inputfile, 'rt') as infile:
        for i, line in enumerate(infile):
            line = line.strip().split()
            if ("All" in keep_contigs) or (line[0]  in keep_contigs):
                pos = int(line[1]) - 1
                nucs = line[-2].split(",")
                ncounts = line[-1].split(":")[1:]
                counts = np.zeros(4, dtype=float)
                for nuc, c1, c2 in zip(
                    nucs, ncounts[0].split(","), ncounts[1].split(",")
                ):
                    c1 = int(c1)
                    c2 = int(c2)
                    if (nuc not in npos) or (line[2] not in npos):
                        continue
                    if require_both_strands:
                        if (c1 == 0) or (c2 == 0):
                            c1 = c2 = 0
                    counts[npos[nuc]] = c1 + c2
                
                all_counts[line[0]][pos,:] = counts/np.sum(counts)

    return all_counts

def plot_pairwise_comparison(count_file_A, count_file_B, outfile, require_both_strands=True, min_freq=0.01, keep_contigs='All'):

    # Column names
    columns = ['A', 'C', 'G', 'T']

    # preprocess files to check contig names match and to generate count matrices

    contig_length_A = Counter()
    contig_length_B = Counter()

    # count number of entries in file
    logging.info("Counting entries in pileup files...")
    with gzip.open(count_file_A, 'rt') as infile:
        for i, line in enumerate(infile):
            line = line.strip().split()
            if ("All" in keep_contigs) or (line[0]  in keep_contigs):
                if contig_length_A[line[0]] < int(line[1]):
                    contig_length_A[line[0]] = int(line[1])

    with gzip.open(count_file_B, 'rt') as infile:
        for i, line in enumerate(infile):
            line = line.strip().split()
            if ("All" in keep_contigs) or (line[0]  in keep_contigs):
                if contig_length_B[line[0]] < int(line[1]):
                    contig_length_B[line[0]] = int(line[1])

    # check contig names match
    if len(set(contig_length_A.keys()).intersection(set(contig_length_B.keys())))==0:
        raise ValueError("No contig names do not match!")
    
    # find lengths of contigs
    contig_length = Counter()
    for c in contig_length_A:
        if c in contig_length_B:
            contig_length[c] = max(contig_length_A[c], contig_length_B[c])
        else:
            contig_length[c] = contig_length_A[c]
    
    for c in contig_length_B:
        if c not in contig_length_A:
            contig_length[c] = contig_length_B[c]

    # generate frequency matrices
    logging.info("Generating frequency matrices...")
    fA = read_pileup(count_file_A, contig_length, require_both_strands=require_both_strands, keep_contigs=keep_contigs)
    fB = read_pileup(count_file_B, contig_length, require_both_strands=require_both_strands, keep_contigs=keep_contigs)


    logging.info("Computing pairwise comparisons...")
    # generate pairwise comparison
    allmismatches = {}
    for contig in fA:
        allmismatch = ((fA[contig]>0) & (fB[contig]>0)).sum(axis=1) == 0
        allmismatches[contig] = allmismatch & (np.sum(fA[contig], axis=1)>0) & (np.sum(fB[contig], axis=1)>0)

    variablesites = {}
    for contig in fA:
        variablesites[contig] = ((fA[contig] + fB[contig])>0).sum(axis=1) > 1

    matches = {}
    for contig in fA:
        temp = pd.DataFrame((fA[contig]>0) & (fB[contig]>0), columns=columns)
        matches[contig] = np.array(temp.melt(ignore_index=False, var_name='allele', value_name='match')['match'])


    # generate pandas dataframe
    logging.info("Generating pandas dataframe...")

    # Iterate over the contigs
    pdf = pd.DataFrame(columns = ['position', 'allmismatch', 'allele', 'frequency', 'sample', 'contig'])
    for sample, f in zip([count_file_A, count_file_B], [fA, fB]):
        for contig in f:
            # Convert the NumPy matrix to a pandas DataFrame
            df = pd.DataFrame(f[contig], columns=columns)
            df['allmismatch'] = allmismatches[contig]
            df['variable'] = variablesites[contig]

            # Melt the DataFrame into a long format
            long_df = df.melt(ignore_index=False, var_name='allele', value_name='frequency', id_vars=['allmismatch', 'variable'])

            long_df['match'] = matches[contig]
            
            long_df = long_df[long_df['frequency'] > min_freq]
            long_df = long_df[(long_df['frequency'] < 1) | long_df['variable'] | long_df['allmismatch']]
            long_df.reset_index(inplace=True)
            long_df.rename(columns={'index': 'position'}, inplace=True)
            long_df['position'] += 1 # convert to 1-based coordinates
            long_df['sample'] = os.path.basename(sample).replace('.txt.gz', '')
            long_df['contig'] = contig

            pdf = pd.concat([pdf, long_df], ignore_index=True)


    # generate plotly figure
    fig = px.scatter(pdf, x="position", y="frequency", facet_col="contig", facet_row="sample", color="match",
                     color_discrete_map={
                        False: "#d73027",
                        True: "#313695"},
                     hover_data=['allele', 'frequency', 'match'])
    fig.update_layout(yaxis_range=[-0.05,1.05])
    
    if len(keep_contigs)>1:
        fig.update_xaxes(matches=None)

    offline.plot(fig, filename=outfile, auto_open=True)

    return

def plots(args):

    logging.basicConfig(format='%(asctime)s - %(message)s', 
                        datefmt='%d-%b-%y %H:%M:%S',
                        level=logging.INFO)

    plot_pairwise_comparison(args.input_files[0], 
                             args.input_files[1], 
                             args.output_file,
                             require_both_strands=args.require_both_strands,
                             min_freq=args.min_freq,
                             keep_contigs=args.contigs)

    return

def main():
    # set up and parse arguments
    parser = argparse.ArgumentParser()
    parser = plots_parser(parser)
    args = parser.parse_args()

    # run plot command
    args.func(args)

    return


if __name__ == "__main__":
    main()
