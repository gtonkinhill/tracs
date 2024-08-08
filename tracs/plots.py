import os
import sys
import argparse
import logging
import plotly.graph_objects as go
import pandas as pd
import plotly.express as px
import plotly.offline as offline
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list
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
        "-p",
        "--prefix",
        dest="output_file",
        required=True,
        help="prefix of output file",
        type=os.path.abspath,
    )

    io_opts.add_argument(
        "--type",
        dest="plot_type",
        required=True,
        help="Type of plot (scatter, line, heatmap)",
        choices=["scatter", "line", "heatmap"],
        type=str,
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

    distance = parser.add_argument_group("Transmission distance options")
    
    distance.add_argument(
            "--column-name",
            dest="column_name",
            help="Column name in distance matrix to use (default='SNP distance')",
            type=str,
            default="SNP distance",
        )

    distance.add_argument(
            "--threshold",
            dest="threshold",
            help="threshold to filter transmission distances (default=None)",
            type=float,
            default=None,
        )

    plot = parser.add_argument_group("Plot options")
    
    plot.add_argument(
        "--alpha",
        dest="alpha",
        help="alpha value for plotting (default=0.1)",
        type=float,
        default=0.1,
    )

    plot.add_argument(
        "--height",
        dest="height",
        help="height value for plotting (default=7)",
        type=float,
        default=7,
    )

    plot.add_argument(
        "--width",
        dest="width",
        help="width value for plotting (default=10)",
        type=float,
        default=10,
    )

    parser.set_defaults(func=plots)

    return parser

def open_file(file_path):
    if file_path.endswith('.gz'):
        return gzip.open(file_path, 'rt')
    else:
        return open(file_path, 'r')

def plot_heatmap(distance_file, outfile, column='SNP distance', threshold=None, height=7, width=10):

    df = pd.read_csv(distance_file)

    # filter by threshold
    if threshold is not None:
        df = df[df[column] <= threshold]

    # Pivot the DataFrame to create a matrix for the heatmap
    heatmap_data = df.pivot(index="sampleA", columns="sampleB", values="SNP distance")
    samples = sorted(set(df['sampleA']).union(set(df['sampleB'])))

    # Fill NaN values with a large number to indicate no similarity
    heatmap_data = heatmap_data.combine_first(heatmap_data.T)

    temp_data = heatmap_data.fillna(heatmap_data.max().max() + 100)

    # Perform hierarchical clustering using single linkage
    linkage_matrix = linkage(temp_data, method='single')
    ordered_indices = leaves_list(linkage_matrix)
    ordered_samples = [samples[i] for i in ordered_indices]

    # Reorder the DataFrame based on the clustering results
    ordered_heatmap_data = heatmap_data.reindex(index=ordered_samples, columns=ordered_samples)

    # Plot the heatmap
    heatmap_matrix = ordered_heatmap_data.values    
    fig, ax = plt.subplots(figsize=(width, height))
    cax = ax.matshow(heatmap_matrix, cmap='viridis')

    # Add color bar
    cbar = fig.colorbar(cax)
    cbar.set_label(column)

    # Set the tick labels
    ax.set_xticks(np.arange(len(ordered_samples)))
    ax.set_yticks(np.arange(len(ordered_samples)))
    ax.set_xticklabels(ordered_samples, rotation=90)
    ax.set_yticklabels(ordered_samples)

    # Add title
    plt.title('Heatmap of ' + column)

    # Save the plot
    plt.savefig(outfile + ".png", dpi=300, bbox_inches='tight')

    return


def read_pileup(inputfile, contig_length, require_both_strands=True, keep_contigs='All'):
    npos = {"A": 0, "C": 1, "G": 2, "T": 3}

    # generate count matrix
    all_counts = {}
    for contig in contig_length:
        all_counts[contig] = np.zeros((contig_length[contig], 4), dtype=float)

    with open_file(inputfile) as infile:
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
                
                all_counts[line[0]][pos,:] = counts/max(1, np.sum(counts))

    return all_counts

def plot_pairwise_scatter(count_file_A, count_file_B, outfile, require_both_strands=True, min_freq=0.01, keep_contigs='All'):

    # Column names
    columns = ['A', 'C', 'G', 'T']

    # preprocess files to check contig names match and to generate count matrices

    contig_length_A = Counter()
    contig_length_B = Counter()

    # count number of entries in file
    logging.info("Counting entries in pileup files...")
    with open_file(count_file_A) as infile:
        for i, line in enumerate(infile):
            line = line.strip().split()
            if ("All" in keep_contigs) or (line[0]  in keep_contigs):
                if contig_length_A[line[0]] < int(line[1]):
                    contig_length_A[line[0]] = int(line[1])

    with open_file(count_file_B) as infile:
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
        variablesites[contig] = ((fA[contig] + fB[contig])>min_freq).sum(axis=1) > 1

    matches = {}
    for contig in fA:
        temp = pd.DataFrame((fA[contig]>0) & (fB[contig]>0), columns=columns)
        matches[contig] = np.array(temp.melt(ignore_index=False, var_name='allele', value_name='match')['match'])


    # generate pandas dataframe
    logging.info("Generating pandas dataframe...")

    # Iterate over the contigs
    # pd.DataFrame(columns = ['position', 'allmismatch', 'allele', 'frequency', 'sample', 'contig'])
    pdf = None
    for sample, f in zip([count_file_A, count_file_B], [fA, fB]):
        for contig in f:
            # Convert the NumPy matrix to a pandas DataFrame
            df = pd.DataFrame(f[contig], columns=columns)
            df['allmismatch'] = allmismatches[contig]
            df['variable'] = variablesites[contig]

            # Melt the DataFrame into a long format
            long_df = df.melt(ignore_index=False, var_name='allele', value_name='frequency', id_vars=['allmismatch', 'variable'])

            long_df['match'] = matches[contig]
            
            long_df = long_df[long_df['frequency'] >= min_freq]
            long_df = long_df[(long_df['frequency'] <= 1-min_freq) | long_df['variable'] | long_df['allmismatch']] #   
            long_df.reset_index(inplace=True)

            if long_df.shape[0] < 1:
                continue

            long_df.rename(columns={'index': 'position'}, inplace=True)
            long_df['position'] += 1 # convert to 1-based coordinates
            long_df['sample'] = os.path.basename(sample).replace('.txt.gz', '')
            long_df['contig'] = contig
            
            if pdf is None:
                pdf = long_df
            else:
                pdf = pd.concat([pdf, long_df], ignore_index=True)


    # generate plotly figure
    fig = px.scatter(pdf, x="position", y="frequency", facet_col="contig", facet_row="sample", color="allele", symbol="match", opacity=0.7,
                     symbol_sequence= ['circle', 'circle-open'],
                     color_discrete_map={
                        'A': '#e41a1c',
                        'C': '#377eb8',
                        'G': '#4daf4a',
                        'T': '#984ea3'},
                     hover_data=['allele', 'frequency', 'match'])
    fig.update_layout(yaxis_range=[-0.05,1.05])
    
    if len(keep_contigs)>1:
        fig.update_xaxes(matches=None)

    # Save data used in the plot
    pdf.to_csv(outfile + ".csv", index=False)

    offline.plot(fig, filename=outfile + ".html", auto_open=False)

    return

def plot_pairwise_line(count_file_A, count_file_B, outfile, keep_contigs='All', require_both_strands=True, min_freq=0.01, 
                       alpha=0.1, height=7, width=10):

    # Column names
    columns = ['A', 'C', 'G', 'T']

    # preprocess files to check contig names match and to generate count matrices

    contig_length_A = Counter()
    contig_length_B = Counter()

    # count number of entries in file
    logging.info("Counting entries in pileup files...")
    with open_file(count_file_A) as infile:
        for i, line in enumerate(infile):
            line = line.strip().split()
            if ("All" in keep_contigs) or (line[0]  in keep_contigs):
                if contig_length_A[line[0]] < int(line[1]):
                    contig_length_A[line[0]] = int(line[1])

    with open_file(count_file_B) as infile:
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


    consensus_diff = {}
    for contig in fA:
        consensus_diff[contig] = np.argmax(fA[contig], axis=1) != np.argmax(fB[contig], axis=1)

    # generate pandas dataframe
    logging.info("Generating pandas dataframe...")

    # Iterate over the contigs
    # pd.DataFrame(columns = ['position', 'allmismatch', 'allele', 'frequency', 'sample', 'contig'])
    pdf = None
    for sample, f in zip([count_file_A, count_file_B], [fA, fB]):
        for contig in f:
            # Convert the NumPy matrix to a pandas DataFrame
            df = pd.DataFrame(f[contig][consensus_diff[contig], :], columns=columns)

            # Melt the DataFrame into a long format
            long_df = df.melt(ignore_index=False, var_name='allele', value_name='frequency')            
            long_df = long_df[long_df['frequency'] >= min_freq]
            long_df = long_df[(long_df['frequency'] <= 1-min_freq)]
            long_df.reset_index(inplace=True)

            if long_df.shape[0] < 1:
                continue

            long_df.rename(columns={'index': 'position'}, inplace=True)
            long_df['position'] += 1 # convert to 1-based coordinates
            long_df['sample'] = os.path.basename(sample).replace('.txt.gz', '')
            long_df['contig'] = contig
            
            if pdf is None:
                pdf = long_df
            else:
                pdf = pd.concat([pdf, long_df], ignore_index=True)

    logging.info("Generating figure...")

    # Group by 'allele', 'contig', and 'position' to plot each combination as a separate line
    pdf['sample_code'] = (pdf['sample']==os.path.basename(count_file_A).replace('.txt.gz', ''))*1
    pdf = pdf.groupby(['allele', 'contig', 'position']).filter(lambda x: len(x) > 1)
    
    # filter our singletons
    groups = pdf.groupby(['allele', 'contig', 'position'])

    # Prepare data for LineCollection
    # Prepare the LineCollection data
    lines = [np.column_stack([group['sample_code'], group['frequency']]) for name, group in groups]

    # Create a LineCollection from the segments
    lc = LineCollection(lines, linewidths=0.5, alpha=alpha)

    fig, ax = plt.subplots(figsize=(width, height))
    # Add the collection to the plot
    ax.add_collection(lc)

    # Set x-ticks and their labels based on the 'samples' vector
    ax.set_xticks([0,1])  # Set integer ticks
    snames = [os.path.basename(n).split('.')[0] for n in [count_file_A, count_file_B]]
    ax.set_xticklabels(snames, rotation=90)  # Set labels based on 'samples' vector

    # Add labels and title
    ax.set_xlabel('Sample')
    ax.set_ylabel('Frequency')
    ax.set_title('Minor allele frequency by sample')

    # Save the plot
    plt.savefig(outfile + ".png", dpi=300, bbox_inches='tight')

    # Save data used in the plot
    pdf.to_csv(outfile + ".csv", index=False)

    return

def plots(args):

    logging.basicConfig(format='%(asctime)s - %(message)s', 
                        datefmt='%d-%b-%y %H:%M:%S',
                        level=logging.INFO)

    if args.plot_type == "scatter":
        plot_pairwise_scatter(args.input_files[0], 
                                args.input_files[1], 
                                args.output_file,
                                require_both_strands=args.require_both_strands,
                                min_freq=args.min_freq,
                                keep_contigs=args.contigs)
    elif args.plot_type == "line":
        plot_pairwise_line(args.input_files[0], 
                                args.input_files[1], 
                                args.output_file,
                                require_both_strands=args.require_both_strands,
                                min_freq=args.min_freq,
                                keep_contigs=args.contigs,
                                alpha=args.alpha,
                                height=args.height,
                                width=args.width)
        
    elif args.plot_type == "heatmap":
        plot_heatmap(args.input_files[0], 
                    args.output_file,
                    column=args.column_name,
                    threshold=args.threshold,
                    height=args.height,
                    width=args.width)

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
