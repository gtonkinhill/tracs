import os, sys
import argparse
import numpy as np
from tqdm import tqdm
import gzip
import pysam
import pysamstats
import pyfastx
from numba import jit
from .dirichlet_multinomial import find_dirichlet_priors
from .__init__ import __version__

def preprocess_pileup(alnfiles,
                reference,
                outputfile,
                expected_freq_threshold=None,
                keep_all=True,
                require_both_strands=True,
                min_contig_length=1000,
                max_depth=8000,
                min_mapq=60,
                min_baseq=13,
                quiet=False):

    # load pileups as output by pysamstats
    if not quiet:
        print('Loading reference...')

    # load reference contig lengths
    contig_lengths = []
    total_length = 0
    for name, seq in pyfastx.Fasta(reference, build_index=False):
        l = len(seq)
        if l < min_contig_length: continue
        contig_lengths.append((name, l))
        total_length += l

    attr = [(0, 'A_fwd', 'A_rev'), (1, 'C_fwd', 'C_rev'), (2, 'G_fwd', 'G_rev'), (3, 'T_fwd', 'T_rev')]

    if not quiet:
        print('Loading alignments...')

    # 0 (no evidence) /  1 (found) / 9 (unknown)
    try:
        os.remove(outputfile)
    except OSError:
        pass

    with gzip.open(outputfile, 'ab') as outfile:
        for alnfile in alnfiles:
            fs = os.path.basename(alnfile).split('.')
            ext = fs[-1]
            prefix = '.'.join(fs[:-1])
            if ext=='sam':
                samfile = pysam.AlignmentFile(alnfile, "r")
            elif ext=='bam':
                samfile = pysam.AlignmentFile(alnfile, "rb")
            elif ext=='cram':
                samfile = pysam.AlignmentFile(alnfile, "rc")
            else:
                raise ValueError('File extension is not supported!')

            all_counts = np.zeros((total_length, 4), dtype=float)

            current = 0
            for name, l in contig_lengths:
                pile = pysamstats.load_variation_strand(samfile, 
                    chrom=name, 
                    start=1,
                    end=l,
                    truncate=False,
                    fafile=reference,
                    max_depth=max_depth,
                    min_mapq=min_mapq,
                    min_baseq=min_baseq)
                for i, f, r in attr:
                    counts = pile[f] + pile[r]
                    if require_both_strands:
                        counts[pile[f]<=0] = 0
                        counts[pile[r]<=0] = 0
                    all_counts[current+pile['pos'], i] = counts
                current += l

            alphas = find_dirichlet_priors(all_counts)
            a0 = np.sum(alphas)

            if expected_freq_threshold is None:
                expected_freq_threshold = alphas[1]/(np.sum(alphas) + 50)
            
            if not quiet:
                print('Using frequency threshold:', expected_freq_threshold)
                print('Calculating expected frequency estimates...')

            calculate_posteriors(all_counts, alphas, keep_all, expected_freq_threshold)

            # Filter out those below the threshold
            all_counts[all_counts<expected_freq_threshold] = 0

            # save output to file
            if not quiet:
                print("saving to file...")
            outfile.write(bytes(prefix + ',', 'utf-8'))
            np.savetxt(outfile, all_counts.reshape((1, all_counts.size)), delimiter=',', newline='\n', fmt='%0.5f')
            outfile.write(b"\n")

    return

@jit(nopython=True)
def calculate_posteriors(counts, alphas, keep_all, expected_freq_threshold):
    a0 = np.sum(alphas)
    amin = alphas[1]/a0
    for i in range(counts.shape[0]):
        if keep_all:
            k = counts[i,:]>0
        
        if np.sum(counts[i,:]) == 0:
            counts[i,:] = amin
        else:
            denom = np.sum(counts[i,:]) + a0
            for a, c in enumerate(np.sort(np.unique(counts[i,:]))[::-1]):
                a = alphas[a]
                for j in range(4):
                    if counts[i,j]==c:
                        counts[i,j] = (counts[i,j] + a)/denom

        if keep_all:
            counts[i,:][(counts[i,:]<expected_freq_threshold) & k] = expected_freq_threshold

    return



def get_options(args):

    description = 'Preprocess read alignments ready for mtrans'
    parser = argparse.ArgumentParser(description=description,
                                     prog='bampileup')

    io_opts = parser.add_argument_group('Input/output')
    io_opts.add_argument(
        "-i",
        "--input",
        dest="input_files",
        required=True,
        help="input sam/bam/cram files",
        type=str,
        nargs='+')
    io_opts.add_argument(
        "-r",
        "--ref",
        dest="reference",
        required=True,
        help="input reference genome in fasta format",
        type=str)
    io_opts.add_argument("-o",
                         "--out",
                         dest="output_file",
                         required=True,
                         help="location of an output file",
                         type=str)


    #Pileup options
    pileup = parser.add_argument_group('Pileup')
    pileup.add_argument("-c",
                          "--threshold",
                          dest="threshold",
                          default=None,
                          help=("Minimum posterior read frequency threshold." +
                                " The default is set that a variant at a " +
                                "location is discounted if it is not found " +
                                "with a coverage of ~100x"),
                          type=float)
    
    pileup.add_argument("--min-contig-length",
                          dest="min_contig_length",
                          default=1000,
                          help="Minimum contig length (default=1000)",
                          type=int)

    pileup.add_argument("--max-read-depth",
                          dest="max_read_depth",
                          default=8000,
                          help="Maximum read depth considered (default=8000)",
                          type=int)
    
    pileup.add_argument("--min-mapq",
                          dest="min_mapq",
                          default=60,
                          help="Minimum mapping quality (default=60)",
                          type=int)

    pileup.add_argument("--min-baseq",
                          dest="min_baseq",
                          default=13,
                          help="Minimum base quality (default=13)",
                          type=int)
    
    pileup.add_argument("--filter-all",
                          dest="keep_all",
                          help="turns on filtering of variants with support below the posterior frequency threshold",
                          action='store_false',
                          default=True)

    pileup.add_argument("--single-strand",
                          dest="require_both_strands",
                          help="turns off the requirement that a variant is supported by both strands",
                          action='store_false',
                          default=True)

    # Other options
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

    args = parser.parse_args(args)
    return (args)


def main():
    args = get_options(sys.argv[1:])

    # update file extension if needed
    args.output_file = args.output_file.split('.')[0] + '.csv.gz'

    preprocess_pileup(
            alnfiles = args.input_files,
            reference = args.reference,
            outputfile = args.output_file,
            expected_freq_threshold = args.threshold,
            require_both_strands = args.require_both_strands,
            keep_all = args.keep_all,
            min_contig_length = args.min_contig_length,
            max_depth = args.max_read_depth,
            min_mapq = args.min_mapq,
            min_baseq = args.min_baseq,
            quiet = args.quiet)

    return


if __name__ == '__main__':
    main()
