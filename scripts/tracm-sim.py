import os, sys
import argparse
import tempfile
import subprocess
import pyfastx
import random
import numpy as np
import shutil

def generate_genome_pair(genome, distance, outputdir):

    sample_nt = {'A':['C','G','T'],
        'C':['A','G','T'],
        'G':['C','A','T'],
        'T':['C','G','A']}

    # load fasta and calculate genome length
    prefix = os.path.splitext(os.path.basename(genome))[0]
    contigs = []
    l = 0
    sequences = []
    for name, seq in pyfastx.Fasta(genome, build_index=False):
        seq = seq.upper().replace('N','')
        sequences.append((name, seq))
        l += len(seq)

    # generate SNP positions
    snp_locs = np.random.choice(l, size=distance, replace=False)

    l = 0
    with open(outputdir + prefix + '_A.fasta', 'w') as outfile:
        with open(outputdir + prefix + '_B.fasta', 'w') as outfile:
            for name, seq in sequences:
                seqA = list(seq)
                seqB = list(seq)
                temp_locs = snp_locs[(snp_locs>=l) & (snp_locs<(l + len(seq))) ]
                for loc in temp_locs:
                    if random.randint(0, 1)==0:
                        seqA[loc - l] = random.sample(sample_nt[seqA[loc - l]], 1)[0]
                    else:
                        seqB[loc - l] = random.sample(sample_nt[seqB[loc - l]], 1)[0]
                outfileA.write('>' + name + '\n' + ''.join(seqA) + '\n')
                outfileB.write('>' + name + '\n' + ''.join(seqB) + '\n')
                l += len(seq)

    return (outputdir + prefix + '_A.fasta', outputdir + prefix + '_B.fasta')


def generate_reads(genomes, depths, platform, prefix, outputdir, quiet=False):

    # Create temporary directory
    temp_dir = os.path.join(tempfile.mkdtemp(dir=outputdir), "")

    for genome, depth in zip(genomes, depths):
        if depth<=0.01: continue
        temp_prefix = os.path.splitext(os.path.basename(genome))[0]
        if platform=='illumina':
            cmd = 'art_illumina'
            cmd += ' -ss HS25'
            cmd += ' -i ' + genome
            cmd += ' -f ' + str(depth)
            cmd += ' -p -l 150 -m 200 -s 10 --noALN'
            cmd += ' -o ' + temp_dir + temp_prefix
        else:
            cmd = 'badread simulate'
            cmd += ' --reference '+ genome
            cmd += ' --quantity ' + str(depth) + 'x'
            cmd += ' gzip > ' + temp_dir + temp_prefix + '.fastq.gz'

        if not quiet:
            print("running cmd: " + cmd)

        subprocess.run(cmd, shell=True, check=True)


    # concatenate into a single set of read file(s)
    if platform=='illumina':
        for r in ['1', '2']:
            cmd = "cat "
            cmd += temp_dir + '/*'+ r + '.fq'
            cmd += ' | gzip > ' + outputdir + prefix + r + '.fastq.gz'
            if not quiet:
                print("running cmd: " + cmd)
            subprocess.run(cmd, shell=True, check=True)
    else:
        cmd = "zcat "
        cmd += temp_dir + '/*'+ r + '.fastq.gz'
        cmd += ' | gzip > ' + outputdir + prefix + '.fastq.gz'
        if not quiet:
            print("running cmd: " + cmd)
        subprocess.run(cmd, shell=True, check=True)

    # clean up
    shutil.rmtree(temp_dir)

    return

def main():

    parser = argparse.ArgumentParser(description = "Simulates pairs of metagenomic samples with a specified snp distance seperating a given strain")

    io_opts = parser.add_argument_group("Input/output")

    io_opts.add_argument(
        "--genomes",
        dest="fasta_files",
        required=True,
        help="Input fasta files (one for each reference genome)",
        type=os.path.abspath,
        nargs="+",
    )

    io_opts.add_argument(
        "-n",
        dest="mean_genomes",
        required=True,
        help="Mean number of genomes per sample",
        type=int
    )

    io_opts.add_argument(
        "-o",
        "--outdir",
        dest="output_dir",
        required=True,
        help="Output directory",
        type=os.path.abspath,
    )

    io_opts.add_argument(
        "--prefix",
        dest="prefix",
        required=True,
        help="Filename prefix",
        type=str,
    )

    readsim = parser.add_argument_group("Read simulation options")

    readsim.add_argument(
        "--platform",
        dest="platform",
        help="The sequencing platform to simulate reads for. One of 'illumina' or 'nanopore'",
        choices=['illumina', 'nanopore'],
        default='illumina',
        type=str,
    )

    readsim.add_argument(
        "--depth",
        dest="depth",
        help="The cumulative read depth (across all genomes) to simulate (default=500)",
        default=500,
        type=int,
    )

    snpdist = parser.add_argument_group("SNP distance options")

    snpdist.add_argument(
        "--tran-genome",
        dest="tran_genome",
        help="A genome which was involved in the transmission between the two samples",
        type=os.path.abspath,
        required=True,
    )

    snpdist.add_argument(
        '-d',
        "--distance",
        dest="distance",
        help="SNP distance to simulate. (Sites will be chosen randomly)",
        type=int,
        required=True,
    )

    snpdist.add_argument(
        '-b',
        "--background-distance",
        dest="background_distance",
        help="The mean SNP distance of isolates not separated by direct transmission",
        type=int,
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

    args = parser.parse_args()

    # get working directory and create temp directory
    # create directory if it isn't present already
    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)
    # make sure trailing forward slash is present
    args.output_dir = os.path.join(args.output_dir, "")

    # sample genomes ensuring the transmitted genome is in both
    tindex = args.fasta_files.index(args.tran_genome)
    rindex = list(range(len(args.fasta_files)))
    rindex.remove(tindex)

    ngenomesA = max(1, np.random.poisson(args.mean_genomes, size=1)[0])
    ngenomesB = max(1, np.random.poisson(args.mean_genomes, size=1)[0])

    genomesA = np.append(np.random.choice(np.array(rindex), ngenomesA-1), tindex)
    genomesB = np.append(np.random.choice(np.array(rindex), ngenomesB-1), tindex)

    proportionsA = np.random.dirichlet([1] * len(args.fasta_files), size=1)[0]
    proportionsA[np.logical_not(np.in1d(range(len(args.fasta_files)), genomesA))] = 0
    proportionsA = proportionsA/np.sum(proportionsA)

    proportionsB = np.random.dirichlet([1] *len(args.fasta_files), size=1)[0]
    proportionsB[np.logical_not(np.in1d(range(len(args.fasta_files)), genomesB))] = 0
    proportionsB = proportionsB/np.sum(proportionsB)

    sim = []
    distances = []
    for genome in args.fasta_files:
        if genome==args.tran_genome:
            distances.append(args.distance)
            sim.append(generate_genome_pair(genome, args.distance, args.output_dir))
        else:
            d = np.random.poisson(args.background_distance, size=1)[0]
            distances.append(d)
            sim.append(generate_genome_pair(genome, d, args.output_dir))

    # generate reads for sample A
    generate_reads([x[0] for x in sim], 
                    depths = proportionsA*args.depth, 
                    platform = args.platform, 
                    prefix = args.prefix + '_A', 
                    outputdir = args.output_dir, 
                    quiet = args.quiet)

    # generate reads for sample B
    generate_reads([x[1] for x in sim], 
                    depths = proportionsB*args.depth, 
                    platform = args.platform, 
                    prefix = args.prefix + '_B', 
                    outputdir = args.output_dir, 
                    quiet = args.quiet)

    # write out proportions and pairwise distances
    with open(args.output_dir + args.prefix + '_dist_props.csv', 'w') as outfile:
        outfile.write('# simulation command: ' + ' '.join(sys.argv) + '\n')
        outfile.write('genome,propA,propB,distance\n')
        for genome, pA, pB, d in zip(args.fasta_files, proportionsA, proportionsB, distances):
            genome = os.path.splitext(os.path.basename(genome))[0]
            outfile.write(','.join([genome, str(pA), str(pB), str(d)])+'\n')


    return



if __name__ == "__main__":
    main()