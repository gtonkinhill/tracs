import os
import sys
import argparse
import subprocess
import tempfile
from collections import defaultdict, Counter

import os
import subprocess
import tempfile
import pyfastx as fx


def align_and_pileup_composite(
    references,
    outdir,
    prefix,
    r1,
    r2=None,
    aligner="minimap2",
    minimap_preset="sr",
    minimap_params=None,
    Q = 0, #minimum base quality
    q = 0, #minimum mapping quality
    l = 0, #minimum query length
    S = 0, #minimum supplementary alignment length
    V = 1, #ignore queries with per-base divergence >FLOAT [1]
    T = 0, #ignore bases within INT-bp from either end of a read [0]
    n_cpu=1,
    quiet=False,
    lowdisk=False
):

    # load pileups generated using the bampileup function
    if not quiet:
        print("Generating alignment and pileup...")

    # Build composite reference
    with open(outdir + 'composite_reference.fasta', 'w') as outfile:
        with open(outdir + 'composite_reference.txt', 'w') as binoutfile:
            for i, ref in enumerate(references):
                for name, seq in fx.Fasta(references[ref], build_index=False):
                    outfile.write('>' + ref + '@' + name + '\n' + seq + '\n')
                    binoutfile.write(ref + '@' + name + '\t' + ref +'\n')

    # run aligner
    temp_file = outdir + prefix + "_composite_aln.bam"

    if os.path.exists(temp_file):
        return(temp_file)

    if aligner == "minimap2":
        cmd = "minimap2"
        cmd += " -t " + str(n_cpu)
        cmd += " -p 1 -N 10"
    else:
        raise ValueError("Minimap2 is the only currently supported aligner!")

    if minimap_params is not None:
        cmd += " " + minimap_params
    else:
        cmd += " -ax " + minimap_preset

    cmd += " " + outdir + 'composite_reference.fasta'
    cmd += " " + r1

    if r2 is not None:
        cmd += " " + r2

    if lowdisk:
        cmd += ' | samtools view -S -b --threads ' + str(n_cpu) + " - | samtools sort --threads " + str(n_cpu) + " - > " + temp_file
    else:
        cmd += " > " +  outdir + "read_aln.sam"

    if not quiet:
        print("running cmd: " + cmd)

    subprocess.run(cmd, shell=True, check=True)

    if not lowdisk:
        cmd = 'samtools view -S -b --threads ' + str(n_cpu) + " " + outdir + "read_aln.sam | samtools sort --threads " + str(n_cpu) + " - > " + temp_file
        if not quiet:
            print("running cmd: " + cmd)
        subprocess.run(cmd, shell=True, check=True)

    if not quiet:
        print("running cmd: " + cmd)

    # run pileup
    # cmd = "htsbox pileup -C -s 0"
    # cmd += ' -f ' +  outdir + 'composite_reference.fasta'
    # cmd += ' -Q ' + str(Q)
    # cmd += ' -q ' + str(q)
    # cmd += ' -l ' + str(l)
    # cmd += ' -S ' + str(S)
    # cmd += ' -V ' + str(V)
    # cmd += ' -T ' + str(T)
    # cmd += ' ' + temp_file.name
    # cmd += " > " + outdir + "composite_pileup.txt"

    # if not quiet:
    #     print("running cmd: " + cmd)

    # subprocess.run(cmd, shell=True, check=True)

    # # split into references
    # for ref in references:
    #     with open(prefix + "_ref_" + str(ref) + "_pileup.txt", 'w') as outfile:
    #         with open(outdir + "composite_pileup.txt", 'r') as infile:
    #             for line in infile:
    #                 line = line.strip().split('@')
    #                 if line[0]==ref:
    #                     outfile.write("@".join(line[1:]) + '\n')

    return (temp_file)



def main():
    
    parser = argparse.ArgumentParser()
    parser.description = "Uses sourmash to identify reference genomes within a read set and then aligns reads to each reference using minimap2"

    io_opts = parser.add_argument_group("Input/output")

    io_opts.add_argument(
        "--inputA",
        dest="input_files_A",
        required=True,
        help="path to fastq files in sample A",
        type=os.path.abspath,
        nargs="+",
    )

    io_opts.add_argument(
        "--inputB",
        dest="input_files_B",
        required=True,
        help="path to fastq files in sample B",
        type=os.path.abspath,
        nargs="+",
    )

    io_opts.add_argument(
        "--references",
        dest="references",
        help="path to reference fastas",
        required=True,
        type=os.path.abspath,
        nargs="+",
    )

    io_opts.add_argument(
        "-o",
        "--output",
        dest="output_dir",
        required=True,
        help="location of an output directory",
        type=os.path.abspath,
    )

    io_opts.add_argument(
        "-p",
        "--prefix",
        dest="prefix",
        default=None,
        help="prefix to describe the input sample read files",
        type=str,
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

    args = parser.parse_args()


    # run alignment
    
    # set up references and subdirectories
    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)
    # make sure trailing forward slash is present
    args.output_dir = os.path.join(args.output_dir, "")

    dirA = args.output_dir + "sampleA/"
    if not os.path.exists(dirA):
        os.mkdir(dirA)

    dirB = args.output_dir + "sampleB/" 
    if not os.path.exists(dirB):
        os.mkdir(dirB)

    references = {}
    for ref in args.references:
        print (ref)
        print (os.path.splitext(os.path.basename(ref))[0])
        references[os.path.splitext(os.path.basename(ref))[0]] = ref


    # this uses the same commands as in tracs but with a composite reference

    # sampleA 
    bamA = align_and_pileup_composite(
        references,
        dirA,
        'sampleA',
        r1=args.input_files_A[0],
        r2=args.input_files_A[1],
        n_cpu=args.n_cpu)

    # sampleB
    bamB = align_and_pileup_composite(
        references,
        dirB,
        'sampleB',
        r1=args.input_files_B[0],
        r2=args.input_files_B[1],
        n_cpu=args.n_cpu)

    # profile with inStrain
    for bam, dirT in [(bamA, dirA), (bamB, dirB)]:
        cmd = "inStrain profile"
        cmd += " " + bam
        cmd += " " + dirA + 'composite_reference.fasta'
        cmd += " --stb " + dirA + 'composite_reference.txt'
        cmd += " -o " + dirT
        cmd + " --skip_plot_generation"
        cmd += " -p " + str(args.n_cpu)
        # cmd += " --pairing_filter all_reads"
        print("running cmd: " + cmd)
        subprocess.run(cmd, shell=True, check=True)

    # distance estimation
    cmd = "inStrain compare"
    cmd += " -i" + " " + dirA + " " + dirB
    cmd += " --min_cov 2"
    cmd += " -o " + args.output_dir + "compare_instrain"
    cmd += " -p " + str(args.n_cpu)
    print("running cmd: " + cmd)
    subprocess.run(cmd, shell=True, check=True)

    # summarise by genomes
    # scaffold        name1   name2   coverage_overlap        compared_bases_count    percent_genome_compared length  consensus_SNPs  population_SNPs popANI  conANI
    counts = defaultdict(lambda: [0,0,0,0])
    with open(args.output_dir + "/compare_instrain/output/compare_instrain_comparisonsTable.tsv", 'r') as infile:
        h = next(infile).strip().split()
        for line in infile:
            line = line.strip().split()
            ref = line[0].split('@')[0]
            for i, col in enumerate([4,6,7,8]):
                counts[ref][i] += int(line[col])

    with open(args.output_dir + "compare_instrain_summary.tsv", 'w') as outfile:
        outfile.write("\t".join(h[k] for k in [0,4,6,7,8]) + '\n')
        for ref in counts:
            outfile.write("\t".join([ref] + [str(k) for k in counts[ref]]) + '\n')

    return


if __name__ == '__main__':
    main()
