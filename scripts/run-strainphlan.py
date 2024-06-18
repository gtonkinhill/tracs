import os
import sys
import argparse
import subprocess
import glob
import tempfile
import shutil
from collections import defaultdict, Counter
import pyfastx as fx




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
        "--refDB",
        dest="refDB",
        help="path to Metaphlan database",
        # required=True,
        default="/data1/gerryt/tracm-data/simulations/metaphlan/metaphlanDB",
        type=os.path.abspath
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

    # Create temporary directory
    temp_dir = os.path.join(tempfile.mkdtemp(dir=args.output_dir), "")

    dirA = args.output_dir + "sampleA/"
    if not os.path.exists(dirA):
        os.mkdir(dirA)

    dirB = args.output_dir + "sampleB/" 
    if not os.path.exists(dirB):
        os.mkdir(dirB)

    # align with metaphlan
    for n, dirT in [('A', dirA), ('B', dirB)]:
        if os.path.exists(dirT + n + "_metagenome.bowtie2.bz2"):
            cmd = "metaphlan --input_type bowtie2out"
            cmd += " " + dirT + n + "_metagenome.bowtie2.bz2"
        else:
            # metaphlan metagenome.fastq --bowtie2out metagenome.bowtie2.bz2 --nproc 5 --input_type fastq -o profiled_metagenome.txt
            cmd = "metaphlan --input_type fastq"
            if n=='A':
                cmd += " " + ','.join(args.input_files_A)
            else:
                cmd += " " + ','.join(args.input_files_B)
            cmd += " --bowtie2out " + dirT + n + "_metagenome.bowtie2.bz2"
        cmd += " --bowtie2db " + args.refDB
        cmd += " -s " + dirT + n + "_metagenome.sam.bz2"
        cmd += " -o " + dirT + "profiled_metagenome.txt"
        cmd += " --nproc " + str(args.n_cpu)
        print("running cmd: " + cmd)
        subprocess.run(cmd, shell=True, check=True)

    # run consensus calling script
    dirCM = args.output_dir + "consensus_markers/"
    if not os.path.exists(dirCM):
        os.mkdir(dirCM)
    cmd = "sample2markers.py" 
    cmd += " -d " + args.refDB
    cmd += " -i " + args.output_dir + "sample*/*metagenome.sam.bz2"
    cmd += " -o " + dirCM
    cmd += " --tmp " + temp_dir
    cmd += " -n " + str(args.n_cpu)
    print("running cmd: " + cmd)
    subprocess.run(cmd, shell=True, check=True)

    # add additional marker files to trick strainphlan into running the MSA
    shutil.copyfile(args.output_dir + "consensus_markers/A_metagenome.pkl", 
                    args.output_dir + "consensus_markers/temp1_A_metagenome.pkl")
    shutil.copyfile(args.output_dir + "consensus_markers/A_metagenome.pkl", 
                    args.output_dir + "consensus_markers/temp2_A_metagenome.pkl")

    
    species_name = {}
    for dirT in [dirA, dirB]:
        with open(dirT + "profiled_metagenome.txt", 'r') as inputfile:
            for line in inputfile:
                if line[0]=="#": continue
                spec = line.split()[0].split("|")[-1]
                if "t__" in spec:
                    species_name[spec] = line.split()[0].split("|")[-2]

    # identify species to pull out
    emDB = glob.glob(args.refDB + "/*.pkl")[0]
    dirPA = args.output_dir + "strainphlan_pa/"
    if not os.path.exists(dirPA):
        os.mkdir(dirPA)

    cmd = "strainphlan"
    cmd += " -s " + args.output_dir + "consensus_markers/*.pkl"
    cmd += " -o " + dirPA
    cmd += " --print_clades_only"
    cmd += " -d " + emDB
    # cmd += " --marker_in_n_samples 100 --sample_with_n_markers 1"
    print("running cmd: " + cmd)
    subprocess.run(cmd, shell=True, check=True)

    specs = []
    with open(dirPA + "print_clades_only.tsv", 'r') as infile:
        next(infile)
        for line in infile:
            if "t__SGB" in line:
                specs.append(line.split()[0])

    print("specs:", specs)

    # extract markers
    dirMarkers = args.output_dir + "db_markers/"
    if not os.path.exists(dirMarkers):
        os.mkdir(dirMarkers)
    
    for spec in specs:
        cmd = "extract_markers.py"
        cmd += " -c " + spec
        cmd += " -o " + dirMarkers
        cmd += " -d " + emDB
        print("running cmd: " + cmd)
        subprocess.run(cmd, shell=True, check=True)

    # run strainphlan
    dirMA = args.output_dir + "strainphlan_output/"
    if not os.path.exists(dirMA):
        os.mkdir(dirMA)
    for spec in specs:
        if not os.path.exists(dirMA + spec):
            os.mkdir(dirMA + spec)
        cmd = "strainphlan"
        cmd += " -s " + args.output_dir + "consensus_markers/*.pkl"
        cmd += " -m " + args.output_dir + "db_markers/" + spec + ".fna"
        cmd += " -o " + dirMA + spec
        cmd += " -c " + spec
        cmd += " -d " + args.refDB
        # cmd += " --marker_in_n_samples 100 --sample_with_n_markers 1"
        print("running cmd: " + cmd)
        subprocess.run(cmd, shell=True, check=True)

    # run pairsnp and summarise
    with open(args.output_dir + "strainphlan_distances.tsv", 'w') as outfile:
        outfile.write("species,sample1,sample2,snp_dist\n")
        for spec in specs:
            cmd = "pairsnp -c -s "
            cmd += dirMA + spec + "/" + spec + ".StrainPhlAn4_concatenated.aln"
            cmd += " > " + dirMA + spec + "/" + spec + "_dist.csv"
            print("running cmd: " + cmd)
            subprocess.run(cmd, shell=True, check=True)

            with open(dirMA + spec + "/" + spec + "_dist.csv", 'r') as infile:
                for line in infile:
                    if "temp" in line: continue
                    line = line.replace("A_metagenome", "A_" + spec)
                    line = line.replace("B_metagenome", "B_" + spec)
                    outfile.write(species_name[spec] + "," + line)

    shutil.rmtree(temp_dir)

    return


if __name__ == '__main__':
    main()
