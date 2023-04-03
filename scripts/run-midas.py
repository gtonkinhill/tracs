import os
import sys
import argparse
import subprocess
import glob
import tempfile
import shutil
from collections import defaultdict, Counter

import os
import subprocess
import tempfile
import lz4.frame



def main():
    
    parser = argparse.ArgumentParser()
    parser.description = "Runs Midas2 on a pair of samples to infer a SNP distance"

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
        help="path to midas2 database",
        # required=True,
        default="/data1/gerryt/tracm-data/simulations/my_midasdb_gtdb",
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

    dirA = args.output_dir + "sampleA/"
    if not os.path.exists(dirA):
        os.mkdir(dirA)

    dirB = args.output_dir + "sampleB/" 
    if not os.path.exists(dirB):
        os.mkdir(dirB)

    # identify species
    for n, dirT, ifiles in [('A', dirA, args.input_files_A), ('B', dirB, args.input_files_B)]:
        cmd = "midas2 run_species"
        cmd += " --sample_name sample" + n
        cmd += " -1 " + ifiles[0]
        cmd += " -2 " + ifiles[1]
        cmd += " --midasdb_name gtdb"
        cmd += " --midasdb_dir " + args.refDB
        cmd += " --num_cores " + str(args.n_cpu)
        cmd += " " + args.output_dir

        print("running cmd: " + cmd)
        subprocess.run(cmd, shell=True, check=True)
   

    # align with midas2
    for n, dirT, ifiles in [('A', dirA, args.input_files_A), ('B', dirB, args.input_files_B)]:
        cmd = "midas2 run_snps"
        cmd += " --sample_name sample" + n
        cmd += " -1 " + ifiles[0]
        cmd += " -2 " + ifiles[1]
        cmd += " --midasdb_name gtdb"
        cmd += " --midasdb_dir " + args.refDB
        cmd += " --num_cores " + str(args.n_cpu)
        cmd += " " + args.output_dir

        print("running cmd: " + cmd)
        subprocess.run(cmd, shell=True, check=True)

    # create sample manifest file
    with open(args.output_dir + "midas_list_of_samples.tsv", 'w') as outfile:
        outfile.write("sample_name\tmidas_outdir\n")
        outfile.write("sampleA\t" + args.output_dir + "\n")
        outfile.write("sampleB\t" + args.output_dir + "\n")

    # run midas pop SNV calls
    cmd = "midas2 merge_snps "
    cmd += "--samples_list " + args.output_dir + "midas_list_of_samples.tsv"
    cmd += " --midasdb_name gtdb"
    cmd += " --midasdb_dir " + args.refDB
    cmd += " --site_ratio 100"
    cmd += " --genome_coverage 0.1"
    cmd += " --genome_depth 0.01"
    cmd += " --snv_type rare"
    cmd += " --site_prev 1"
    cmd += " --snp_pooled_method abundance"
    cmd += " --num_cores " + str(args.n_cpu)
    cmd += " " + args.output_dir + "merge"

    print("running cmd: " + cmd)
    subprocess.run(cmd, shell=True, check=True)

    # Load species IDs
    spIDs = {}
    with open(args.refDB + "/metadata.tsv", 'r') as infile:
        for line in infile:
            line = line.split('\t')
            spIDs[line[1]] = line[4]

    # Load coverage
    cov = {}
    with open(args.output_dir + "merge/snps/snps_summary.tsv", 'r') as infile:
        for line in infile:
            line = line.strip().split()
            cov[(line[0], line[1])] = line[-2:]

    # count snps and write sumamry report
    with open(args.output_dir + "summary_snp_dist.tsv", 'w') as outfile:
        outfile.write("speciesID\tsnp_dist\tspecies_name\tfcovA\tmcovA\tfcovB\tmcovB\n")
        snp_files = glob.glob(args.output_dir + "merge/snps/*/*snps_info.tsv.lz4")
        for sf in snp_files:
            species = os.path.basename(sf).replace(".snps_info.tsv.lz4", "")
            with lz4.frame.open(sf, mode='r') as fp:
                next(fp)
                nsnp = 0
                for line in fp:
                    line = line.strip().split()
                    is_snp=1
                    tot = 0
                    for i in range(8,12):
                        c = int(line[i])
                        tot += c
                        if c >= 2: 
                            is_snp=0
                    if tot>1:
                        nsnp+=is_snp
            outfile.write("\t".join([species, str(nsnp), spIDs[species]] + 
                cov[("sampleA", species)] + 
                cov[("sampleB", species)]) + '\n')

    return


if __name__ == '__main__':
    main()
