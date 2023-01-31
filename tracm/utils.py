import os
import sys
import argparse
import subprocess
import random
import gzip
import pyfastx as fx


def run_sketch(
    input_files,
    prefix,
    output,
    ksize=51,
    scaled=10000):

    cmd = "sourmash sketch dna"
    cmd += " --merge " + prefix
    cmd += " -p " + f"scaled={scaled},k={ksize},noabund"
    cmd += " -o " + output
    cmd += " " + " ".join(input_files)

    print(f"sketching input files...")
    print(f"command: {cmd}")
    subprocess.run(cmd, shell=True, check=True)

    return



def run_gather(
    input_files,
    databasefile,
    output,
    temp_dir,
    ksize=51,
    scaled=10000,
    threshold_bp=50000,
    max_hits=99999,
    p_match=0.1,
    cache_size=0):

    # Hash query
    run_sketch(
        input_files=input_files,
        prefix='query',
        output=temp_dir + "query.sig",
        ksize=ksize,
        scaled=scaled
    )

    # Run Sourmash Gather
    cmd = "sourmash gather"
    cmd += " -o " + output + ".csv"
    cmd += " --threshold-bp " + str(threshold_bp)
    cmd += " --ignore-abundance"
    cmd += " " + temp_dir + "query.sig"
    cmd += " " + databasefile

    print(f"finding references...")
    print(f"command: {cmd}")
    subprocess.run(cmd, shell=True, check=True)

    # Process results
    references = []
    potential = []
    with open(output + ".csv", "r") as infile:
        # outfile.write("query,reference,f_unique_to_query,f_match_orig\n")
        next(infile)
        for line in infile:
            line = line.strip().split(',')
            line[2] = float(line[2])
            line[0] = float(line[0])
            potential.append(line)
    
    potential = sorted(potential, reverse=True)

    prev=True
    pcov = potential[0][0]
    for line in potential:
        if (line[2] >= p_match) or (prev and (line[0]/pcov >= 0.98)):
            print(f"Using reference: {line[8]}")
            references.append(line[9])
        else:
            prev=False
        pcov = line[0]

    return references


def check_positive_int(value):
    ivalue = int(value)
    if ivalue <= 0:
        raise argparse.ArgumentTypeError("%s is an invalid positive int value" % value)
    return ivalue


def check_positive_float(value):
    ivalue = float(value)
    if ivalue <= 0:
        raise argparse.ArgumentTypeError(
            "%s is an invalid positive float value" % value
        )
    return ivalue

def generate_reads(fasta, outputfile, coverage=10, read_length=300):
    with gzip.open(outputfile, 'wt') as outfile:
        for seq in fx.Fasta(fasta):
            seq_length = len(seq)
            forward = str(seq.seq)
            reverse = str(seq.antisense)
            nreads = max(coverage+10, int((seq_length/read_length) * coverage + 1))
            for i in range(nreads):
                start = random.randint(0, max(0, seq_length - read_length))
                if i%2==0:
                    r = forward[start:(start+read_length)]
                else:
                    r = reverse[start:(start+read_length)]
                outfile.write(f">{seq.name}_read{i}\n{r}\n")
           
    return