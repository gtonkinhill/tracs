import os
import sys
import argparse

import subprocess


def run_gather(
    input_files,
    databasefile,
    output,
    temp_dir,
    ksize=51,
    scaled=10000,
    threshold_bp=50000,
    max_hits=99999,
    cache_size=0):

    # Hash query
    cmd = "sourmash sketch dna"
    cmd += " --merge query"
    cmd += " -p " + f"scaled={scaled},k={ksize},noabund"
    cmd += " -o " + temp_dir + "query.sig"
    cmd += " " + " ".join(input_files)

    print(f"sketching input files...")
    print(f"command: {cmd}")
    subprocess.run(cmd, shell=True, check=True)

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
    with open(output + ".csv", "r") as infile:
        # outfile.write("query,reference,f_unique_to_query,f_match_orig\n")
        next(infile)
        for line in infile:
            line = line.strip().split(',')
            references.append(line[9])

    return references
