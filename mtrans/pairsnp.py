import os
import subprocess


def run_pairsnp(msa, snp_threshold, outputfile, ncpu=1):
    # runs pairsnp and reads result into a numpy array

    cmd = "pairsnp"
    cmd += " -d " + str(snp_threshold)
    cmd += " -t " + str(ncpu)
    cmd += " -s"
    cmd += " " + msa
    cmd += " > " + outputfile

    subprocess.run(cmd, shell=True, check=True)

    # load output
    distances = []
    with open(outputfile, 'r') as infile:
        next(infile)
        for line in infile:
            distances.append([int(t) for t in line.strip().split()])

    return distances