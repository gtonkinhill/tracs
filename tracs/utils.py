import os
import sys
import argparse
import logging
import subprocess
import random
import gzip
import pysylph
import pyfastx as fx


def run_sylph_profile(
    input_files,
    databasefile,
    output=None,
    ksize=31,
    c=200,
    min_eff_cov=0.1,       # Used as sequence abundance threshold
    min_ani=0.95       # Added: Sylph excels at ANI-based filtering
):
    """
    Profiles sequencing files against a Sylph database in memory.
    """
    # 1. Load database
    logging.info(f"Loading Sylph database: {databasefile}")
    database = pysylph.Database.load(databasefile) 

    # 2. Sketch input files
    sample = sylph_sketch(input_files, prefix="query", ksize=ksize, c=c)

    # 3. Profile with Sylph
    logging.info("Profiling query against database...")
    profiler = pysylph.Profiler()
    results = profiler.profile(sample, database)

    # Writing a CSV output if specified
    if output:
        with open(f"{output}.csv", "w") as out_csv:
            out_csv.write("reference,ani,sequence_abundance,eff_cov\n")
            for res in results:
                # `res` object attributes map to Sylph's output columns
                out_csv.write(f"{res.genome},{res.ani},{res.seq_abund},{res.eff_cov}\n")

    # 4. Filter results
    references = []
    
    for res in results:
        if res.ani >= min_ani and res.eff_cov >= min_eff_cov:
            logging.info(
                f"Using reference: {res.genome} "
                f"(ANI: {res.ani:.3f}, Abund: {res.seq_abund:.3f}, Eff_Cov: {res.eff_cov:.2f})"
            )
            references.append(res.genome)

    return references


def sylph_sketch(input_files, prefix, ksize=31, c=200):
    """
    Sketches and merges multiple sequencing files into a single sylph sample.
    
    Parameters:
        input_files (list): List of file paths (FASTQ or FASTA)
        prefix (str): The name for the merged sample
        ksize (int): K-mer size (sylph default is 31)
        c (int): Min-spacing compression parameter (sylph default is 200)
    """
    logging.info(f"Sketching {len(input_files)} files into merged sample '{prefix}'...")
    
    # 1. Initialize the sketcher with Sylph-specific parameters
    sketcher = pysylph.Sketcher(k=ksize, c=c)

    # 2. Create a chained generator to stream all files without loading into memory
    def stream_all_reads():
        for file in input_files:
            for record in fx.Fastx(file):
                yield record[1]

    # 3. Sketch with sylph using the streamed reads
    sample = sketcher.sketch_single(name=prefix, reads=stream_all_reads())

    return sample

def run_sketch(input_files, prefix, output, ksize=51, scaled=10000):
    cmd = "sourmash sketch dna"
    cmd += " --merge " + prefix
    cmd += " -p " + f"scaled={scaled},k={ksize},noabund"
    cmd += " -o " + output
    cmd += " " + " ".join(input_files)

    logging.info(f"sketching input files...")
    logging.info(f"command: {cmd}")
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
    cache_size=0,
):
    # Hash query
    run_sketch(
        input_files=input_files,
        prefix="query",
        output=temp_dir + "query.sig",
        ksize=ksize,
        scaled=scaled,
    )

    # Run Sourmash Gather
    cmd = "sourmash gather"
    cmd += " -o " + output + ".csv"
    cmd += " --threshold-bp " + str(threshold_bp)
    cmd += " --ignore-abundance"
    cmd += " " + temp_dir + "query.sig"
    cmd += " " + databasefile

    logging.info(f"finding references...")
    logging.info(f"command: {cmd}")
    subprocess.run(cmd, shell=True, check=True)

    # Process results
    references = []
    potential = []
    with open(output + ".csv", "r") as infile:
        # outfile.write("query,reference,f_unique_to_query,f_match_orig\n")
        next(infile)
        for line in infile:
            line = line.strip().split(",")
            line[2] = float(line[2])
            line[0] = float(line[0])
            potential.append(line)

    potential = sorted(potential, reverse=True)

    prev = True
    pcov = potential[0][0]
    for line in potential:
        if (line[2] >= p_match) or (prev and (line[0] / pcov >= 0.98)):
            logging.debug(line)
            logging.info(f"Using reference: {line[8]}")
            references.append(line[9])
        else:
            prev = False
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
    with gzip.open(outputfile, "wt") as outfile:
        for seq in fx.Fasta(fasta):
            seq_length = len(seq)
            forward = str(seq.seq)
            reverse = str(seq.antisense)
            nreads = max(coverage + 10, int((seq_length / read_length) * coverage + 1))
            for i in range(nreads):
                start = random.randint(0, max(0, seq_length - read_length))
                if i % 2 == 0:
                    r = forward[start : (start + read_length)]
                else:
                    r = reverse[start : (start + read_length)]
                outfile.write(f">{seq.name}_read{i}\n{r}\n")

    return
