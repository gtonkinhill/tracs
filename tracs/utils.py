import os
import sys
import argparse
import logging
import subprocess
import random
import gzip
import csv
import pyfastx as fx
import sourmash
from sourmash.exceptions import SourmashError


# extensions stripped from sylph's Genome_file so names match build_db.genome_prefix
GENOME_EXTENSIONS = (".fasta.gz", ".fna.gz", ".fa.gz", ".fasta", ".fna", ".fa")


def strip_genome_extension(path):
    name = os.path.basename(path)
    for ext in GENOME_EXTENSIONS:
        if name.endswith(ext):
            return name[: -len(ext)]
    return name


def run_sylph_profile(
    input_files,
    databasefile,
    output_dir,
    prefix="query",
    min_eff_cov=0.1,
    min_ani=95.0,  # sylph reports ANI as a percentage (e.g. 99.20)
    n_cpu=1,
):
    """Profile sequencing files against a sylph database.

    Returns the reference names (genome file name without extension) that pass
    the ANI and effective coverage thresholds; an empty list if none do.
    """
    if not input_files:
        raise ValueError("At least one input file is required.")

    out_tsv = os.path.join(output_dir, f"{prefix}_sylph_profile.tsv")

    cmd = ["sylph", "profile", str(databasefile), "-t", str(n_cpu)]

    if len(input_files) == 2:
        cmd.extend(["-1", str(input_files[0]), "-2", str(input_files[1])])
    else:
        cmd.extend(map(str, input_files))

    cmd.extend(["-o", out_tsv])

    logging.info("Profiling query against database with sylph...")
    logging.info("Command: %s", subprocess.list2cmdline(cmd))
    subprocess.run(cmd, check=True)

    references = []
    with open(out_tsv, newline="") as infile:
        for row in csv.DictReader(infile, delimiter="\t"):
            genome = row.get("Genome_file")
            try:
                ani = float(row.get("Adjusted_ANI") or 0)
                eff_cov = float(row.get("Eff_cov") or 0)
                seq_abund = float(row.get("Sequence_abundance") or 0)
            except (TypeError, ValueError):
                continue

            if not genome:
                continue

            name = strip_genome_extension(genome)
            if ani >= min_ani and eff_cov >= min_eff_cov:
                logging.info(
                    "Using reference: %s (ANI: %.2f, Abund: %.3f, Eff_Cov: %.2f)",
                    name, ani, seq_abund, eff_cov,
                )
                references.append(name)
            else:
                logging.info(
                    "Skipping sylph hit below thresholds: %s (ANI: %.2f, Eff_Cov: %.2f)",
                    name, ani, eff_cov,
                )

    return list(dict.fromkeys(references))


def is_valid_sourmash_db(filepath):
    try:
        db = sourmash.load_file_as_index(filepath)
        manifest = db.manifest
        if manifest is None:
            return False
        n = len(db)  # forces manifest/signature enumeration
        if n == 0:
            return False
        return True
    except Exception:
        # Catch broadly: zipfile.BadZipFile, json.JSONDecodeError,
        # KeyError, ValueError, SourmashError, FileNotFoundError, etc.
        # can all surface once parsing is forced.
        return False

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

    # sourmash does not write an output file if it finds no matches;
    # the caller warns and skips the sample
    if not os.path.isfile(output + ".csv"):
        logging.info(
            f"sourmash found no reference genomes sharing >= {threshold_bp}bp with the query."
        )
        return []

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

    if not potential:
        return []

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
