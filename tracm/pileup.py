import os
import subprocess
import tempfile
import gzip
import logging
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
    lowdisk=False
):

    # load pileups generated using the bampileup function
    logging.info("Generating alignment and pileup...")

    # Build composite reference
    with open(outdir + 'composite_reference.fasta', 'w') as outfile:
        for ref in references:
            for name, seq in fx.Fasta(references[ref], build_index=False):
                outfile.write('>' + ref + '@' + name + '\n' + seq + '\n')

    # run aligner
    temp_file = tempfile.NamedTemporaryFile(delete=False, dir=outdir)
    temp_file.close()

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
        cmd += " | htsbox samview -S -b - | htsbox samsort -t " + str(n_cpu) + " - > " + temp_file.name
    else:
        cmd += " > " +  outdir + "read_aln.sam"

    logging.info("running cmd: %s", cmd)

    subprocess.run(cmd, shell=True, check=True)

    if not lowdisk:
        cmd = "htsbox samview -S -b " + outdir + "read_aln.sam | htsbox samsort -t " + str(n_cpu) + " - > " + temp_file.name
        logging.info("running cmd: %s", cmd)
        subprocess.run(cmd, shell=True, check=True)

    # run pileup
    cmd = "htsbox pileup -C -s 0"
    cmd += ' -f ' +  outdir + 'composite_reference.fasta'
    cmd += ' -Q ' + str(Q)
    cmd += ' -q ' + str(q)
    cmd += ' -l ' + str(l)
    cmd += ' -S ' + str(S)
    cmd += ' -V ' + str(V)
    cmd += ' -T ' + str(T)
    cmd += ' ' + temp_file.name
    cmd += " > " + outdir + "composite_pileup.txt"

    logging.info("running cmd: %s", cmd)

    subprocess.run(cmd, shell=True, check=True)

    # split into references
    for ref in references:
        with gzip.open(prefix + "_ref_" + str(ref) + "_pileup.txt", 'wt') as outfile:
            with open(outdir + "composite_pileup.txt", 'r') as infile:
                for line in infile:
                    line = line.strip().split('@')
                    if line[0]==ref:
                        outfile.write("@".join(line[1:]) + '\n')

    # clean up
    os.remove(temp_file.name)

    return


def align_and_pileup(
    reference,
    outdir,
    prefix,
    r1,
    r2=None,
    aligner="minimap2",
    minimap_preset="sr",
    minimap_params=None,
    max_div=1,
    Q = 0, #minimum base quality
    q = 0, #minimum mapping quality
    l = 0, #minimum query length
    S = 0, #minimum supplementary alignment length
    V = 1, #ignore queries with per-base divergence >FLOAT [1]
    T = 0, #ignore bases within INT-bp from either end of a read [0]
    n_cpu=1,
    lowdisk=False
):

    # load pileups generated using the bampileup function
    logging.info("Generating alignment and pileup...")

    # run aligner
    temp_file = tempfile.NamedTemporaryFile(delete=False, dir=outdir)
    temp_file.close()

    if aligner == "minimap2":
        cmd = "minimap2"
        cmd += " -t " + str(n_cpu)
        cmd += " -p 1 -N 10"

    if minimap_params is not None:
        cmd += " " + minimap_params
    else:
        cmd += " -ax " + minimap_preset

    cmd += " " + reference
    cmd += " " + r1

    if r2 is not None:
        cmd += " " + r2

    if lowdisk:
        cmd += ' | samtools view -S -b --threads ' + str(n_cpu) + ' --input-fmt-option "filter=[de] < ' + str(max_div) + '" - | samtools sort --threads ' + str(n_cpu) + " - > " + temp_file.name
    else:
        cmd += " > " +  outdir + "read_aln.sam"

    logging.info("running cmd: %s",  cmd)

    subprocess.run(cmd, shell=True, check=True)

    if not lowdisk:
        cmd = 'samtools view -S -b --threads ' + str(n_cpu) + ' --input-fmt-option "filter=[de] < ' + str(max_div) + '" ' + outdir + "read_aln.sam | samtools sort --threads " + str(n_cpu) + " - > " + temp_file.name
        logging.info("running cmd: " + cmd)
        subprocess.run(cmd, shell=True, check=True)

    logging.info("running cmd: %s", cmd)

    subprocess.run(cmd, shell=True, check=True)

    # run pileup
    cmd = "htsbox pileup -C -s 0"
    cmd += ' -f ' + reference
    cmd += ' -Q ' + str(Q)
    cmd += ' -q ' + str(q)
    cmd += ' -l ' + str(l)
    cmd += ' -S ' + str(S)
    cmd += ' -V ' + str(V)
    cmd += ' -T ' + str(T)
    cmd += ' ' + temp_file.name
    cmd += " > " + prefix + "_pileup.txt"

    logging.info("running cmd: %s", cmd)

    subprocess.run(cmd, shell=True, check=True)


    cmd = "gzip " + prefix + "_pileup.txt"
    logging.info("running cmd: %s", cmd)

    subprocess.run(cmd, shell=True, check=True)

    # clean up
    os.remove(temp_file.name)

    return