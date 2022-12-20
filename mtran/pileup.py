import os
import subprocess
import tempfile


def align_and_pileup(
    reference,
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
):

    # load pileups generated using the bampileup function
    if not quiet:
        print("Generating alignment and pileup...")

    # run aligner
    temp_file = tempfile.NamedTemporaryFile(delete=False, dir=outdir)
    temp_file.close()

    if aligner == "minimap2":
        cmd = "minimap2"
        cmd += " -t " + str(n_cpu)

    if minimap_params is not None:
        cmd += " " + minimap_params
    else:
        cmd += " -ax " + minimap_preset

    cmd += " " + reference
    cmd += " " + r1

    if r2 is not None:
        cmd += " " + r2

    cmd += " | htsbox samview -S -b - | htsbox samsort - > " + temp_file.name

    if not quiet:
        print("running cmd: " + cmd)

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

    if not quiet:
        print("running cmd: " + cmd)

    subprocess.run(cmd, shell=True, check=True)

    # clean up
    os.remove(temp_file.name)

    return
