from cgi import print_form
import os
import sys
import argparse
import time
from .__init__ import __version__

# import sourmash stuff
import screed
import sourmash
from sourmash import MinHash
from sourmash import sourmash_args
from sourmash.logging import error
from sourmash.search import GatherDatabases, format_bp
from sourmash.command_compute import add_seq, set_sig_name, save_siglist
from sourmash.command_sketch import _signatures_for_sketch_factory
from sourmash.signature import SourmashSignature


def gather(
    input_files,
    databasefile,
    outdir,
    prefix,
    ksize=51,
    threshold_bp=50000,
    min_match=0.5,
    prefetch=True,
    max_hits=99999,
    cache_size=0,
):

    # create hash from reads
    signatures_factory = _signatures_for_sketch_factory(
        ["k=51,scaled=1000,noabund"], "dna", False
    )
    sigs = signatures_factory()

    total_seq = 0
    for filename in input_files:
        # consume & calculate signatures
        print(f"... reading sequences from {filename}")

        n = None
        for n, record in enumerate(screed.open(filename)):
            if n % 10000 == 0 and n:
                print(f"\r... {filename} {n}")
            add_seq(sigs, record.sequence, False, False)
            # if n % 20000 == 0 and n: break

        if n is not None:
            print(f"... {filename} {n + 1} sequences")
            total_seq += n + 1
        else:
            print(f"no sequences found in '{filename}'?!")

    set_sig_name(sigs, prefix, prefix)
    print(
        f"calculated 1 signature for {total_seq} sequences taken from {len(input_files)} files"
    )

    # save and load hash. This could do with improving
    save_siglist(sigs, outdir + prefix + ".sig")
    # load the query signature & figure out all the things
    query = sourmash_args.load_query_signature(
        outdir + prefix + ".sig", ksize=ksize, select_moltype="DNA", select_md5=None
    )

    # check query is not empty
    if not len(query.minhash):
        error("no query hashes!? exiting.")
        sys.exit(-1)

    # set up the search databases
    if cache_size == 0:
        cache_size = None
    databases = sourmash_args.load_dbs_and_sigs([databasefile], query, False)

    if not len(databases):
        error("Nothing found to search!")
        sys.exit(-1)

    # run prefetch
    if prefetch:  # note: on by default!
        print("Starting prefetch sweep across databases.")
        prefetch_query = query.copy()
        prefetch_query.minhash = prefetch_query.minhash.flatten()
        noident_mh = prefetch_query.minhash.to_mutable()
        # set up prefetch CSV output, write headers, etc.
        prefetch_csvout_fp = None
        prefetch_csvout_w = None

        counters = []
        for db in databases:
            counter = db.counter_gather(prefetch_query, threshold_bp)

            # subtract found hashes as we can.
            for found_sig in counter.siglist:
                noident_mh.remove_many(found_sig.minhash)

            counters.append(counter)
    else:
        counters = databases
        # we can't track unidentified hashes w/o prefetch
        noident_mh = None

    # Now run gather
    found = []
    weighted_missed = 1
    orig_query_mh = query.minhash
    gather_iter = GatherDatabases(
        query,
        counters,
        threshold_bp=threshold_bp,
        ignore_abundance=True,
        noident_mh=noident_mh,
    )

    for result, weighted_missed in gather_iter:
        if not len(found):  # first result? print header.
            print("")
            print("overlap\t\tp_query\tp_match")
            print("-------\t\t-------\t-------")

        # print interim result & save in `found` list for later use
        pct_query = "{:.1f}%".format(result.f_unique_weighted * 100)
        pct_genome = "{:.1f}%".format(result.f_match * 100)
        name = result.match._display_name(40)

        print(f"{format_bp(result.intersect_bp)}\t\t{pct_query}\t{pct_genome}\t{name}")
        found.append(result)

        if len(found) >= max_hits:
            break

    # report on thresholding -
    if gather_iter.query:
        # if still a query, then we failed the threshold.
        print(f"found less than {format_bp(threshold_bp)} in common. => exiting")

    # basic reporting:
    print(f"\nfound {len(found)} matches total;")
    if len(found) == max_hits:
        print(f"(truncated gather because the maximum hit limit was reached!)")

    p_covered = (1 - weighted_missed) * 100
    print(f"the recovered matches hit {p_covered:.1f}% of the query (unweighted)\n")
    if gather_iter.scaled != query.minhash.scaled:
        print(
            f"WARNING: final scaled was {gather_iter.scaled}, vs query scaled of {query.minhash.scaled}"
        )

    with open(outdir + "sourmash_hits.csv", "w") as outfile:
        outfile.write("query,reference,f_unique_to_query,f_match_orig\n")
        for result in found:
            outfile.write(
                ",".join(
                    [
                        result.query_name,
                        result.name,
                        str(result.f_unique_to_query),
                        str(result.f_match_orig),
                    ]
                )
                + "\n"
            )

    return


def main():

    parser = argparse.ArgumentParser(
        description="Uses sourmash to identify reference genomes within a read set",
        prog="gather",
    )

    parser.add_argument(
        "-i",
        "--input",
        dest="input_files",
        required=True,
        help="path to query signature",
        type=os.path.abspath,
        nargs="+",
    )

    parser.add_argument(
        "--database",
        dest="database",
        help="path to database signatures",
        type=os.path.abspath,
        default=None,
    )

    parser.add_argument(
        "-o",
        "--out_dir",
        dest="output_dir",
        required=True,
        help="location of an output directory",
        type=os.path.abspath,
    )

    parser.add_argument(
        "-p",
        "--prefix",
        dest="prefix",
        default=None,
        help="prefix to describe the input sample read files",
        type=str,
    )

    parser.add_argument(
        "-t",
        "--threads",
        dest="n_cpu",
        help="number of threads to use (default=1)",
        type=int,
        default=1,
    )

    parser.add_argument(
        "--quiet",
        dest="quiet",
        help="turns off some console output",
        action="store_true",
        default=False,
    )

    parser.add_argument(
        "--version", action="version", version="%(prog)s " + __version__
    )

    args = parser.parse_args()

    # create directory if it isn't present already
    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)
    # make sure trailing forward slash is present
    args.output_dir = os.path.join(args.output_dir, "")

    if args.prefix is None:
        args.prefix = os.path.basename(args.input_files[0]).split(".")[0]

    gather(
        input_files=args.input_files,
        databasefile=args.database,
        outdir=args.output_dir,
        prefix=args.prefix,
    )

    return


if __name__ == "__main__":
    main()
