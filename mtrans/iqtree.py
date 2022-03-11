import os
import subprocess
import pyfastx
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

def run_iqtree_index_cases(msa, clusters, dates, outdir, ncpu=1):
    # runs iqtree on the index case for each transmission cluster

    # identify index case
    cluster_indexes = {}
    for i in dates:
        if (clusters[i] not in cluster_indexes) or (
                cluster_indexes[clusters[i]][1] > dates[i][1]):
            cluster_indexes[clusters[i]] = (i, dates[i][1])
    cluster_indexes = set([t[0] for t in cluster_indexes.values()])

    # create fasta file
    index_fasta = outdir + "index_seqs.fasta"
    date_file = outdir + "index_dates.tab"
    with open(index_fasta, 'w') as outfile, open(date_file, 'w') as datefile:
        for i, seq in enumerate(pyfastx.Fasta(msa, build_index=False)):
            if i in cluster_indexes:
                outfile.write(">cluster_" + str(clusters[i]+1) + "_" + seq[0] +
                              "\n" + seq[1] + "\n")
                datefile.write("cluster_" + str(clusters[i]+1) + "_"  + seq[0] + 
                              "\t" + dates[i][0] + "\n")

    # run iqtree
    cmd = "iqtree2"
    cmd += " -s " + index_fasta
    cmd += " --date " + date_file
    cmd += " -nt " + str(ncpu)
    cmd += " -st DNA"
    cmd += " -m GTR+G"
    cmd += " -quiet"
    cmd += " --date-root 2020-01-01"
    cmd += " -redo"

    subprocess.run(cmd, shell=True, check=True)

    return


def run_iqtree_mrca_cases(msa, clusters, dates, sparse_snp_dist, outdir, ncpu=1):
    # runs iqtree on the index case for each transmission cluster

    # identify unique sequences
    nsamples = len(dates)
    row_ind = []
    col_ind = []
    data = []
    for i, j, d in sparse_snp_dist:
        if d<1:
            row_ind.append(i)
            col_ind.append(i)
            data.append(1)
    sparse_dist_matrix = csr_matrix((data, (row_ind, col_ind)),
                                    shape=(nsamples, nsamples))
    n_components, labels = connected_components(csgraph=sparse_dist_matrix,
                                                directed=False,
                                                return_labels=True)
    index_to_cluster = {}
    for index, cluster in enumerate(labels):
        index_to_cluster[index] = cluster

    print(n_components)

    # identify index case
    cluster_indexes = {}
    for i in dates:
        if (index_to_cluster[i] not in cluster_indexes) or (
                cluster_indexes[index_to_cluster[i]][1] > dates[i][1]):
            cluster_indexes[index_to_cluster[i]] = (i, dates[i][1])
    cluster_indexes = set([t[0] for t in cluster_indexes.values()])

    # create fasta file
    index_fasta = outdir + "unique_index_seqs.fasta"
    date_file = outdir + "unique_index_dates.tab"
    with open(index_fasta, 'w') as outfile, open(date_file, 'w') as datefile:
        datefile.write("name\tdate\n")
        for i, seq in enumerate(pyfastx.Fasta(msa, build_index=False)):
            if i in cluster_indexes:
                outfile.write(">cluster_" + str(index_to_cluster[i]+1) + "_" + seq[0] +
                              "\n" + seq[1] + "\n")
                datefile.write(seq[0] + "\t" + dates[i][0] + "\n")

    # run iqtree
    cmd = "iqtree2"
    cmd += " -s " + index_fasta
    cmd += " --date " + date_file
    cmd += " -nt " + str(ncpu)
    cmd += " -st DNA"
    cmd += " -m GTR+G"
    cmd += " -quiet"
    cmd += " -redo"

    subprocess.run(cmd, shell=True, check=True)

    # collapse to mrca for each cluster
    # tree = Phylo.read(index_fasta + ".treefile", 'newick')


    return