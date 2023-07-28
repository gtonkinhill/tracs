from TRACM import pairsnp

def test_pairsnp(datafolder):
    # test the interface with the c++ version of pairsnp
    distances = pairsnp(fasta=[datafolder + 'ambig.aln'], n_threads=1, dist=10, filter=False)

    assert distances[0] == [0, 0, 0, 0, 1, 1, 1, 2, 2, 3]
    assert distances[1] == [1, 2, 3, 4, 2, 3, 4, 3, 4, 4]
    assert distances[2] == [0, 2, 1, 1, 2, 2, 2, 3, 3, 0]

    return


def test_pairsnp_filt(datafolder):
    # test the interface with the c++ version of pairsnp and recombination filtering
    distances = pairsnp(fasta=[datafolder + 'long_filt.aln'], n_threads=1, dist=1000, filter=True)

    assert distances[0] == [0, 0, 1]
    assert distances[1] == [1, 2, 2]
    assert distances[4] == [2.0, 2.0, 4.0]

    return