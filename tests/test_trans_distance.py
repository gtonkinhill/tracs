from tracm.distance import main

import tempfile
import os, sys
import numpy as np

def test_trans_distance(datafolder):

    with tempfile.TemporaryDirectory() as tmpoutdir:
        
        # run distance script
        sys.argv = [
            "", "--msa", datafolder + "ambig.aln",
            "--meta", datafolder + 'dates_ambig.csv', 
            "-o", tmpoutdir + "distances.csv", "-K", "10", "--snp_threshold", "5"
        ]
        main()

        # read output file
        with open(tmpoutdir + 'distances.csv', 'r') as infile:
            lines = infile.readlines()

        # check output
        print(lines)
        line1 = lines[1].strip().split(',')
        line2 = lines[2].strip().split(',')

        # check date difference
        assert abs(float(line1[2])-0.002737907006988508) < 1e-6
        assert abs(float(line2[2])-0.002737907006988508) < 1e-6

        # check snp distance
        assert int(line1[3]) == 0
        assert int(line2[3]) == 2
        
        # check direct transmission distance
        assert abs(float(line1[4])-0.23794988406662973) < 1e-6
        assert abs(float(line2[4])-0.024467137572328577) < 1e-6

        # check expected transmission distance (K)
        assert abs(float(line1[5])-2.6335200453700187) < 1e-6
        assert abs(float(line2[5])-7.315670110063259) < 1e-6

    return

# TODO: Manually check these
# sampleA,sampleB,date difference,SNP distance,transmission distance,expected K
# seq1,seq2,0.002737907006988508,0,0.23794988406662973,2.6335200453700187
# seq1,seq3,0.002737907006988508,2,0.024467137572328577,7.315670110063259