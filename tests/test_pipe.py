from tracm.pipe import main

import tempfile
import os, sys

def test_pipe(datafolder):

    with tempfile.TemporaryDirectory() as tmpoutdir:

        # rewrite input csv to add temporary paths
        with open(tmpoutdir + "input_data.tsv", 'w') as outfile:
            with open(datafolder + "input_data.tsv", 'r') as infile:
                outfile.write(next(infile))
                for line in infile:
                    line = line.strip().split()
                    line = [line[0]] + [datafolder + p for p in line[1:]]
                    outfile.write(" ".join(line) + '\n')
        
        # run distance script
        sys.argv = [
            "", 
            "-i", tmpoutdir + "input_data.tsv",
            "--database", datafolder + "pneumoexampleDB.zip",
            "-t", "4",
            "-o", tmpoutdir + "test_pipe"
        ]
        main()

        # read output file
        with open(tmpoutdir + 'test_pipe/transmission_clusters.csv', 'r') as infile:
            lines = infile.read()

        # check output
        print(lines)
        assert "sero_1_4_18C_19F_SRR9998163_WGS_of_strep_pneumoniae_Serotype_19F_assembly,0" in lines
        assert "sero_1_19F_SRR9998163_WGS_of_strep_pneumoniae_Serotype_19F_assembly,0" in lines
        assert "sero_1_4_18C_19F_SRR9998186_WGS_of_strep_pneumoniae_Serotype_1_assembly,1" in lines
        assert "sero_1_19F_SRR9998186_WGS_of_strep_pneumoniae_Serotype_1_assembly,1" in lines

    return
