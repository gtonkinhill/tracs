import os
import subprocess
import tempfile
import argparse
from .__init__ import __version__

def align_and_pileup(reference,
                outdir,
                prefix,
                r1, 
                r2=None,
                aligner='minimap2',
                minimap_preset='sr',
                minimap_params=None,
                n_cpu=1,
                quiet=False):

    # load pileups generated using the bampileup function
    if not quiet:
        print('Generating alignment and pileup...')

    # run aligner
    temp_file = tempfile.NamedTemporaryFile(delete=False, dir=outdir)
    temp_file.close()

    if aligner=='minimap2':
        cmd = "minimap2" 
        cmd += " -t " + str(n_cpu)

    if minimap_params is not None:
        cmd += ' ' + minimap_params
    else:
        cmd += ' -ax ' + minimap_preset
    
    cmd += " " + reference
    cmd += " " + r1
    
    if r2 is not None:
        cmd += ' ' + r2
    
    cmd += ' | htsbox samview -S -b - | htsbox samsort - > ' + temp_file.name

    if not quiet:
        print("running cmd: " + cmd)

    subprocess.run(cmd, shell=True, check=True)


    # run pileup
    cmd = "htsbox pileup -C "
    cmd += temp_file.name
    cmd += ' > ' + prefix + '_pileup.txt'

    if not quiet:
        print("running cmd: " + cmd)

    subprocess.run(cmd, shell=True, check=True)

    # clean up
    os.remove(temp_file.name)

    return 


def main():

    parser = argparse.ArgumentParser(description='Generates pileup using minimap2 & htsbox',
                                     prog='pileup')

    parser.add_argument("--r1",
                        dest="r1",
                        required=True,
                        help="fasta/q file 1",
                        type=os.path.abspath)
    
    parser.add_argument("--r2",
                        dest="r2",
                        help="fasta/q file 2 (optional)",
                        type=os.path.abspath,
                        default=None)
    
    parser.add_argument("--reference",
                        dest="reference",
                        required=True,
                        help="reference file in fasta format",
                        type=os.path.abspath)

    parser.add_argument("-o",
                         "--out_dir",
                         dest="output_dir",
                         required=True,
                         help="location of an output directory",
                         type=os.path.abspath)
    
    parser.add_argument("--prefix",
                         dest="prefix",
                         help="prefix to use",
                         default='test',
                         type=str)
    
    parser.add_argument("--minimap_preset",
                         dest="minimap_preset",
                         help="minimap preset to use - one of 'sr' (default), 'map-ont' or 'map-pb'",
                         default='sr',
                         type=str)

    parser.add_argument("-t",
                        "--threads",
                        dest="n_cpu",
                        help="number of threads to use (default=1)",
                        type=int,
                        default=1)

    parser.add_argument("--quiet",
                        dest="quiet",
                        help="turns off some console output",
                        action='store_true',
                        default=False)

    parser.add_argument('--version',
                        action='version',
                        version='%(prog)s ' + __version__)

    args = parser.parse_args()

    # create directory if it isn't present already
    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)
    # make sure trailing forward slash is present
    args.output_dir = os.path.join(args.output_dir, "")
    

    align_and_pileup(
                reference=args.reference,
                outdir=args.output_dir,
                prefix=args.prefix,
                r1=args.r1, 
                r2=args.r2,
                aligner='minimap2',
                minimap_preset=args.minimap_preset,
                minimap_params=None,
                n_cpu=args.n_cpu,
                quiet=args.quiet)

    return


if __name__ == '__main__':
    main()
