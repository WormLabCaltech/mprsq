"""
A script to process RNA-seq fastas into pseudoaligned reads using Kallisto.

authors: David Angeles-Albores and Paul W. Sternberg
contact: dangeles at caltech edu
"""
# import os
# import numpy as np
import argparse
# import pandas as pd
# import shutil

# argparse
parser = argparse.ArgumentParser()
parser.add_argument("-d", "--directory", type=str,
                    help="Directory where fasta reads are kept.")
parser.add_argument("-nmin", "--number_min", nargs='?', type=int,
                    help="Minimum number of fasta files to use.")
parser.add_argument("-nmax", "--number_max", nargs='?', type=int,
                    help="Maximum number of fasta files to use.")
parser.add_argument("-b", "--by", nargs='?', type=int,
                    help="n_min to n_max by n.")
parser.add_argument("-i", "--iterations", nargs='?', type=int,
                    help="Iterations per downsampling number.")
parser.add_argument('-ok', '--output_kallisto', type=str,
                    help='kallisto output directory')
parser.add_argument('-sd', '--sleuth_dir', type=str,
                    help='kallisto output directory')
args = parser.parse_args()

# params
if args.directory:
    directory = args.directory
else:
    directory = 'input/rawseq/'

kallisto_dir = args.output_kallisto
sleuth_dir = args.sleuth_dir

if args.number_min:
    nmin = args.number_min
else:
    nmin = -1

if args.number_max:
    nmax = args.number_max
else:
    nmax = -1

if args.by:
    b = args.by
else:
    b = 2

if args.iterations:
    iterations = args.iterations
else:
    iterations = 1

#  python downsampler.py -d input/rawseq/ -nmin 1
# -nmax 6 -b 2 -ok input/kallisto_sampler/ -sd sleuth_sampler/
# --iterations 10

# python src/kallisto_script_generator.py -f kallisto_downsampler_-1.sh
# -d input/kallisto_downsampler/
# -n -1 -o kallisto_sampler/

# python src/r_script_downsampler.py -k input/kallisto_sampler
# -s sleuth_sampler -f sleuth_downsampler.sh
# -g input/library_genotype_mapping.txt

# echo "backend: TkAgg" >> ~/.matplotlib/matplotlibrc

def make_command(n):
    """Write a string with all the correct commands."""
    # run kallisto downsampler.py
    comment = '# kallisto commands\n'
    s1 = 'python src/kallisto_script_generator.py '
    # kallisto fname in __main__
    fname = '-f kallisto_downsampler_{0}.sh '.format(n)
    # directory where fasta reads are
    d = '-d {0} '.format(directory)
    n_ = '-n {0} '.format(n)
    o = '-o {0}\n'.format(kallisto_dir)
    kallisto_py = comment + s1 + fname + d + n_ + o
    # chmod +x kallisto.sh
    kallisto_chmod = 'chmod +x {0}\n'.format(fname[3:])
    # run the kallisto command
    kallisto_sh = 'sh {0}\n'.format(fname[3:])
    kallisto = kallisto_py + kallisto_chmod + kallisto_sh

    # TODO: Sleuth scripts are not in the right directories
    # run the r_script downsampler.py
    comment = '# sleuth prep and analysis\n'
    s1 = 'python src/r_script_downsampler.py '
    k = '-k {0} '.format(kallisto_dir)
    s_dir = '-s {0} '.format(sleuth_dir)
    g = '-g input/library_genotype_mapping.txt '
    fname = '-f sleuth_downsampler.sh\n'
    sleuth_py = comment + s1 + k + s_dir + g + fname
    # chmod +x sleuth.sh
    sleuth_chmod = 'chmod +x {0}\n'.format(sleuth_dir + fname[3:])
    # run the sleuth command
    if sleuth_dir[-1:] == '/':
        f = sleuth_dir + fname[3:]
    else:
        f = sleuth_dir + '/' + fname[3:]
    sleuth_sh = 'sh {0}\n'.format(f)
    sleuth = sleuth_py + sleuth_chmod + sleuth_sh

    # run the python analysis
    analysis_py = 'python src/DownSamplingScript.py '
    n_ = '-n {0} '.format(n)
    k = '-k {0} '.format(kallisto_dir)
    s_dir = '-s {0} '.format(sleuth_dir)
    analysis_py += n_ + k + s_dir + '\n'
    # rm -rf everything
    rmrf = 'rm -rf {0} {1} {2}\n'.format(kallisto_dir, sleuth_dir,
                                         sleuth_dir+fname[3:])

    return kallisto + sleuth + analysis_py + rmrf


def make_commands():
    """Write a bash script with all the right commands."""
    with open('../downsampler.sh', 'w') as f:
        if nmin < 0:
            s = make_command(-1, 0)
            f.write(s)
        else:
            for n in range(nmin, nmax, b):
                for i in range(0, iterations):
                    s = make_command(n)
                    f.write(s)

make_commands()
