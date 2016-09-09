"""
A script to process RNA-seq fastas into pseudoaligned reads using Kallisto.

authors: David Angeles-Albores and Paul W. Sternberg
contact: dangeles at caltech edu
"""
import os
import numpy as np
import argparse

# get the read number from terminal
parser = argparse.ArgumentParser()
parser.add_argument("-d", "--directory", type=str,
                    help="Directory where fasta reads are kept.")
parser.add_argument("-l", "--length", nargs='?', type=int,
                    help="Fragment length.")
parser.add_argument("-s", "--sigma", nargs='?', type=int,
                    help="Standard deviation of fragment.")
parser.add_argument("-btstrp", "--bootstrap", nargs='?', type=int,
                    help="Number of bootstraps to perform.")
parser.add_argument("-t", "--threads", nargs='?', type=int,
                    help="Number of threads to use.")
parser.add_argument("-o", "--output_dir", nargs='?', type=str,
                    help="Directory to place output in.")
args = parser.parse_args()

directory = args.directory

# params
directory = args.directory

if args.length:
    length = args.length
else:
    length = 180

if args.sigma:
    sigma = args.sigma
else:
    sigma = 60

if args.bootstrap:
    btstrp = args.bootstrap
else:
    btstrp = 200

if args.threads:
    thrds = args.threads
else:
    thrds = 4

if args.output_dir:
    output_dir = args.output_dir
else:
    output_dir = '../input/kallisto_all/'

# sequences:
seqs = next(os.walk(directory))[1]


def explicit_kallisto(directory, files, res_dir):
    """TODO: Make a function systematically set up parameters for all runs."""
    if type(directory) is not str:
        raise ValueError('directory must be a str')
    if type(files) is not list:
        raise ValueError('files must be a list')

    print('This sequence file contains a Kallisto_Info file\
            and cannot be processed at the moment.')
    return '# {0} could not be processed'.format(res_dir), ''


def implicit_kallisto(directory, files, res_dir):
    """A function to write a Kallisto command with standard parameter setup."""
    if type(directory) is not str:
        raise ValueError('directory must be a str')
    if type(files) is not list:
        raise ValueError('files must be a list')

    # parts of each kallisto statement

    # information
    info = '# kallisto command for {0}'.format(directory)
    # transcript file location:
    k_head = 'kallisto quant -i input/transcripts.idx -o '

    # output file location
    k_output = 'input/kallisto_all/' + res_dir + '/kallisto '
    # parameter info:
    k_params = '--single -s {0} -l {1} -b {2} -t {3}'.format(sigma, length,
                                                             btstrp, thrds)

    # what files to use:
    k_files = ''
    # go through each file and add it to the command
    # unless it's a SampleSheet.csv file, in which
    # case you should ignore it.
    for y in files:
        if y != 'SampleSheet.csv':
            if directory[:3] == '../':
                d = directory[3:]
            else:
                d = directory[:]
            k_files += ' ' + d + '/' + y
    # all together now:
    kallisto = k_head + k_output + k_params + k_files + ';'
    return info, kallisto


def random_kallisto(directory, files, res_dir, n=1):
    """A function to write a Kallisto command with standard parameter setup."""
    if type(directory) is not str:
        raise ValueError('directory must be a str')
    if type(files) is not list:
        raise ValueError('files must be a list')

    # parts of each kallisto statement

    # information
    info = '# kallisto command for {0}'.format(directory)
    # transcript file location:
    k_head = 'kallisto quant -i input/transcripts.idx -o '

    # output file location
    k_output = 'input/kallisto_all/' + res_dir + '/kallisto '
    # parameter info:
    k_params = '--single -s {0} -l {1} -b {2} -t {3}'.format(sigma, length,
                                                             btstrp, thrds)

    # what files to use:
    k_files = ''

    # remove the 'SampleSheet.csv' entry from files:
    if 'SampleSheet.csv' in files:
        files.remove('SampleSheet.csv')

    # randomly select n files:
    selected = np.random.choice(files, n)
    for i, y in enumerate(selected):
        if y != 'SampleSheet.csv':
            if directory[:3] == '../':
                d = directory[3:]
            else:
                d = directory[:]
            k_files += ' ' + d + '/' + y
    # all together now:
    kallisto = k_head + k_output + k_params + k_files + ';'
    return info, kallisto


def walk_seq_directories(directory, output_dir='../input/kallisto_all/'):
    """
    Given a directory, find all the rna-seq folders and make kallisto commands.

    Directory format is predefined and must follow my rules.
    """
    kallisto = ''
    # directory contains all the projects, walk through it:
    for x in os.walk(directory):
        # first directory is always parent
        # if it's not the parent, move forward:
        if x[0] != directory:
            # cut the head off and get the project name:
            res_dir = x[0][len(directory)+1:]

            # if this project has attributes explicitly written in
            # use those parameter specs:
            if 'Kallisto_Info.csv' in x[2]:
                info, command = explicit_kallisto(x[0], x[2], res_dir)
                continue

            # otherwise, best guesses:
            info, command = implicit_kallisto(x[0], x[2], res_dir)
            kallisto += info + '\n' + command + '\n'

            if not os.path.exists(output_dir + res_dir):
                os.makedirs(output_dir + res_dir)
    return kallisto

with open('../kallisto_commands.sh', 'w') as f:
    f.write('#!/bin/bash\n')
    f.write('# make transcript index\n')
    s1 = 'kallisto index -i '
    s2 = 'input/transcripts.idx input/c_elegans_WBcel235.rel79.cdna.all.fa;\n'
    f.write(s1+s2)
    kallisto = walk_seq_directories(directory, output_dir)
    f.write(kallisto)
