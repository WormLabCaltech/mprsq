"""
A script to process RNA-seq fastas into pseudoaligned reads using Kallisto.

authors: David Angeles-Albores and Paul W. Sternberg
contact: dangeles at caltech edu
"""
import os
import numpy as np
import argparse
import pandas as pd
import shutil

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

kallisto_loc = '../input/kallisto_all/'
genmap = pd.read_csv('../input/library_genotype_mapping.txt', comment='#')
genmap.genotype = genmap.genotype.apply(str)
# make sure everything is always in lowercase
genmap.genotype = genmap.genotype.apply(str.lower)


# Make all possible combinations of WT, X
combs = []
for gene in genmap.genotype.unique():
    if gene != 'wt':
        combs += [['WT', gene]]

# Make all the folders required for sleuth processing
sleuth_loc = '../sleuth/'

if not os.path.exists(sleuth_loc):
    os.makedirs(sleuth_loc)

WTnames = genmap[genmap.genotype=='wt'].project_name.values

# For each combination, make a folder
for comb in combs:
    current = sleuth_loc + comb[0]+'_'+comb[1]

    if not os.path.exists(current):
        os.makedirs(current)


    # copy the right files into the new directory
    # inside a folder called results
    def copy_cat(src_folder, dst_folder, names):
        """
        A function that copies a set of directories from one place to another.
        """
        for name in names:
            print('The following file was created:', dst_folder+name)
            shutil.copytree(src_folder + name, dst_folder + name)

    # copy WT files into the new directory
    copy_cat(kallisto_loc, current+'/results/', WTnames)

    # copy the MT files into the new directory
    MTnames = genmap[genmap.genotype == comb[1]].project_name.values
    copy_cat(kallisto_loc, current+'/results/', MTnames)


def matrix_design(name, factor, df, a, b, directory):
    """
    A function that makes the matrix design file for sleuth.

    This function can only make single factor design matrices.

    This function requires a folder 'results' to exist within
    'directory', and the 'results' folder in turn must contain
    files that are named exactly the same as in the dataframe.

    name - a string
    factor - list of factors to list in columns
    df - a dataframe containing the list of project names and the value for each factor
    i.e. sample1, wt, pathogen_exposed.
    a, b - conditions to slice the df with, i.e: a=WT, b=MT1
    directory - the directory address to place file in folder is in.
    """
    with open(directory + name, 'w') as f:
        f.write('# Sleuth design matrix for {0}-{1}\n'.format(a, b))
        f.write('experiment {0}\n'.format(factor))

        # walk through the results directory and get each folder name
        # write in the factor value by looking in the dataframe
        names = next(os.walk(directory+'/results/'))[1]
        for name in names:
            fval = df[df.project_name == name][factor].values[0]

            # add a if fval is WT or z otherwise
            # this is to ensure sleuth does
            # the regression as WT --> MT
            # but since sleuth only works alphabetically
            # simply stating WT --> MT doesn't work
            if fval == 'wt':
                fval = 'a' + fval
            else:
                fval = 'z' + fval
            line = name + ' ' + fval + '\n'
            f.write(line)


# Now make the matrix for each combination
# I made this separately from the above if loop
# because R is stupid and wants the files in this
# folder to be in the same order as they are
# listed in the matrix design file....
for comb in combs:
    current = sleuth_loc + comb[0]+'_'+comb[1] + '/'

    # write a file called rna_seq_info for each combination
    matrix_design('rna_seq_info.txt', 'genotype', genmap,
                  comb[0], comb[1], current)
