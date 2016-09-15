"""
A script to process RNA-seq fastas into pseudoaligned reads using Kallisto.

authors: David Angeles-Albores and Paul W. Sternberg
contact: dangeles at caltech edu
"""
import os
import argparse
import pandas as pd
import shutil

# get the read number from terminal
parser = argparse.ArgumentParser()
parser.add_argument("-k", "--kallisto_dir", nargs='?', type=str,
                    help="Directory where kallisto files are kept.")
parser.add_argument("-sl", "--sleuth_dir", nargs='?', type=str,
                    help="Directory to place output in.")
parser.add_argument("-f", "--fname", nargs='?', type=str,
                    help="Name of the sh script to execute.")
parser.add_argument("-g", "--genmap", nargs='?', type=str,
                    help="Location of the genmap file.")
args = parser.parse_args()

# params
if args.kallisto_dir:
    kallisto_dir = args.kallisto_dir
else:
    kallisto_dir = '../input/kallisto_all/'

if args.sleuth_dir:
    sleuth_dir = args.sleuth_dir
else:
    sleuth_dir = '../sleuth/'

if args.fname:
    fname = args.fname
else:
    fname = sleuth_dir + 'sleuth_commands.sh'

if args.genmap:
    genmap = args.genmap
else:
    genmap = '../input/library_genotype_mapping.txt'

if kallisto_dir[-1:] != '/':
    kallisto_dir += '/'
if sleuth_dir[-1:] != '/':
    sleuth_dir += '/'

# load genotype information
genmap = pd.read_csv(genmap, comment='#')
genmap.genotype = genmap.genotype.apply(str)
# make sure everything is always in lowercase
genmap.genotype = genmap.genotype.apply(str.lower)


def matrix_design(name, factor, df, a, b, directory):
    """
    A function that makes the matrix design file for sleuth.

    This function can only make single factor design matrices.

    This function requires a folder 'results' to exist within
    'directory', and the 'results' folder in turn must contain
    files that are named exactly the same as in the dataframe.

    name - a string
    factor - list of factors to list in columns
    df - dataframe with the list of project names and the value for each factor
    i.e. sample1, wt, mt.
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


def sleuth_analysis(main, directory, genovar):
    """A function to write the diff_exp_analyzer batch command."""
    s = 'Rscript {0}diff_exp_analyzer.R -d {1} --genovar {2}'
    heart = s.format(main, main+directory, genovar)
    return heart


def walk_sleuth_directories(directory):
    """Given a directory, generate all necessary sleuth commands."""
    sleuth = ''
    # directory contains all the projects, walk through it:
    current, dirs, files = next(os.walk(directory))
    for d in dirs:
        # genovar always begins with a z:
        genovar = 'z' + d[-1:]
        chmod = 'chmod 755 {0}{1}/kallisto/abundance.h5\n'.format(sleuth_dir,
                                                                  d)
        message = chmod + '# Sleuth analysis command for {0}\n'.format(d)
        command = sleuth_analysis(sleuth_dir, d, genovar) + '\n'
        sleuth += message
        sleuth += command
    return sleuth

###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################

# Make all possible combinations of WT, X
combs = []
for gene in genmap.genotype.unique():
    if gene != 'wt':
        combs += [['WT', gene]]

# make the sleuth directory
if not os.path.exists(sleuth_dir):
    os.makedirs(sleuth_dir)

# sort the groups by batches, then do the comparisons by batch
grouped = genmap.groupby('batch')


# do the comparison by batches
for name, group in grouped:
    WTnames = group[group.genotype == 'wt'].project_name.values
    # For each combination, make a folder
    for comb in combs:
        current = sleuth_dir + comb[0]+'_'+comb[1]
        MTnames = group[group.genotype == comb[1]].project_name.values
        if len(MTnames) == 0:
            continue
        # make the 'current' directory if it doesn't exist already
        if not os.path.exists(current):
            os.makedirs(current)

        # copy the right files into the new directory
        # inside a folder called results
        def copy_cat(src_folder, dst_folder, names):
            """A function to copy directories from one place to another."""
            for name in names:
                shutil.copytree(src_folder + name, dst_folder + name)

        # copy WT files into the new directory
        copy_cat(kallisto_dir, current+'/results/', WTnames)

        # copy the MT files into the new directory
        copy_cat(kallisto_dir, current+'/results/', MTnames)

# Now make the matrix for each combination
# I made this separately from the above if loop
# because R is stupid and wants the files in this
# folder to be in the same order as they are
# listed in the matrix design file....
for comb in combs:
    current = sleuth_dir + comb[0]+'_'+comb[1] + '/'

    # write a file called rna_seq_info for each combination
    matrix_design('rna_seq_info.txt', 'genotype', genmap,
                  comb[0], comb[1], current)

# Now write the sleuth commands!
# sequences:
# analysis = next(os.walk(sleuth_dir))[1]
# TODO
if sleuth_dir[-1:] != '/':
    sleuth_dir += '/'

shutil.copyfile('sleuth/diff_exp_analyzer.R',
                sleuth_dir + 'diff_exp_analyzer.R')


with open(sleuth_dir + fname, 'w') as f:
    f.write('#!/bin/bash\n')
    f.write('# Bash commands for diff. expression analysis using Sleuth.\n')
    sleuth_command = walk_sleuth_directories(sleuth_dir)
    f.write(sleuth_command)
