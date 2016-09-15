"""
A script to investigate how low we can go with the downsampling procedure.

authors: David Angeles-Albores and Paul W. Sternberg
contact: dangeles at caltech.. edu
"""
import pandas as pd
import os
import morgan as morgan
import argparse
import numpy as np

# get the read number from terminal
parser = argparse.ArgumentParser()
parser.add_argument("-n", "--number", type=int,
                    help="Number of files that were processed in this sample")
parser.add_argument("-k", "--kallisto_dir", nargs='?', type=str,
                    help="Directory where kallisto files are kept.")
parser.add_argument("-sl", "--sleuth_dir", nargs='?', type=str,
                    help="Directory to place output in.")

args = parser.parse_args()
if args.number:
    number = args.number
else:
    number = -1
if args.kallisto_dir:
    kallisto_dir = args.kallisto_dir
else:
    kallisto_dir = '../input/kallisto_all/'

if args.sleuth_dir:
    sleuth_dir = args.sleuth_dir
else:
    sleuth_dir = '../sleuth/'


# load parameters:
single_mutants = ['b', 'c', 'd', 'e', 'g']
double_mutants = {'a': 'bd', 'f': 'bc'}

# initialize thomas object and load all params
thomas = morgan.hunt('target_id', 'b', 'tpm', 'qval')
thomas.add_genmap('input/library_genotype_mapping.txt', comment='#')
thomas.add_single_mutant(single_mutants)
thomas.add_double_mutants(['a', 'f'], ['bd', 'bc'])
thomas.set_qval()

# Add the tpm files:
kallisto_loc = args.kallisto_dir
thomas.add_tpm(kallisto_loc, '/kallisto/abundance.tsv', '')

# Make all possible combinations of WT, X
combs = {}
for gene in thomas.genmap.genotype.unique():
    if gene != 'wt':
        combs[gene] = 'WT_'+gene+'/'

# load all the beta values for each genotype:
sleuth_loc = args.sleuth_dir
thomas.add_betas(sleuth_loc, 'betas.csv', combs)
thomas.filter_data(0, 0.1)

# extract relevant information
mat = np.matrix([['all', thomas.beta['a'].ens_gene.unique().shape[0],
                number, 'total_genes'],
                ['all', thomas.beta_filtered['a'].ens_gene.unique().shape[0],
                number, 'filtered_genes']])

for key in thomas.beta_filtered.keys():
    df = thomas.beta_filtered[key]
    sig = df[df.qval < 0.1].ens_gene.unique().shape[0]
    row1 = [key, df.ens_gene.unique().shape[0], number, 'total_genes']
    row2 = [key, sig, number, 'significant_genes']
    mat = np.vstack((mat, row1))
    mat = np.vstack((mat, row2))

df = pd.DataFrame(mat, columns=['pair', 'a_value', 'analysis', 'reads'])

# perform analyses
fast = morgan.brenner('single_mutant_analysis', thomas)
bayes = morgan.mcclintock('single_mutant_analysis', thomas, progress=False)


# super melt function
def supermelt(df, analysis, reads, col='test_value'):
    """A function to melt a dataframe into a predetermined format."""
    new = pd.melt(df, id_vars='corr_with', var_name='genotype', value_name=col)
    new = new[new[col].abs() > 0]
    new['pair'] = new.genotype + new.corr_with
    new = new[new.genotype != new.corr_with]
    new['analysis'] = analysis
    new['reads'] = reads
    return new[['pair', col, 'analysis', 'reads']]

# melt everything
rho = supermelt(fast.rho, 'spearman', number, 'a_value')
hyper_plus = supermelt(fast.hyper_plus, 'hyper_plus', number, 'a_value')
hyper_minus = supermelt(fast.hyper_minus, 'hyper_minus', number, 'a_value')
bayes_primary = supermelt(bayes.robust_slope, 'bayes_primary', number,
                          'a_value')

bayes.secondary_slope['corr_with'] = thomas.single_mutants
bayes_secondary = supermelt(bayes.secondary_slope, 'bayes_secondary', number,
                            'a_value')

# save it to a spreadsheet
if os.path.isfile('input/downsampling.csv') is False:
    rho.to_csv('input/downsampling.csv', index=False)
    hyper_plus.to_csv('input/downsampling.csv', mode='a',
                      index=False, header=False)
    hyper_minus.to_csv('input/downsampling.csv', mode='a',
                       index=False, header=False)
    bayes_primary.to_csv('input/downsampling.csv', mode='a',
                         index=False, header=False)
    bayes_secondary.to_csv('input/downsampling.csv', mode='a',
                           index=False, header=False)
    df.to_csv('input/downsampling.csv', mode='a',
              index=False, header=False)
else:
    rho.to_csv('input/downsampling.csv', mode='a',
               index=False, header=False)
    hyper_plus.to_csv('input/downsampling.csv', mode='a',
                      index=False, header=False)
    hyper_minus.to_csv('input/downsampling.csv', mode='a',
                       index=False, header=False)
    bayes_primary.to_csv('input/downsampling.csv', mode='a',
                         index=False, header=False)
    bayes_secondary.to_csv('input/downsampling.csv', mode='a',
                           index=False, header=False)
    df.to_csv('input/downsampling.csv', mode='a', index=False, header=False)
