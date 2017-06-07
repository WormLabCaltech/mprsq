# important stuff:
import os
import pandas as pd
import numpy as np

# TEA and morgan
import tissue_enrichment_analysis as tea
import morgan as morgan

# Graphics
import matplotlib as mpl
import matplotlib.ticker as plticker
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patheffects as path_effects

import genpy

from matplotlib import rc
rc('text', usetex=True)
rc('text.latex', preamble=r'\usepackage{cmbright}')
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})

# JB's favorite Seaborn settings for notebooks
rc = {'lines.linewidth': 2,
      'axes.labelsize': 18,
      'axes.titlesize': 18,
      'axes.facecolor': 'DFDFE5'}
sns.set_context('notebook', rc=rc)
sns.set_style("dark")

ft = 35  # title fontsize


def pathify(title, xlabel, ylabel, xticks=True, yticks=True, **kwargs):
    """
    A function to pathify the labels, titles and ticks in a plot.

    Params:
    """
    labelsize = kwargs.pop('labelsize', 20)
    titlesize = kwargs.pop('titlesize', 25)

    # make the labels and title into paths
    effect = [path_effects.Normal()]
    plt.ylabel(ylabel,
               fontsize=labelsize).set_path_effects(effect)
    plt.xlabel(xlabel,
               fontsize=labelsize).set_path_effects(effect)
    plt.title(title,
              fontsize=titlesize).set_path_effects(effect)

    ax = plt.gca()
    # go through each xtick or ytick and make
    # it a path if user specified to do so.
    if xticks is True:
        for i, label in enumerate(ax.get_xticklabels()):
            ax.get_xticklabels()[i].set_path_effects(effect)
    if yticks is True:
        for i, label in enumerate(ax.get_yticklabels()):
            ax.get_yticklabels()[i].set_path_effects(effect)


def plot_by_term(term, df, kind='go', q=0.1, swarm=True,
                 x='genotype', y='b', gene='ens_gene'):
    """
    Plot ontology terms by a given column.

    Params:
    term - term to look for in melted_df
    df - a tidy dataframe with columns x and y
    kind - the ontology to use
    q - q-value for statistical significance. defaults to 0.1
    swarm - if True, plots a swarmplot. Else, plots a violinplot.
    x - column to plot on x axis
    y - column to plot on y axis
    gene - column in the given df where gene WBIDs are provided

    Output:
    ax - an axis object containing a graph
    genes - a list of genes obtained from the melted df
    """
    if type(kind) is not str:
        raise ValueError('`kind` variable must be a string.')

    if kind.lower() not in ['tissue', 'phenotype', 'go']:
        raise ValueError('`kind` must be one of `tissue`, `phenotype` or `go`')

    if type(term) is not str:
        raise ValueError('`term` must be a string.')

    if kind.lower() == 'tissue':
        onto_df = tea.fetch_dictionary()
    elif kind.lower() == 'phenotype':
        onto_df = pd.read_csv('../input/phenotype_ontology.csv')
    else:
        onto_df = pd.read_csv('../input/go_dictionary.csv')

    # melt the df:
    melted_df = pd.melt(onto_df, id_vars='wbid', var_name='term',
                        value_name='expressed')
    melted_df = melted_df[melted_df.expressed == 1]

    # warnings and bells:
    sel = melted_df.term.str.contains(term)
    if len(melted_df[sel].term.unique()) > 1:
        print('Warning: Provided term matches more than one ontological term.')

    genes = melted_df[sel].wbid

    if len(genes) == 0:
        raise ValueError('Provided term is not in ontology dictionary')

    ind = (df.qval < q) & (df[gene].isin(genes))
    fig, ax = plt.subplots()
    if swarm:
        ax = sns.swarmplot(x='genotype', y='b', data=df[ind])
    else:
        ax = sns.violinplot(x='genotype', y='b', data=df[ind])
    return ax, genes
