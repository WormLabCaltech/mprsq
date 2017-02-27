"""
A collection of useful functions for genetic analysis of transcriptomic data.

Author: David Angeles-Albores
contact: dangeles@caltech.edu

Citation: Forthcoming
"""
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import pymc3 as pm
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
import gvars
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

mpl.rcParams['xtick.labelsize'] = 18
mpl.rcParams['ytick.labelsize'] = 18
mpl.rcParams['legend.fontsize'] = 14

genvar = gvars.genvars()
#################################################################
# PCA functions:
#################################################################


def robust_regress(data):
    """A regression using a StudentT distribution instead of Gaussian."""
    with pm.Model() as model_robust:
        # pick the family to do regression with
        family = pm.glm.families.StudentT()
        # specify we want a generalized linear
        pm.glm.glm('y ~ x', data, family=family)
        # find the MAP as a good starting point
        start = pm.find_MAP()
        # simulate
        step = pm.NUTS(scaling=start)
        trace_robust = pm.sample(2000, step, progressbar=True)
        return trace_robust


# a function to rank order the data
def find_rank(df):
    """A function to find the rank values of a variable."""
    # make a copy of the dataframe, then sort it inplace
    d = df.copy()
    d.sort_values('b', inplace=True)
    # make a rank vector and append it to the sorted dataframe
    rank = np.linspace(0, len(d)-1, len(d))
    d['r'] = rank
    # sort by isoform name again and return the modified df
    d.sort_values('target_id', inplace=True)
    return d


# find inliers and outliers (see text description below)
def find_inliers(ovx, ovy, trace):
    """A function to identify inliers and outliers in a distribution."""
    # calculate some important numbers
    intercept, slope = trace.Intercept.mean(), trace.x.mean()
    distribution = ovy.r - intercept - ovx.r*slope
    mean, std = distribution.mean(), distribution.std()

    # variances add
    sigma = np.sqrt(trace.x.std()**2 + trace.Intercept.std()**2 + std**2)
    # find outliers
    sel = np.abs(distribution - mean)/sigma < 1.5

    # get the outliers and inliers
    distribution_outliers = distribution[~sel]

    # get the gene names of the outliers
    inverse = distribution_outliers + intercept + ovx.r*slope
    outliers = ovy[ovy.r.isin(inverse)].target_id

    return outliers


# Define a plotting function to plot only a triangular heat map
def tri_plot(matrix, xlabels, ylabels=[]):
    """Given a matrix, draw a triangle plot."""
    # Minimum and maximum for colormap
    vmin = matrix.min().min()
    vmax = np.max(matrix).max()

    # assume ylabels are the same as xlabels
    if len(ylabels) == 0:
        ylabels = xlabels

    # make the lower triangle of the matrix,
    # Also, remove the diagonal
    mask = np.zeros_like(matrix)
    mask[np.tril_indices_from(mask)] = True

    # draw and adjust xtick size
    with sns.axes_style("white"):
        ax = sns.heatmap(matrix, xticklabels=xlabels,
                         yticklabels=ylabels, cmap='RdBu',
                         mask=mask, square=True, vmin=vmin,
                         vmax=vmax)
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)

    return ax


def make_genetic_graph(robust_slope, id_vars='corr_with',
                       var_name='geno1', value_name='correl',
                       **kwargs):
    """A function to generate a graph of genetic relationships."""
    w = kwargs.pop('w', 2)

    # initialize an empty graph:
    G = nx.Graph()
    df = pd.melt(robust_slope, id_vars='corr_with', var_name='geno1',
                 value_name='correl')
    df = df[df.correl != 0]

    # add edges between nodes
    for g1 in df.geno1.unique():
        for g2 in df.corr_with.unique():
            if g1 == g2:
                continue

            # extract the correlation coefficient
            sel = (df.geno1 == g1) & (df.corr_with == g2)
            if len(df[sel]):
                r = df[sel].correl.values[0]
            # only add an edge if a correlation coefficient exists:
            if r:
                # add the edge
                G.add_edge(genvar.mapping[g1],
                           genvar.mapping[g2], weight=r)

    # parameterize the edge and width and color of the graphs
    elarge = [(u, v) for (u, v, d) in G.edges(data=True)]
    # set the width
    # width will be proportional to the log(correl/smallest_abs_corr)
    width = [w*np.log(np.abs(d['weight']/df.correl.abs().min()))  # log
             for (u, v, d) in G.edges(data=True)]  # width
    # extract the weights
    weights = [d['weight'] for (u, v, d) in G.edges(data=True)]

    return G, width, weights, elarge


# a qPCR barplot
def qPCR_plot(df, plotting, colors, **kwargs):
    """
    A function to make prettified qPCR barplots.

    Takes as entry a dataframe as output by qPCR_prep method
    Params:
    ------
    df -- df as output by qPCR_prep
    plotting -- a dictionary of plotting order for each gene;
                keys must be in dataframe column 'plotting_group'
    colors -- color to be used for each gene
    kwargs -- clustering - the name of a column within the dataframe,
              bars grouped within the same cluster are given the same color;
              plotting_group - a string that must be a column within the
              dataframe, bars belonging to the same plotting group are
              plotted adjacent to each other;
              alpha (transparency, float);
              q (stat. sig. thresh, float);
              save (string to save as)
              rotation;
              title

    outputs:
    -------
    a Seaborn barchart
    """
    clustering = kwargs.pop('clustering', 'ext_gene')
    plotting_group = kwargs.pop('plotting_group', 'genotype')
    alpha = kwargs.pop('alpha', 0.7)
    q = kwargs.pop('q', 0.1)
    rotation = kwargs.pop('rotation', 45)

    index = np.linspace(0, df[plotting_group].unique().shape[0]-1,
                        df[plotting_group].unique().shape[0])

    # error bars
    error_config = {'ecolor': '0.2'}

    # groupby gene name if it exists:
    grouped = df.groupby(clustering)

    bar_width = 1/(len(grouped)+1)

    # go through each gene
    for name, group in grouped:
        # figure out where each bar goes:
        if name not in plotting.keys():
            print(name, 'not in plotting.keys()')
            where = max(plotting.keys(),
                        key=lambda k: plotting[k])
            val = plotting[where]
            plotting[name] = val + 1

        add = plotting[name]*bar_width
        # figure out what color to give:
        if name in colors.keys():
            # add the bar:
            plt.bar(index + add + bar_width/2, group.b.values,
                    bar_width, alpha=alpha,
                    yerr=group.se_b.values,
                    error_kw=error_config, label=name,
                    color=colors[name])
        else:
            # add the bar but don't specify color
            plt.bar(index + add, group.b.values,
                    bar_width, alpha=alpha,
                    yerr=group.se_b.values,
                    error_kw=error_config, label=name)

        # significance threshold:
        sig = group.qval < q
        k = group[sig].order - 1

        # plot stars on top of stat. sig. results
        plt.plot(k + add + bar_width/2,
                 group[sig].b.values + group[sig].se_b.values + 0.20,
                 r'*', color='k')

    # shade in every other bar group for clarity:
    grouped2 = df.groupby(plotting_group)
    k = 0
    col = '#CFCFCF'

    ymin, ymax = plt.gca().get_ylim()
    for name, group in grouped2:
        if k % 2 == 0:
            xmin = k - bar_width*0.5
            xmax = k + bar_width*(len(grouped) + 0.5)

            plt.fill_between([xmin, xmax], ymin, color=col)
            plt.fill_between([xmin, xmax], ymax, color=col)
        k += 1

    # fix the xlims and tick params etc...
    if (k - 1) % 2 == 0:
        plt.xlim(0, xmax)
    else:
        plt.xlim(0, plt.gca().get_xlim()[1] - 3/2*bar_width)

    # plt.tick_params(axis='y', which='major', labelsize=18)
    plt.tick_params(axis='y', which='major')

    fancy_names = []
    fancy_names = []
    for label in df[plotting_group].unique():
        genename = df[df[plotting_group] == label].ens_gene.values[0]
        n_isoforms = df[df.ens_gene == genename].shape[0]/len(df.code.unique())
        if n_isoforms == 1:
            name = df[df[plotting_group] == label].ext_gene.values[0]
        else:
            name = label
        if r'\emph' not in name:
            fancy_names += [r'\emph{' + name + r'}']
        else:
            fancy_names += name

    plt.yticks(fontsize=18)
    plt.xticks(index + bar_width*len(df[clustering].unique())/2,
               fancy_names, rotation=rotation, fontsize=18,
               fontname='Helvetica')
    # plt.xticks(index + bar_width*len(df[clustering].unique())/2,
    #            fancy_names, rotation=rotation)

    # pathify(title, '', r'Regression Coefficient, $\beta$')
    plt.ylabel(r'Regression Coefficient, $\beta$', fontsize=20)
    legend = plt.legend(title='Genotype', loc=(1.02, 0.5), fontsize=18)
    plt.setp(legend.get_title(), fontsize=18)
    plt.ylim(ymin, ymax)
