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


# def make_matrix(hash_of_df, map_df, map_col, sel_col):
#     """
#     A function to concatenate columns from different dataframes into a matrix.
#
#     hash_of_df - a hash of dfs each containing a col. 'sel column'.
#     map_df - the df that contains the keys for each df in hash_of_df
#     map_col - the column that contains the keys for hash_of_df in map_df
#     sel_col - the column to concatenate.
#     """
#     array = np.array([])
#     for name in map_df[map_col].unique():
#         x = hash_of_df[name][sel_col].values
#         if len(array) == 0:
#             array = np.array(x)
#         else:
#             array = np.vstack((array, x))
#     return array


# def pca(matrix, v=True):
#     """
#     A thin wrapper around sklearn.decomposition.PCA.
#
#     matrix - the matrix to do PCA on
#     v - verbose mode (prints out a message)
#     alpha - percent of the variance that should be explained by n components
#
#     returns:
#     sklearn_pca - the sklearn fit model of pca
#     n - the number of components that explain alpha% of the var.
#     """
#     sklearn_pca = sklearn.decomposition.PCA()
#     sklearn_pca.fit(matrix)
#     n = np.max(np.where(np.cumsum(
#                         sklearn_pca.explained_variance_ratio_) < 0.9))
#
#     if v:
#         message = 'The first {} principal components explain >=90% of the data'
#         print(message.format(n))
#
#     return sklearn_pca, n


# def tidy_pca(array, n):
#     """
#     A function to tidy up a sklearn pca into a tidy df.
#
#     array - array to to feed into pca
#     n - components that you desire to keep.
#
#     Returns:
#     df_nD - a pandas df
#     """
#     # Perform the PCA again retaining only the top 'where' components
#     sklearn_pca = sklearn.decomposition.PCA(n_components=n)
#     sklearn_pca.fit(array)
#
#     # Project the data into this 'where'D space and convert it
#     # back to a tidy dataframe
#     cols = []
#     for i in np.arange(1, n+1):
#         cols += ['PCA{0}'.format(i)]
#
#     df_nD = pd.DataFrame(sklearn_pca.transform(array),
#                          columns=cols)
#     return df_nD

#################################################################
# Probabilistic Genetic Analysis
#################################################################


# def find(x, y, col='target_id'):
#     """
#     Thin wrapper around isin() function.
#
#     Given a dataframe 'x' with a column 'col',
#     and a set of target_id names 'y', find the entries
#     in x that are present in y.
#
#     Params:
#     x - dataframe
#     y - pandas series or list-like
#
#     output:
#     x[col].isin(y) - a pandas series object
#     """
#     return x[col].isin(y)


# def overlap(x, y, q=0.1, col='target_id', qval='qval'):
#     """
#     Given two dataframes and a qvalue, returns the overlap between sig. genes.
#
#     Both x and y must contain columns 'col' and 'qval'.
#
#     Params:
#     x - pandas df
#     y - pandas df
#     q - qvalue threshold
#     col - column containing gene IDs (must be present in both x and y)
#     qval - column containing qvalues (must be present in both x and y)
#
#     output:
#     overlap - a numpy ndarray
#     """
#     x_sig = x[x[qval] < q]
#     y_sig = y[y[qval] < q]
#
#     # overlap between x, y
#     overlap = x_sig[x_sig[col].isin(y_sig[col].values)][col].values
#     return overlap




# def overlap_prob(ref_df, test_df, sign, q=0.1, qval='qval', genes='target_id',
#                  change='b'):
#     """
#     A function to calculate the hypergeom probability of a particular overlap.
#
#     Note: this function can't distinguish the direction of the network, so
#     a ---> b and b ----> a are the same for the purposes of this function.
#
#     Params:
#     ref_df - reference dataframe (must be df with largest # of sig. genes)
#     test_df - test dataframe
#     sign - "+", "-" or "p". + means activating relationship, "-" means
#             inhibitory relationship, and "p" is otherwise (not yet implemented)
#     qval - name of the column containing qvalues
#     genes - name of column containing gene IDs
#
#     Output:
#     prob_overlap - probability of overlap as measured by hypergeometric fnctn
#     overlap_fraction - fraction of overlap genes between dfs (A and B)/(A U B)
#     expected_frac - expected fraction of overlap for these gene sets
#     polarized_ids - gene ids of the overlapping set
#     """
#     # find overlap:
#     ovrlp_ids = overlap(ref_df, test_df, q=q, col=genes, qval=qval)
#
#     # sig genes:
#     ref_sig_df = ref_df[ref_df[qval] < q].copy()
#     test_sig_df = test_df[test_df[qval] < q].copy()
#
#     # for overlapping ids, check what number of them satisfy the condition
#     # specified by 'sign'
#     # call them polarized bc they have a sign
#     polarized_ids = polarize(test_df, ref_df, ovrlp_ids, sign,
#                              col=genes, change=change)
#
#     # turn this into a scalar:
#     n_polar = len(polarized_ids)
#
#     # genes that changed significantly in either comparison (non-redundant)
#     g_sig_ref = list(set(pd.Series.tolist(ref_sig_df[genes].values)))
#     g_sig_test = list(set(pd.Series.tolist(test_sig_df[genes].values)))
#
#     # convert those lists to scalars:
#     n_sig_ref = len(g_sig_ref)
#     n_sig_test = len(g_sig_test)
#
#     # to calculate the overlap fraction, we need to know how
#     # many genes in A U B (A or B)
#     total_genes_changed = len(list(set(g_sig_ref+g_sig_test)))
#
#     # total genes measured in this experiment pair:
#     genes_ref = pd.Series.tolist(ref_df[genes].values)
#     genes_test = pd.Series.tolist(test_df[genes].values)
#     total_genes_measured = len(list(set(genes_ref + genes_test)))
#
#     # calculate prob of overlap:
#     prob_overlap = stats.hypergeom.cdf(n_polar, total_genes_measured,
#                                        n_sig_ref, n_sig_test)
#
#     # overlap fraction
#     overlap_fraction = n_polar/total_genes_changed
#
#     # expected:
#     expected = stats.hypergeom.mean(total_genes_measured,
#                                     n_sig_ref, n_sig_test)
#     expected_frac = expected/total_genes_changed
#
#     return prob_overlap, overlap_fraction, expected_frac, polarized_ids


# def a_interacts_b(a_df, b_df, q=0.1, sign='+', qval='qval', genes='target_id',
#                   change='b'):
#     """
#     A function to test whether a interacts with b in some way.
#
#     a --> b or b --> a are the same thing for this function.
#
#     Params:
#     a_df - pandas dataframe
#     b_df - pandas dataframe
#     q - qvalue cutoff
#     sign - '+', '-', 'p'. + tests activation, - tests inhibition, p not yet
#             implemented.
#     qval - name of column containing qvalues
#     genes - name of column containing gene IDs
#
#     Ouput:
#     overlap_p
#     overlap_f
#     expected
#     ids_overlap
#     """
#     a_sig_genes = (a_df.qval < q)
#     b_sig_genes = (b_df.qval < q)
#
#     # Check which has more upregulated genes
#     if len(a_df[a_sig_genes]) > len(b_df[b_sig_genes]):
#         results = overlap_prob(a_df, b_df, sign, q=q, qval=qval, genes=genes,
#                                change=change)
#     else:
#         results = overlap_prob(b_df, a_df, sign, q=q, qval=qval, genes=genes,
#                                change=change)
#     overlap_p, overlap_f, expected_frac, ids_overlap = results
#     return overlap_p, overlap_f, expected_frac, ids_overlap


# def single_mutant_analysis(single_mutants, df_hash, genes='target_id',
#                            analysis='spearmanr', qval='qval', q=0.1,
#                            change='b', alpha=10**-4):
#     """
#     A function to perform single mutant analyses on a dataset of mutants.
#
#     Params:
#     single_mutants - list-like iterable containing the mutants to be correlated
#     df_hash - hash with keys 'single_mutants'. contains the data to be used
#     genes - name of the column that contains the genes ids
#     analysis - one of spearmanr or interaction, defines the kind of analysis
#                to be implemented
#     qval - name of the column that contains the qvalues for each gene
#     q - q value threshold for significance
#     change - name of the column that contains the fold-changes or regression
#             values for each gene
#     alpha - significance value for spearmanr correlation
#
#     Outputs:
#     res_dict - a hash containing the results of the analysis.
#     """
#     def lind(x, col=qval):
#         return (x[col] < 0.1)
#
#     acceptable = ['spearmanr', 'interaction']
#     if analysis not in acceptable:
#         raise ValueError('analysis must be one of spearmanr or interaction')
#     s = len(single_mutants)
#
#     if analysis == 'spearmanr':
#
#         rho_matrix = np.empty(shape=(s, s))
#         res_dict = {'rho': rho_matrix}
#     else:
#         prob_plus_matrix = np.zeros(shape=(s, s))
#         overlap_plus_matrix = np.zeros(shape=(s, s))
#         expected_plus_matrix = np.zeros(shape=(s, s))
#
#         prob_minus_matrix = np.zeros(shape=(s, s))
#         overlap_minus_matrix = np.zeros(shape=(s, s))
#         expected_minus_matrix = np.zeros(shape=(s, s))
#
#         res_dict = {
#             'prob_pos': prob_plus_matrix,
#             'prob_minus': prob_minus_matrix,
#             'overlap_pos': overlap_plus_matrix,
#             'overlap_minus': overlap_minus_matrix,
#             'expected_pos': expected_plus_matrix,
#             'expected_minus': expected_minus_matrix,
#             'ids': {},
#         }
#     l = 0
#     ids = {}
#     for i in single_mutants:
#         m = 0
#         for j in single_mutants:
#             x = df_hash[i]
#             y = df_hash[j]
#
#             if analysis == 'spearmanr':
#                 # find overlap in stat.sig.genes between both lists:
#                 ovx = x[lind(x)]
#                 ovy = y[lind(y) & y[genes].isin(ovx[genes])]
#                 ovx = x[lind(x) & x[genes].isin(ovy[genes])]
#
#                 # spearman analysis
#                 rho = stats.spearmanr(ovx[change], ovy.b)
#                 if rho[1] < alpha:
#                     res_dict['rho'][l, m] = rho[0]
#                 else:
#                     res_dict['rho'][l, m] = 0
#             elif analysis == 'interaction':
#                 # store the results from genpy.a_interacts_b in an array
#                 # called results, but remember it has 4 elements:
#                 # overlap prob, overlap frac, expected frac, ids overlapped
#                 results = a_interacts_b(x, y, sign='+', q=q,
#                                         qval=qval, genes=genes, change=change)
#                 results2 = a_interacts_b(x, y, sign='-', q=q, qval=qval,
#                                          genes=genes, change=change)
#
#                 # artificially set i,i entries for overlap fraction to zero,
#                 # this allows better discrimination of interactions for
#                 # heatmaps
#                 if i == j:
#                     res_dict['overlap_pos'][l, m] = 0
#                     res_dict['overlap_minus'][l, m] = 0
#                 else:
#                     res_dict['overlap_pos'][l, m] = results[1]
#                     res_dict['overlap_minus'][l, m] = results2[1]
#
#                 res_dict['prob_pos'][l, m] = results[0]
#                 res_dict['prob_minus'][l, m] = results2[0]
#                 res_dict['expected_pos'][l, m] = results[2]
#                 res_dict['expected_minus'][l, m] = results2[2]
#
#                 ids[('plus', i, j)] = results[3]
#                 ids[('minus', i, j)] = results2[3]
#
#             m += 1
#         l += 1
#
#     if analysis == 'interaction':
#         res_dict['ids'] = ids
#
#     return res_dict
#
#
# def double_mutant_corr_analysis(double_muts, df_hash, genes='target_id',
#                                 change='b', q=0.1, qval='qval'):
#     """
#     """
#     def lind(x, col=qval):
#         return (x[col] < 0.1)
#
#     rho_matrix_doubles = np.zeros(shape=(2, 4))
#
#     l = 0
#     for key in double_muts:
#         m = l*2
#         for j in double_muts[key]:
#             x = df_hash[key]
#             y = df_hash[j]
#
#             ovx = x[lind(x)]
#             ovy = y[lind(y) & y[genes].isin(ovx[genes])]
#             ovx = x[lind(x) & x[genes].isin(ovy[genes])]
#
#             rho = stats.spearmanr(ovx[change], ovy[change])
#             rho_matrix_doubles[l, m] = rho[0]
#             print(key, j, '{0:.2g}'.format(rho[0]))
#             m += 1
#         l += 1
#
#     return rho_matrix_doubles
#

# def double_mutant_analysis(double_mutants, df_hash, genes='target_id',
#                            qval='qval', q=0.1, change='b'):
#     """
#     """
#     def lind(x, col=qval):
#         return (x[col] < 0.1)
#
#     res_dict = {}
#
#     for number, i in enumerate(double_mutants):
#         j, k = double_mutants[i]
#
#         x = df_hash[i]
#         y = df_hash[j]
#         z = df_hash[k]
#
#         # store the results from genpy.a_interacts_b in an array
#         # called results, but remember it has 4 elements:
#         # overlap prob, overlap frac, expected frac, ids overlapped
#         results_xy = a_interacts_b(x, y, sign='+', q=q, qval=qval,
#                                    genes=genes, change=change)
#         results_xz = a_interacts_b(x, z, sign='+', q=q, qval=qval,
#                                    genes=genes, change=change)
#
#         # results_xy2 = a_interacts_b(x, y, sign='-', q=q, qval=qval,
#         #                             genes=genes, change=change)
#         # results_xz2 = a_interacts_b(x, z, sign='-', q=q, qval=qval,
#         #                             genes=genes, change=change)
#
#         # artificially set the i,i entries for overlap fraction to zero,
#         # this allows better discrimination of interactions for
#         # heatmaps
#         # print(i, results_xy[0])
#         # print(i, results_xz[0])
#         # print(i, results_xy2[0])
#         # print(i, results_xz2[0])
#         # res_dict[(i, 'prob_pos')] = results_xy[0]*results_xz[0]
#         # res_dict[(i, 'prob_minus')] = results_xy2[0]*results_xz2[0]
#
#         # test epistasis:
#         # log (OR) of probabilities...
#         pxy = results_xy[0]
#         pxz = results_xz[0]
#         res_dict[(i, 'epistasis_odds')] = np.log(pxy) - np.log(pxz)
#         res_dict[(i, 'epistasis_magnitude')] = np.maximum(pxy, pxz)
#
#     return res_dict
