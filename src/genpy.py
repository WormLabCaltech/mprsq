"""
A collection of useful functions for genetic analysis of transcriptomic data.

Author: David Angeles-Albores
contact: dangeles@caltech.edu

Citation: Forthcoming
"""
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import sklearn.decomposition
from scipy import stats

#################################################################
# PCA functions:
#################################################################


def make_matrix(hash_of_df, map_df, map_col, sel_col):
    """
    A function to concatenate a set of columns into a matrix.

    hash_of_df - a hash of dfs each containing a col. 'sel column'.
    map_df - the df that contains the keys for each df in hash_of_df
    map_col - the column that contains the keys for hash_of_df in map_df
    sel_col - the column to concatenate.
    """
    array = np.array([])
    for name in map_df[map_col].unique():
        x = hash_of_df[name][sel_col].values
        if len(array) == 0:
            array = np.array(x)
        else:
            array = np.vstack((array, x))
    return array


def pca(matrix, v=True):
    """
    A thin wrapper around sklearn.decomposition.PCA.

    matrix - the matrix to do PCA on
    v - verbose mode (prints out a message)
    alpha - percent of the variance that should be explained by n components

    returns:
    sklearn_pca - the sklearn fit model of pca
    n - the number of components that explain alpha% of the var.
    """
    sklearn_pca = sklearn.decomposition.PCA()
    sklearn_pca.fit(matrix)
    n = np.max(np.where(np.cumsum(
                        sklearn_pca.explained_variance_ratio_) < 0.9))

    if v:
        message = 'The first {} principal components explain >=90% of the data'
        print(message.format(n))

    return sklearn_pca, n


def tidy_pca(array, n):
    """
    A function to tidy up a sklearn pca into a tidy df.

    array - array to to feed into pca
    n - components that you desire to keep.

    Returns:
    df_nD - a pandas df
    """
    # Perform the PCA again retaining only the top 'where' components
    sklearn_pca = sklearn.decomposition.PCA(n_components=n)
    sklearn_pca.fit(array)

    # Project the data into this 'where'D space and convert it
    # back to a tidy dataframe
    cols = []
    for i in np.arange(1, n+1):
        cols += ['PCA{0}'.format(i)]

    df_nD = pd.DataFrame(sklearn_pca.transform(array),
                         columns=cols)
    return df_nD

#################################################################
# Probabilistic Genetic Analysis
#################################################################


def find(x, y, col='target_id'):
    """
    Thin wrapper around isin() function.

    Given a dataframe 'x' with a column 'col',
    and a set of target_id names 'y', find the entries
    in x that are present in y.

    Params:
    x - dataframe
    y - pandas series or list-like

    output:
    x[col].isin(y) - a pandas series object
    """
    return x[col].isin(y)


def overlap(x, y, q=0.1, col='target_id', qval='qval'):
    """
    Given two dataframes and a qvalue, returns the overlap between sig. genes.

    Both x and y must contain columns 'col' and 'qval'.

    Params:
    x - pandas df
    y - pandas df
    q - qvalue threshold
    col - column containing gene IDs (must be present in both x and y)
    qval - column containing qvalues (must be present in both x and y)

    output:
    overlap - a numpy ndarray
    """
    x_sig = x[x[qval] < q]
    y_sig = y[y[qval] < q]

    # overlap between x, y
    overlap = x_sig[x_sig[col].isin(y_sig[col].values)][col].values
    return overlap


def polarize(test_df, ref_df, ids, sign, col='target_id', change='b'):
    """
    Return the list of polarized genes between two dataframes.

    Given two dataframes, a list of overlapping genes and a sign,
    returns the list of genes that satisfy the sign condition.

    Params:
    test_df - pandas dataframe
    ref_df - pandas dataframe
    ids - list of ids that overlap between test_df and ref_df
          (from 'overlap' function)
    sign - one of '+', '-', 'p'. '+' means that genes must be positively
           correlated between test_df and ref_df, '-' requires anticorrelation,
           and 'p' means that they must have 0 correlation (function not yet
           implemented)
    col - name of column containing the gene names
    change - name of column containing fold change, regression value or some
             related metric for change relative to the WT. Absolute value
             doesn't matter, but it must not be absolute change!

    Output:
    g_overlap - numpy.ndarray of overlapping gene names
    """
    test = test_df[find(test_df, ids, col=col)].copy()
    ref = ref_df[find(ref_df, ids, col=col)]

    if ~(test[col].values == ref[col].values).all():
        test.sort_values(col, inplace=True)
        ref.sort_values(col, inplace=True)

    test['polarity'] = test[change]*ref[change]

    if sign == '+':
        g_overlap = test[test.polarity > 0][col].values
    elif sign == '-':
        g_overlap = test[test.polarity < 0][col].values
    elif sign == 'p':
        g_overlap = 0
        print('Unfinished business!')
    else:
        raise ValueError('sign must be one of +, -, or p')
    return g_overlap


def overlap_prob(ref_df, test_df, sign, q=0.1, qval='qval', genes='target_id',
                 change='b'):
    """
    A function to calculate the hypergeom probability of a particular overlap.

    Note: this function can't distinguish the direction of the network, so
    a ---> b and b ----> a are the same for the purposes of this function.

    Params:
    ref_df - reference dataframe (must be df with largest # of sig. genes)
    test_df - test dataframe
    sign - "+", "-" or "p". + means activating relationship, "-" means
            inhibitory relationship, and "p" is otherwise (not yet implemented)
    qval - name of the column containing qvalues
    genes - name of column containing gene IDs

    Output:
    prob_overlap - probability of overlap as measured by hypergeometric fnctn
    overlap_fraction - fraction of overlap genes between dfs (A and B)/(A U B)
    expected_frac - expected fraction of overlap for these gene sets
    polarized_ids - gene ids of the overlapping set
    """
    # find overlap:
    ovrlp_ids = overlap(ref_df, test_df, q=q, col=genes, qval=qval)

    # sig genes:
    ref_sig_df = ref_df[ref_df[qval] < q].copy()
    test_sig_df = test_df[test_df[qval] < q].copy()

    # for overlapping ids, check what number of them satisfy the condition
    # specified by 'sign'
    # call them polarized bc they have a sign
    polarized_ids = polarize(test_df, ref_df, ovrlp_ids, sign,
                             col=genes, change=change)

    # turn this into a scalar:
    n_polar = len(polarized_ids)

    # genes that changed significantly in either comparison (non-redundant)
    g_sig_ref = list(set(pd.Series.tolist(ref_sig_df[genes].values)))
    g_sig_test = list(set(pd.Series.tolist(test_sig_df[genes].values)))

    # convert those lists to scalars:
    n_sig_ref = len(g_sig_ref)
    n_sig_test = len(g_sig_test)

    # to calculate the overlap fraction, we need to know how
    # many genes in A U B (A or B)
    total_genes_changed = len(list(set(g_sig_ref+g_sig_test)))

    # total genes measured in this experiment pair:
    genes_ref = pd.Series.tolist(ref_df[genes].values)
    genes_test = pd.Series.tolist(test_df[genes].values)
    total_genes_measured = len(list(set(genes_ref + genes_test)))

    # calculate prob of overlap:
    prob_overlap = stats.hypergeom.cdf(n_polar, total_genes_measured,
                                       n_sig_ref, n_sig_test)

    # overlap fraction
    overlap_fraction = n_polar/total_genes_changed

    # expected:
    expected = stats.hypergeom.mean(total_genes_measured,
                                    n_sig_ref, n_sig_test)
    expected_frac = expected/total_genes_changed

    return prob_overlap, overlap_fraction, expected_frac, polarized_ids


def a_interacts_b(a_df, b_df, q=0.1, sign='+', qval='qval', genes='target_id',
                  change='b'):
    """
    A function to test whether a interacts with b in some way.

    a --> b or b --> a are the same thing for this function.

    Params:
    a_df - pandas dataframe
    b_df - pandas dataframe
    q - qvalue cutoff
    sign - '+', '-', 'p'. + tests activation, - tests inhibition, p not yet
            implemented.
    qval - name of column containing qvalues
    genes - name of column containing gene IDs

    Ouput:
    overlap_p
    overlap_f
    expected
    ids_overlap
    """
    a_sig_genes = (a_df.qval < q)
    b_sig_genes = (b_df.qval < q)

    # Check which has more upregulated genes
    if len(a_df[a_sig_genes]) > len(b_df[b_sig_genes]):
        results = overlap_prob(a_df, b_df, sign, q=q, qval=qval, genes=genes,
                               change=change)
    else:
        results = overlap_prob(b_df, a_df, sign, q=q, qval=qval, genes=genes,
                               change=change)
    overlap_p, overlap_f, expected_frac, ids_overlap = results
    return overlap_p, overlap_f, expected_frac, ids_overlap


def single_mutant_analysis(single_mutants, df_hash, genes='target_id',
                           analysis='spearmanr', qval='qval', q=0.1,
                           change='b', alpha=10**-4):
    """
    A function to perform single mutant analyses on a dataset of mutants.

    Params:
    single_mutants - list-like iterable containing the mutants to be correlated
    df_hash - hash with keys 'single_mutants'. contains the data to be used
    genes - name of the column that contains the genes ids
    analysis - one of spearmanr or interaction, defines the kind of analysis
               to be implemented
    qval - name of the column that contains the qvalues for each gene
    q - q value threshold for significance
    change - name of the column that contains the fold-changes or regression
            values for each gene
    alpha - significance value for spearmanr correlation

    Outputs:
    res_dict - a hash containing the results of the analysis.
    """
    def lind(x, col=qval):
        return (x[col] < 0.1)

    acceptable = ['spearmanr', 'interaction']
    if analysis not in acceptable:
        raise ValueError('analysis must be one of spearmanr or interaction')
    s = len(single_mutants)

    if analysis == 'spearmanr':

        rho_matrix = np.empty(shape=(s, s))
        res_dict = {'rho': rho_matrix}
    else:
        prob_plus_matrix = np.zeros(shape=(s, s))
        overlap_plus_matrix = np.zeros(shape=(s, s))
        expected_plus_matrix = np.zeros(shape=(s, s))

        prob_minus_matrix = np.zeros(shape=(s, s))
        overlap_minus_matrix = np.zeros(shape=(s, s))
        expected_minus_matrix = np.zeros(shape=(s, s))

        res_dict = {
            'prob_pos': prob_plus_matrix,
            'prob_minus': prob_minus_matrix,
            'overlap_pos': overlap_plus_matrix,
            'overlap_minus': overlap_minus_matrix,
            'expected_pos': expected_plus_matrix,
            'expected_minus': expected_minus_matrix,
            'ids': {},
        }
    l = 0
    ids = {}
    for i in single_mutants:
        m = 0
        for j in single_mutants:
            x = df_hash[i]
            y = df_hash[j]

            if analysis == 'spearmanr':
                # find overlap in stat.sig.genes between both lists:
                ovx = x[lind(x)]
                ovy = y[lind(y) & y[genes].isin(ovx[genes])]
                ovx = x[lind(x) & x[genes].isin(ovy[genes])]

                # spearman analysis
                rho = stats.spearmanr(ovx[change], ovy.b)
                if rho[1] < alpha:
                    res_dict['rho'][l, m] = rho[0]
                else:
                    res_dict['rho'][l, m] = 0
            elif analysis == 'interaction':
                # store the results from genpy.a_interacts_b in an array
                # called results, but remember it has 4 elements:
                # overlap prob, overlap frac, expected frac, ids overlapped
                results = a_interacts_b(x, y, sign='+', q=q,
                                        qval=qval, genes=genes, change=change)
                results2 = a_interacts_b(x, y, sign='-', q=q, qval=qval,
                                         genes=genes, change=change)

                # artificially set i,i entries for overlap fraction to zero,
                # this allows better discrimination of interactions for
                # heatmaps
                if i == j:
                    res_dict['overlap_pos'][l, m] = 0
                    res_dict['overlap_minus'][l, m] = 0
                else:
                    res_dict['overlap_pos'][l, m] = results[1]
                    res_dict['overlap_minus'][l, m] = results2[1]

                res_dict['prob_pos'][l, m] = results[0]
                res_dict['prob_minus'][l, m] = results2[0]
                res_dict['expected_pos'][l, m] = results[2]
                res_dict['expected_minus'][l, m] = results2[2]

                ids[('plus', i, j)] = results[3]
                ids[('minus', i, j)] = results2[3]

            m += 1
        l += 1

    if analysis == 'interaction':
        res_dict['ids'] = ids

    return res_dict


def double_mutant_corr_analysis(double_muts, df_hash, genes='target_id',
                                change='b', q=0.1, qval='qval'):
    """
    """
    def lind(x, col=qval):
        return (x[col] < 0.1)

    rho_matrix_doubles = np.zeros(shape=(2, 4))

    l = 0
    for key in double_muts:
        m = l*2
        for j in double_muts[key]:
            x = df_hash[key]
            y = df_hash[j]

            ovx = x[lind(x)]
            ovy = y[lind(y) & y[genes].isin(ovx[genes])]
            ovx = x[lind(x) & x[genes].isin(ovy[genes])]

            rho = stats.spearmanr(ovx[change], ovy[change])
            rho_matrix_doubles[l, m] = rho[0]
            print(key, j, '{0:.2g}'.format(rho[0]))
            m += 1
        l += 1

    return rho_matrix_doubles


def double_mutant_analysis(double_mutants, df_hash, genes='target_id',
                           qval='qval', q=0.1, change='b'):
    """
    """
    def lind(x, col=qval):
        return (x[col] < 0.1)

    res_dict = {}

    for number, i in enumerate(double_mutants):
        j, k = double_mutants[i]

        x = df_hash[i]
        y = df_hash[j]
        z = df_hash[k]

        # store the results from genpy.a_interacts_b in an array
        # called results, but remember it has 4 elements:
        # overlap prob, overlap frac, expected frac, ids overlapped
        results_xy = a_interacts_b(x, y, sign='+', q=q, qval=qval,
                                   genes=genes, change=change)
        results_xz = a_interacts_b(x, z, sign='+', q=q, qval=qval,
                                   genes=genes, change=change)

        # results_xy2 = a_interacts_b(x, y, sign='-', q=q, qval=qval,
        #                             genes=genes, change=change)
        # results_xz2 = a_interacts_b(x, z, sign='-', q=q, qval=qval,
        #                             genes=genes, change=change)

        # artificially set the i,i entries for overlap fraction to zero,
        # this allows better discrimination of interactions for
        # heatmaps
        # print(i, results_xy[0])
        # print(i, results_xz[0])
        # print(i, results_xy2[0])
        # print(i, results_xz2[0])
        # res_dict[(i, 'prob_pos')] = results_xy[0]*results_xz[0]
        # res_dict[(i, 'prob_minus')] = results_xy2[0]*results_xz2[0]

        # test epistasis:
        # log (OR) of probabilities...
        pxy = results_xy[0]
        pxz = results_xz[0]
        res_dict[(i, 'epistasis_odds')] = np.log(pxy) - np.log(pxz)
        res_dict[(i, 'epistasis_magnitude')] = np.maximum(pxy, pxz)

    return res_dict
