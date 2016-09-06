
"""
Morgan.

authors: David Angeles-Albores and Paul W. Sternberg
contact: dangeles at caltech edu
"""
import pandas as pd
import warnings as wng
import numpy as np
from scipy import stats
import tissue_enrichment_analysis as tea
import pymc3 as pm
# import theano

###############################################################################
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
###############################################################################


class hunt(object):
    """morgan objects are used for genetic analysis using RNA-seq.

    Each genotype can be associated with two attributes: Read counts
    and log(fold-change). These attributes are provided in two
    different dataframes. If you provide a dataframe with fold-change
    (not log-foldchange) certain functions will not work correctly!

    Attributes:
    ------------------
    gene
    change
    counts
    qval
    q
    """

    def __init__(self, gene, change, counts, qval, q=0.1):
        """
        The initialize function.

        Params:
        gene
        change
        counts
        qval
        q
        """
        if not gene:
            raise ValueError('`gene` cannot be empty')
        if not change:
            raise ValueError('`change` cannot be empty')
        if not counts:
            raise ValueError('`counts` cannot be empty')
        if not qval:
            raise ValueError('`qval` cannot be empty')

        if type(gene) is not str:
            raise ValueError('`gene` must be a string')
        if type(change) is not str:
            raise ValueError('`change` must be a string')
        if type(counts) is not str:
            raise ValueError('`counts` must be a string')
        if type(qval) is not str:
            raise ValueError('`qval` must be a string')
        if type(q) is not float:
            raise ValueError('`q` must be a float')

        if q <= 0 or q >= 1:
            raise ValueError('`q` must be between 0 and 1')

        self.gene = gene
        self.change = change
        self.counts = counts
        self.qval = qval
        self.q = q
        self.single_mutants = []
        self.double_muts = {}

    def add_single_mutant(self, single):
        """
        Add a single mutant to the list.

        Params:
        single - str or listlike

        Note: ALL letter codes are coerced to lowercase!
        """
        if type(single) not in [str, list]:
            raise ValueError('`single` must be a str or list of strings')
        if type(single) is str:
            self.single_mutants += [single.lower()]

        if type(single) is list:
            self.single_mutants += [x.lower() for x in single]

        self.single_mutants = list(set(self.single_mutants))

    def add_double_mutants(self, lettercode, genotype):
        """
        A method that adds double mutants codes to a dictionary.

        Params:
        ---------
        lettercode - str or list, contains the code by which the double
        mutant will be referred to
        genotype - str or list, contains the genotype that lettercode refers
        to
        i.e.
        {a: bc} - the lettercode is a, and the genotype is b(minus)c(minus)
        Output:
        appends to the double mutant dictionary.
        """
        if type(lettercode) != type(genotype):
            raise ValueError('types of lettercode and genotype must match!')

        if type(lettercode) is not str:
            if len(lettercode) != len(genotype):
                raise ValueError('lengths of lettercode\
                                 and genotype must match!')

        if type(lettercode) is str:
            if lettercode.lower() in self.double_muts.keys():
                w = '{0} is already in string\
                     and was replaced'.format(lettercode.lower())
                wng.warn(w)

            self.double_muts[lettercode.lower()] = genotype.lower()
        else:
            for i, letter in enumerate(lettercode):
                if letter.lower() in self.double_muts.keys():
                    w = '{0} is already in string\
                        and was replaced'.format(letter.lower())
                    wng.warn('{0} is already in string!'.format(letter))

                self.double_muts[letter.lower()] = genotype[i].lower()

    def add_genmap(self, genmap_path, sep=',', comment='#'):
        """
        Add a genmap path to this object.

        The genmap file must have exactly two columns:
        project_name - the name of each RNA-seq run
        genotype - typically, each genotype has n replicates
        with n project_name's

        I.e.:
        run1,WT
        run2,WT
        run3,WT
        run4,mut
        run5,mut
        run6,mut

        Params:
        genmap_path - path (including filename) to genmap file
        sep - separator used to make genmap
        comment - if there are comments, marker used to define comments
        """
        self.genmap = pd.read_csv(genmap_path, sep=sep, comment=comment)

        if (self.genmap.columns != ['project_name', 'genotype']).all():
            raise ValueError('genmap is not in the right format!')

        self.genmap.genotype = self.genmap.genotype.apply(str)
        self.genmap.genotype = self.genmap.genotype.apply(str.lower)
        # make sure everything is always in lowercase

    def add_tpm(self, main_path, tpm_fname, folder='', sep='\t'):
        """
        Add tpm files.

        main_path - path where all the tpm files are kept
        tpm_fname - generic name of all tpm files (i.e., tpm.csv)
        folder - if there are any subfolders to get to tpm_fname, go here
        sep - separator used in tpm files

        i.e.:
        main_path -> genmap.project_name[0] -> folder -> tpm_fname
        main_path -> genmap.project_name[1] -> folder -> tpm_fname

        returns:
        self.tpm - a dictionary (project_name, df)
        """
        self.tpm = {}  # initialize an empty hash

        # get tpm for each project
        for prjct in self.genmap.project_name.unique():
            path = main_path + prjct + folder + tpm_fname
            self.tpm[prjct] = pd.read_csv(path, sep=sep)
            self.tpm[prjct].sort_values(self.gene, inplace=True)
            self.tpm[prjct].reset_index(drop=True, inplace=True)

    def add_betas(self, main_path, fc_fname, folders, sep=','):
        """
        Add fold change dfs.

        Params:
        -------------------------
        main_path - str, path to each processed read folder
        folders - dict, where keys are the genotypes and the values
        are the names of the folders the genotype is in
        fc_fname - str, standard name of the fold-change data
        sep - separators between columns

        Output:
        -------------------------
        self.beta - dictionary of dataframes
        """
        if type(folders) is not dict:
            raise ValueError('`folders` must be listlike')

        self.beta = {}  # empty hash

        # get betas for each genotype comparison:
        for genotype in folders.keys():
            path = main_path + folders[genotype] + fc_fname
            self.beta[genotype] = pd.read_csv(path, sep=sep)
            # beta dataframes from sleuth MUST BE SORTED By ID!!!!
            # if you don't do this, nothing works!
            self.beta[genotype].sort_values(self.gene, inplace=True)
            self.beta[genotype].reset_index(drop=True, inplace=True)

    def set_qval(self, q=0.1):
        """A function to set the qvalue parameter."""
        if type(q) is not float:
            raise ValueError('`q` must be a float!')

        if q == 0 or q == 1:
            raise ValueError('`q` must be between 0, 1 noninclusive')

        self.q = q

    def filter_data(self, count_min, count_quantile):
        """
        A function to filter data based on a set of conditions.

        Params:
        count_min - int or float
        count_quantile - float

        outputs:
        filtered_tpm
        filtered_beta
        """
        names = ['']
        for key, df in self.tpm.items():
            # drop anything in the bottom 10%
            quantile = (df.tpm.quantile(count_quantile))
            min_tpm = (df.tpm < quantile)
            # drop anything that has 0 tpm
            zero = (df.tpm <= count_min)
            # find those ids:
            series = df[min_tpm | zero][self.gene]
            names = names + pd.Series.tolist(series)
        # find the set of names (no repetitions) that don't pass the filter
        names = list(set(names))
        names = names[1:]

        na = []
        for key, genotype in enumerate(self.beta):

            # replace beta dfs with filtered values also
            ind = self.beta[genotype][self.gene].isin(names)
            # don't know if this works
            self.beta[genotype] = self.beta[genotype][~ind]

            # b value filtering:
            subset = ['b', 'qval']
            not_dropped = self.beta[genotype].dropna(axis=0,
                                                     subset=subset)[self.gene]
            ind = ~self.beta[genotype][self.gene].isin(not_dropped)
            na_here = self.beta[genotype][ind][self.gene].values

            na = na + pd.Series.tolist(na_here)

        na = list(set(na))
        print('Number of na genes: {0}'.format(len(na)))

        # filter everything that has beta == nan
        self.beta_filtered = {}
        self.tpm_filtered = {}
        for key, name in enumerate(self.tpm):
            # replace tpm dfs with the filtered values:
            filter1 = (~self.tpm[name][self.gene].isin(names))
            temp = self.tpm[name][filter1]
            filter2 = (~temp[self.gene].isin(na))
            temp2 = temp[filter2]
            self.tpm_filtered[name] = temp2.copy()
            self.tpm_filtered[name].reset_index(drop=True, inplace=True)

        for key, genotype in enumerate(self.beta):
            # replace beta dfs with filtered values also
            filter1 = (~self.beta[genotype][self.gene].isin(names))
            temp = self.beta[genotype][filter1]
            filter2 = (~temp[self.gene].isin(na))
            temp2 = temp[filter2]
            self.beta_filtered[genotype] = temp2.copy()
            self.beta_filtered[genotype].reset_index(drop=True, inplace=True)

    def enrichment_analysis(self, x, dictionary='tissue', analysis=''):
        """
        A neatly packaged enrichment analysis tool.

        Params:
        ------------
        dictionary - str or pandas df, at the moment, can only be 'tissue'
        if a pandas dataframe, must be in the same format as the dataframe
        used for TEA.

        Outputs:
        ------------
        enrichment - a hash of the results for each and every dataframe in
        the betas hash
        """
        if dictionary == 'tissue':
            dictionary = tea.fetch_dictionary()
            analysis = 'tissue'

        self.enrichment = {}

        if type(dictionary) is str:
            dictionary = tea.fetch_dictionary()
            self.enrichment[dictionary] = {}
        else:
            if len(analysis) == 0:
                raise ValueError('`analysis` can\'t be empty!')
            else:
                self.enrichment[analysis] = {}

        for key, df in self.beta.items():
            sig_genes = df[df[self.qval] < self.q][x]
            df_enr, unused = tea.enrichment_analysis(sig_genes, dictionary,
                                                     show=False, save=False)
            if len(df_enr) > 0:
                self.enrichment[analysis][key] = df_enr


###############################################################################
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
###############################################################################

class brenner(object):
    """
    An object to perform and hold fast analyses for morgan single-mutants.

    Attributes:
    ------------------
    rho
    hyper_plus
    overlap_plus
    expected_plus
    hyper_minus
    overlap_minus
    expected_minus
    """

    def __init__(self, name, morgan, alpha=10**-4):
        """
        Initialize function.

        Automatically performs spearman analysis and
        hypergeometric interaction tests.
        """
        self.name = name
        self.spearman_analysis(morgan, alpha=alpha)
        self.probabilistic_interaction_analysis(morgan)

    def spearman_analysis(self, morgan, alpha=10**-4):
        """
        A function to perform spearmanr analyses on all single mutants.

        Params:
        alpha - significance value for spearmanr correlation

        Outputs:
        res_dict - a hash containing the results of the analysis.
        """
        s = len(morgan.single_mutants)

        rho_matrix = np.empty(shape=(s, s))

        l = 0
        for i in morgan.single_mutants:
            m = 0
            for j in morgan.single_mutants:
                x = morgan.beta[i]
                y = morgan.beta[j]
                # find overlap in stat.sig.genes between both lists:
                indx = (x[morgan.qval] < 0.1)
                indy = (y[morgan.qval] < 0.1)
                ovx = x[indx]
                ovy = y[indy & y[morgan.gene].isin(ovx[morgan.gene])]
                ovx = x[x[morgan.gene].isin(ovy[morgan.gene])]

                if len(ovx) < 20:
                    continue

                # spearman analysis
                rho = stats.spearmanr(ovx[morgan.change], ovy.b)
                if rho[1] < alpha:
                    # return the spearman corr. weighted by fraction overlap
                    a_and_b = len(ovx)
                    a_or_b = len(x[indx]) + len(y[indy]) - a_and_b
                    rho_matrix[l, m] = rho[0]*a_and_b/a_or_b
                else:
                    rho_matrix[l, m] = 0
                m += 1
            l += 1
        self.rho = pd.DataFrame(data=rho_matrix, columns=morgan.single_mutants)
        self.rho['corr_with'] = morgan.single_mutants

    def find(self, x, y, morgan):
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
        return x[morgan.gene].isin(y)

    def polarize(self, test_df, ref_df, ids, sign, morgan):
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
         change - name of column containing fold change, regression value or a
             related metric for change relative to the WT. Absolute value
             doesn't matter, but it must not be absolute change!

        Output:
        g_overlap - numpy.ndarray of overlapping gene names
        """
        test = test_df[self.find(test_df, ids, morgan)].copy()
        ref = ref_df[self.find(ref_df, ids, morgan)]

        if ~(test[morgan.gene].values == ref[morgan.gene].values).all():
            test.sort_values(morgan.gene, inplace=True)
            test.reset_index(inplace=True)
            ref.sort_values(morgan.gene, inplace=True)
            ref.reset_index(inplace=True)

        test['polarity'] = test[morgan.change]*ref[morgan.change]

        if sign == '+':
            g_overlap = test[test.polarity > 0][morgan.gene].values
        elif sign == '-':
            g_overlap = test[test.polarity < 0][morgan.gene].values
        elif sign == 'p':
            g_overlap = 0
            print('Unfinished business!')
        else:
            raise ValueError('sign must be one of +, -, or p')
        return g_overlap

    def overlap(self, x, y, morgan):
        """
        Given two dataframes, returns the overlap between sig. genes.

        Params:
        x - pandas df
        y - pandas df

        output:
        overlap - a numpy ndarray
        """
        x_sig = x[x[morgan.qval] < morgan.q]
        y_sig = y[y[morgan.qval] < morgan.q]

        # overlap between x, y
        ind = x_sig[morgan.gene].isin(y_sig[morgan.gene].values)
        overlap = x_sig[ind][morgan.gene].values
        return overlap

    def overlap_prob(self, key1, key2, sign, morgan):
        """
        A function to calculate the hypergeom probability of overlap.

        Note: this function can't distinguish the direction of the network, so
        a ---> b and b ----> a are the same for the purposes of this function.

        Params:
        key1 - reference dataframe (must be df with largest # of sig. genes)
        key2 - test dataframe
        sign - "+", "-" or "p". + means activating relationship, "-" means
                inhibitory relationship

        Output:
        prob_overlap - probability of overlap as measured by hypergeometric
        overlap_fraction - fraction of overlap genes (A and B)/(A U B)
        expected_frac - expected fraction of overlap for these gene sets
        polarized_ids - gene ids of the overlapping set
        """
        ref_df = morgan.beta[key1]
        test_df = morgan.beta[key2]

        # find overlap:
        ovrlp_ids = self.overlap(ref_df, test_df, morgan)

        # sig genes:
        ref_sig_df = ref_df[ref_df[morgan.qval] < morgan.q].copy()
        test_sig_df = test_df[test_df[morgan.qval] < morgan.q].copy()

        # for overlapping ids, check what number of them satisfy the condition
        # specified by 'sign'
        # call them polarized bc they have a sign
        polarized_ids = self.polarize(test_df, ref_df, ovrlp_ids, sign,
                                      morgan)

        # turn this into a scalar:
        n_polar = len(polarized_ids)

        # genes that changed significantly in either comparison (non-redundant)
        temp = pd.Series.tolist(ref_sig_df[morgan.gene].values)
        g_sig_ref = list(set(temp))

        temp = pd.Series.tolist(test_sig_df[morgan.gene].values)
        g_sig_test = list(set(temp))

        # convert those lists to scalars:
        n_sig_ref = len(g_sig_ref)
        n_sig_test = len(g_sig_test)

        # to calculate the overlap fraction, we need to know how
        # many genes in A U B (A or B)
        total_genes_changed = len(list(set(g_sig_ref+g_sig_test)))

        # total genes measured in this experiment pair:
        genes_ref = pd.Series.tolist(ref_df[morgan.gene].values)
        genes_test = pd.Series.tolist(test_df[morgan.gene].values)
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

    def a_interacts_b(self, key1, key2, morgan, sign='+'):
        """
        A function to test whether a interacts with b in some way.

        a --> b or b --> a are the same thing for this function.

        Params:
        key1 - pandas dataframe
        key2 - pandas dataframe
        sign - '+', '-', 'p'. + tests activation, - tests inhibition, p not yet
                implemented.

        Ouput:
        overlap_p
        overlap_f
        expected
        ids_overlap
        """
        a_sig_genes = (morgan.beta[key1][morgan.qval] < morgan.q)
        b_sig_genes = (morgan.beta[key1][morgan.qval] < morgan.q)

        # name proxies for ease of reading
        a = morgan.beta[key1][a_sig_genes]
        b = morgan.beta[key2][b_sig_genes]

        # choose the reference set
        if len(a) > len(b):
            results = self.overlap_prob(key1, key2, sign, morgan)
        else:
            results = self.overlap_prob(key2, key1, sign, morgan)
        # unpack results
        overlap_p, overlap_f, expected_frac, ids_overlap = results
        # return
        return overlap_p, overlap_f, expected_frac, ids_overlap

    def probabilistic_interaction_analysis(self, morgan):
        """A function that predicts interactions based on Bayesian modeling."""
        s = len(morgan.single_mutants)

        prob_plus = np.zeros(shape=(s, s))
        overlap_plus = np.zeros(shape=(s, s))
        expected_plus = np.zeros(shape=(s, s))

        prob_minus = np.zeros(shape=(s, s))
        overlap_minus = np.zeros(shape=(s, s))
        expected_minus = np.zeros(shape=(s, s))

        self.overlap_ids_plus = {}
        self.overlap_ids_minus = {}

        l = 0
        for i in morgan.single_mutants:
            m = 0
            for j in morgan.single_mutants:
                # store the results from genpy.a_interacts_b in an array
                # called results, but remember it has 4 elements:
                # overlap prob, overlap frac, expected frac, ids overlapped
                results = self.a_interacts_b(i, j, morgan, sign='+')
                results2 = self.a_interacts_b(i, j, morgan, sign='-')

                # artificially set i,i entries for overlap fraction to zero,
                # this allows better discrimination of interactions for
                # heatmaps
                if i == j:
                    overlap_plus[l, m] = 0
                    overlap_minus[l, m] = 0
                else:
                    overlap_plus[l, m] = results[1]
                    overlap_minus[l, m] = results2[1]

                prob_plus[l, m] = results[0]
                prob_minus[l, m] = results2[0]
                expected_plus[l, m] = results[2]
                expected_minus[l, m] = results2[2]
                if i != j:
                    if len(results[3]):
                        self.overlap_ids_plus[(i, j)] = results[3]
                    if len(results2[3]):
                        self.overlap_ids_minus[(i, j)] = results2[3]
                m += 1
            l += 1

        self.hyper_plus = pd.DataFrame(prob_plus,
                                       columns=morgan.single_mutants)
        self.hyper_plus['corr_with'] = morgan.single_mutants

        self.overlap_plus = pd.DataFrame(overlap_plus,
                                         columns=morgan.single_mutants)
        self.overlap_plus['corr_with'] = morgan.single_mutants

        self.expected_plus = pd.DataFrame(expected_plus,
                                          columns=morgan.single_mutants)
        self.expected_plus['corr_with'] = morgan.single_mutants

        self.hyper_minus = pd.DataFrame(prob_minus,
                                        columns=morgan.single_mutants)
        self.hyper_minus['corr_with'] = morgan.single_mutants

        self.overlap_minus = pd.DataFrame(overlap_minus,
                                          columns=morgan.single_mutants)
        self.overlap_minus['corr_with'] = morgan.single_mutants

        self.expected_minus = pd.DataFrame(expected_minus,
                                           columns=morgan.single_mutants)
        self.expected_minus['corr_with'] = morgan.single_mutants


###############################################################################
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
###############################################################################

class mcclintock(object):
    """
    An object that performs bayesian robust regression on a morgan object.

    For single mutant analysis.

    Attributes:
    ------------------
    name
    robust_slope
    primary_weights
    secondary_slope
    secondary weights
    """

    def __init__(self, name, morgan):
        """
        Initialize function.

        Performs bayesian primary and secondary regression.
        """
        self.name = name
        self.robust_regression_primary(morgan)
        self.robust_regression_secondary(morgan)

    def mcmc_robust(self, data):
        """Bayesian Regression Using PyMC3."""
        # with pm.Model() as model_robust:
        with pm.Model():
            family = pm.glm.families.StudentT()
            pm.glm.glm('y ~ x', data, family=family)
            start = pm.find_MAP()
            step = pm.NUTS(scaling=start)
            trace_robust = pm.sample(2000, step, progressbar=True)

        return trace_robust

    def robust_regression_primary(self, morgan, alpha=10**-4):
        """
        A function to perform robust spearmanr analyses on all single mutants.

        Params:
        alpha - significance value for spearmanr correlation

        Outputs:
        res_dict - a hash containing the results of the analysis.
        """
        def find_rank(morgan, df):
            """A function to find the rank values of a variable."""
            d = df.copy()
            d.sort_values(morgan.change, inplace=True)
            rank = np.linspace(0, len(d)-1, len(d))
            d['r'] = rank
            d.sort_values(morgan.gene, inplace=True)
            return d

        def find_inliers(morgan, ovx, ovy, trace):
            """A function to find inliers from the Bayesian regression."""
            # find the mean and std of the distribution along the line
            mean = np.mean(ovy.r - trace.Intercept.mean() -
                           ovx.r*trace.x.mean())
            std = np.std(ovy.r - trace.Intercept.mean() -
                         ovx.r*trace.x.mean())
            # find the total distribution:
            intercept = trace.Intercept.mean()
            slope = trace.x.mean()
            distribution = ovy.r - intercept - ovx.r*slope

            # call the inliers and outliers.
            # fairly aggressive -- < 1std is inlier, > is outlier
            inliers = (np.abs(distribution - mean)/std < 1)

            # get a list of the gene candidates (genes close to line)
            candidates = ovy[ovy.r.isin(ovy.r[inliers])][morgan.gene]
            return candidates

        def perform_mcmc(morgan, ovx, ovy, mut_a, mut_b):
            """
            A function to perform the robust spearmanr regress.

            Written mainly to avoid running into RAM issues. Not
            entirely meant for public use.
            ovx, ovy -- dataframes to be correlated
            mut_a, mut_b -- genotypes of ovx and ovy
            """
            # rank order:
            ovx = find_rank(morgan, ovx)
            ovy = find_rank(morgan, ovy)

            # place in a dict
            data = dict(x=ovx.r, y=ovy.r)
            # run PyMC3 with student T distribution
            # to minimize impact of outliers
            print('\nstarting comparison of {0}, {1}'.format(i, j))
            trace_robust = self.mcmc_robust(data)

            # find the mean and std of the distribution along the line
            candidates = find_inliers(morgan, ovx, ovy, trace_robust)
            # also get a list of the outliers
            outliers = ovy[~ovy[morgan.gene].isin(candidates)]

            # place in a list:
            temp1 = candidates
            temp2 = outliers
            self.correlated_genes[(mut_a, mut_b)] = {'corr': temp1,
                                                     'outliers': temp2
                                                     }
            # fill in the slope matrix
            tmean = trace_robust.x.mean()
            tstd = trace_robust.x.std()

            # mean value of the slope:
            t = trace_robust.x.mean()/np.abs(trace_robust.x.mean())

            return tmean, tstd, outliers, t

        # begin the main function
        self.correlated_genes = {}
        s = len(morgan.single_mutants)

        slope_matrix = np.zeros(shape=(s, s))
        error_matrix = np.zeros(shape=(s, s))

        for l in range(len(morgan.single_mutants)):
            for m in range(l, len(morgan.single_mutants)):
                if l == m:
                    continue
                # only move forward if l != m
                i = morgan.single_mutants[l]
                j = morgan.single_mutants[m]
                x = morgan.beta[i]
                y = morgan.beta[j]
                # find overlap in stat.sig.genes between both lists:
                indx = (x[morgan.qval] < morgan.q)
                indy = (y[morgan.qval] < morgan.q)

                ovx = x[indx]
                ovy = y[indy & y[morgan.gene].isin(ovx[morgan.gene])]
                ovx = x[x[morgan.gene].isin(ovy[morgan.gene])]

                # don't do anything if overlap is less than 20 items
                if len(ovx) < 20:
                    continue

                results = perform_mcmc(morgan, ovx, ovy, i, j)
                tmean, tstd, outliers, t = results

                # calculate overlap
                # overlap here is considered to be only the inliers to the
                # regression
                a_and_b = len(ovx) - len(outliers)
                a_or_b = len(x[indx]) + len(y[indy]) - a_and_b
                slope_matrix[l, m] = tmean*a_and_b/a_or_b
                error_matrix[l, m] = tstd*a_and_b/a_or_b

        # place in a neat dataframe so the user can see this.
        self.robust_slope = pd.DataFrame(data=slope_matrix,
                                         columns=morgan.single_mutants)
        self.robust_slope['corr_with'] = morgan.single_mutants

    def robust_regression_secondary(self, morgan, frac=0.1):
        """
        A function to find secondary interactions.

        Only use if you have first called robust_regression_primary.
        """
        def find_rank(morgan, df):
            """A function to find the rank values of a variable."""
            d = df.copy()
            d.sort_values(morgan.change, inplace=True)
            rank = np.linspace(0, len(d)-1, len(d))
            d['r'] = rank
            d.sort_values(morgan.gene, inplace=True)
            return d

        s = len(morgan.single_mutants)
        matrix = np.zeros(shape=(s, s))
        for l in range(0, s):
            for m in range(l, s):
                print(morgan.single_mutants[l], morgan.single_mutants[m])
                if l == m:
                    continue
                letter1 = morgan.single_mutants[l]
                letter2 = morgan.single_mutants[m]

                if (letter1, letter2) not in self.correlated_genes.keys():
                    continue

                # find significant overlap between the two genotypes:
                significance1 = morgan.beta_filtered[letter1].qval < morgan.q
                x = morgan.beta_filtered[letter1][significance1]

                significance2 = morgan.beta_filtered[letter2].qval < morgan.q
                inx = morgan.beta_filtered[letter2][morgan.gene].\
                    isin(x[morgan.gene])

                y = morgan.beta_filtered[letter2][(significance2) &
                                                  (inx)].copy()
                # store the sizes of x and y:
                sizex = len(x)
                sizey = len(y)

                x = x[x[morgan.gene].isin(y[morgan.gene])].copy()

                # rank order values
                x = find_rank(morgan, x)
                y = find_rank(morgan, y)

                # get the outliers from the previous run:
                genes = self.correlated_genes[(letter1,
                                              letter2)]['outliers']
                X = x[x[morgan.gene].isin(genes[morgan.gene])]
                Y = y[y[morgan.gene].isin(genes[morgan.gene])]

                # only proceed if there are enough genes to test:
                if len(X) > 20:
                    # do the mcmc
                    data = dict(x=X.r, y=Y.r)
                    trace = self.mcmc_robust(data)
                    a_and_b = len(genes)
                    a_or_b = sizex + sizey - a_and_b
                    matrix[l, m] = trace.x.mean()*a_and_b/a_or_b
                    del(trace)
        # return a non-tidy dataframe
        self.secondary_slope = pd.DataFrame(data=matrix,
                                            columns=morgan.single_mutants)


###############################################################################
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
###############################################################################

class sturtevant(object):
    """
    An object to perform epistasis analysis of double mutants.

    Attributes
    ------------------
    epistasis
    """

    def __init__(self, name):
        """Initialize function."""
        self.name = name

    def robust_regress(self, data):
        """A robust regression using a StudentT instead of a Gaussian model."""
        with pm.Model():
            family = pm.glm.families.StudentT()
            pm.glm.glm('y ~ x', data, family=family)
            start = pm.find_MAP()
            step = pm.NUTS(scaling=start)
            trace_robust = pm.sample(2000, step, progressbar=True)
            return trace_robust

    def epistasis_analysis(self, morgan):
        """
        A function to perform correlational epistatic analysis of mutants.

        Parameters:
        -----------
        morgan
        """
        def find_rank(morgan, df):
            """A function to find the rank values of a variable."""
            d = df.copy()
            d.sort_values('b', inplace=True)
            rank = np.linspace(0, len(d)-1, len(d))
            d['r'] = rank
            d.sort_values(morgan.gene, inplace=True)
            return d

        def find_inliers(morgan, ovx, ovy, trace):
            """A function to find inliers from the Bayesian regression."""
            # find the mean and std of the distribution along the line
            mean = np.mean(ovy.r - trace.Intercept.mean() -
                           ovx.r*trace.x.mean())
            std = np.std(ovy.r - trace.Intercept.mean() -
                         ovx.r*trace.x.mean())
            # find the total distribution:
            intercept = trace.Intercept.mean()
            slope = trace.x.mean()
            distribution = ovy.r - intercept - ovx.r*slope

            # call the inliers and outliers.
            # fairly aggressive -- < 1std is inlier, > is outlier
            inliers = (np.abs(distribution - mean)/std < 1)

            # get a list of the gene candidates (genes close to line)
            candidates = ovy[ovy.r.isin(ovy.r[inliers])][morgan.gene]
            return candidates

        # lambda index function:
        def lind(x):
            """To avoid writing ind = blah < q all the time."""
            return (x.qval < 0.1)

        double_mat = np.zeros(shape=(2, 5))  # TODO Fix this!
        weights = np.zeros(shape=(2, 5))
        l = 0
        cols = []
        size = []
        for key in morgan.double_muts:
            m = 0
            cols += [key]
            for j in morgan.single_mutants:
                x = morgan.beta_filtered[key]
                y = morgan.beta_filtered[j]

                ovx = x[lind(x)]
                ovy = y[lind(y) & y[morgan.gene].isin(ovx[morgan.gene])].copy()
                ovx = x[lind(x) & x[morgan.gene].isin(ovy[morgan.gene])].copy()

                ovx = find_rank(morgan, ovx)
                ovy = find_rank(morgan, ovy)

                # size of the overlap between double and single
                size += [len(ovx)]

                # go to next comparison if size is less than 20
                if len(ovx) < 20:
                    continue

                # mean_slope = stats.spearmanr(ovx.b, ovy.b)[0]
                # Set up MCMC parameters
                data = dict(x=ovx.r, y=ovy.r)
                trace = self.robust_regress(data)
                candidates = find_inliers(morgan, ovx, ovy, trace)
                a_or_b = len(candidates)
                a_and_b = len(x[lind(x)]) + len(y[lind(y)]) - a_or_b
                weight = a_or_b/a_and_b
                mean_slope = trace.x.mean()*weight
                print('\n{0:.2g}'.format(mean_slope))
                double_mat[l, m] = mean_slope
                weights[l, m] = weight
                m += 1
            l += 1

        # place results in dataframe:
        double_corr = pd.DataFrame(double_mat.transpose(), columns=cols)
        double_corr = pd.melt(double_corr, var_name='double_mutant',
                              value_name='correlation')
        double_corr['corr_with'] = morgan.single_mutants*2

        w = pd.DataFrame(weights.transpose(), columns=cols)
        w = pd.melt(w, var_name='double_mutant',
                    value_name='weights')
        # # correct the correlations by the size of the overlap
        # double_corr['s'] = size
        #
        # def weighted_corr(x):
        #     m = len(morgan.beta_filtered[x.double_mutant])
        #     return x.correlation*x.s/m
        #
        # # weight the correlation by the size of the group that was compared,
        # # then normalize to one by dividing by the max
        # double_corr['weighted_corr'] = double_corr.apply(weighted_corr,
        #                                                  axis=1)

        self.epistasis = double_corr
        self.epistasis['weights'] = w['weights']

###############################################################################
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
###############################################################################
