"""
Morgan.

authors: David Angeles-Albores and Paul W. Sternberg
contact: dangeles at caltech edu
"""
import pandas as pd
import warnings as wng
import numpy as np
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
        self.beta = None

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
            return

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
        batch - batch each project belonged to

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
        columns = ['project_name', 'genotype', 'batch']
        if (self.genmap.columns != columns).all():
            raise ValueError('genmap is not in the right format!')

        self.genmap.genotype = self.genmap.genotype.apply(str)
        # make sure everything is always in lowercase
        self.genmap.genotype = self.genmap.genotype.apply(str.lower)

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
            self.beta[genotype].sort_values(self.gene, inplace=True)
            self.beta[genotype].reset_index(drop=True, inplace=True)

    def add_beta(self, fname, key, **kwargs):
        """A function to add a file to the beta dictionary."""
        if self.beta is None:
            self.beta = {}
        self.beta[key] = pd.read_csv(fname, **kwargs)

    def set_qval(self, q=0.1):
        """A function to set the qvalue parameter."""
        if type(q) is not float:
            raise ValueError('`q` must be a float!')
        if q == 0 or q == 1:
            raise ValueError('`q` must be between 0, 1 noninclusive')
        self.q = q

    def filter_data(self):
        """
        A function to filter out NaNs in the beta dataframes.

        Params:
        count_min - int or float
        count_quantile - float

        outputs:
        filtered_tpm
        filtered_beta
        """
        for genotype, df in self.beta.items():
            df.dropna(subset=['b'], inplace=True)

###############################################################################
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
###############################################################################
# some common functions


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


def robust_regress(data, progress=False):
    """A robust regression using a StudentT instead of a Gaussian model."""
    with pm.Model():
        family = pm.glm.families.StudentT()
        pm.glm.glm('y ~ x', data, family=family)
        start = pm.find_MAP()
        step = pm.NUTS(scaling=start)
        trace_robust = pm.sample(2000, step, progressbar=progress)
    return trace_robust


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

    def __init__(self, name, morgan, progress):
        """
        Initialize function.

        Performs bayesian primary and secondary regression.
        """
        self.name = name
        self.progress = progress
        self.robust_regression_primary(morgan, progress)
        self.robust_regression_secondary(morgan, progress)

    def mcmc_robust(self, data, progress=True):
        """Bayesian Regression Using PyMC3."""
        # with pm.Model() as model_robust:
        with pm.Model():
            family = pm.glm.families.StudentT()
            pm.glm.glm('y ~ x', data, family=family)
            start = pm.find_MAP()
            step = pm.NUTS(scaling=start)
            trace_robust = pm.sample(2000, step, progressbar=progress)

        return trace_robust

    def robust_regression_primary(self, morgan, alpha=10**-4, progress=True):
        """
        A function to perform robust spearmanr analyses on all single mutants.

        Params:
        alpha - float, significance value for spearmanr correlation
        progress - Boolean, show progressbar for mcmc

        Outputs:
        res_dict - a hash containing the results of the analysis.
        """
        def perform_mcmc(morgan, ovx, ovy, mut_a, mut_b, progress=True):
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
            trace_robust = robust_regress(data, progress)

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
        #####################################################################
        #####################################################################
        # begin robust regression ###########################################
        #####################################################################
        #####################################################################
        self.correlated_genes = {}
        s = len(morgan.single_mutants)

        slope_matrix = np.zeros(shape=(s, s))
        weights = np.zeros(shape=(s, s))
        error_matrix = np.zeros(shape=(s, s))
        if len(morgan.single_mutants) == 0:
            raise ValueError('morgan single_mutants is empty!')
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

                results = perform_mcmc(morgan, ovx, ovy, i, j, progress)
                tmean, tstd, outliers, t = results

                # calculate overlap
                # overlap here is considered to be only the inliers to the
                # regression
                a_and_b = len(ovx) - len(outliers)
                a_or_b = len(x[indx]) + len(y[indy]) - a_and_b
                # slope_matrix[l, m] = tmean*a_and_b/a_or_b
                # error_matrix[l, m] = tstd*a_and_b/a_or_b
                weight = a_and_b/a_or_b
                slope_matrix[l, m] = tmean*weight
                error_matrix[l, m] = tstd*weight
                weights[l, m] = weight

        # place in a neat dataframe so the user can see this.
        self.robust_slope = pd.DataFrame(data=slope_matrix,
                                         columns=morgan.single_mutants)
        self.robust_slope['corr_with'] = morgan.single_mutants

        self.errors_primary = pd.DataFrame(data=error_matrix,
                                           columns=morgan.single_mutants)
        self.errors_primary['corr_with'] = morgan.single_mutants

        w = pd.DataFrame(weights, columns=morgan.single_mutants)
        self.weights_primary = w
        self.weights_primary['corr_with'] = morgan.single_mutants

    def robust_regression_secondary(self, morgan, frac=0.1):
        """
        A function to find secondary interactions.

        Only use if you have first called robust_regression_primary.
        """
        s = len(morgan.single_mutants)
        matrix = np.zeros(shape=(s, s))
        error_matrix = np.zeros(shape=(s, s))
        weights = np.zeros(shape=(s, s))
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
                significance1 = morgan.beta[letter1].qval < morgan.q
                x = morgan.beta[letter1][significance1]

                significance2 = morgan.beta[letter2].qval < morgan.q
                inx = morgan.beta[letter2][morgan.gene].\
                    isin(x[morgan.gene])

                y = morgan.beta[letter2][(significance2) &
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
                    trace = robust_regress(data)
                    a_and_b = len(genes)
                    a_or_b = sizex + sizey - a_and_b
                    # total = len(morgan.beta[letter1][significance1])
                    weight = a_and_b/a_or_b
                    matrix[l, m] = trace.x.mean()*weight
                    tstd = trace.x.std()
                    error_matrix[l, m] = tstd*weight
                    weights[l, m] = weight
                    del(trace)
        # return a non-tidy dataframe
        self.secondary_slope = pd.DataFrame(data=matrix,
                                            columns=morgan.single_mutants)

        self.errors_secondary = pd.DataFrame(data=error_matrix,
                                             columns=morgan.single_mutants)
        self.errors_secondary['corr_with'] = morgan.single_mutants

        w = pd.DataFrame(weights, columns=morgan.single_mutants)
        self.weights_secondary = w
        self.weights_secondary['corr_with'] = morgan.single_mutants
