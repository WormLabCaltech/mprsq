"""A script that contains all functions to do RNA-seq epistasis analysis."""
# important stuff:
import pandas as pd
import numpy as np

# Graphics
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.odr as odr

# labeller:
import gvars
from scipy.stats import gaussian_kde
from matplotlib import rc

rc('text', usetex=True)
rc('text.latex', preamble=r'\usepackage{cmbright}')
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})

# JB's favorite Seaborn settings for notebooks
rc = {'lines.linewidth': 2,
      'axes.labelsize': 18,
      'axes.titlesize': 18,
      'axes.facecolor': 'DFDFE5'}
sns.set_context('paper', rc=rc)
sns.set_style("dark")

mpl.rcParams['xtick.labelsize'] = 16
mpl.rcParams['ytick.labelsize'] = 16
mpl.rcParams['legend.fontsize'] = 16

genvar = gvars.genvars()
epi = gvars.epistasis()


def label(code1, code2):
    """A function to make epistasis labels given two code-letters."""
    return '{0} $>$ {1}'.format(genvar.fancy_mapping[code1],
                                genvar.fancy_mapping[code2])


def find_overlap(genotypes, df, q=0.1, col='code'):
    """Given a list of genotypes, df and a q-value, find DEG common to all."""
    # find only DE genes:
    sig = df[(df[col].isin(genotypes)) & (df.qval < q)]
    grouped = sig.groupby('target_id')
    genes = []
    for target, group in grouped:
        # make sure the group contains all desired genotypes
        all_in = (len(group.code.unique()) == len(genotypes))
        if all_in:
            genes += [target]
    return genes


def find_additive(single_muts, double_mut, df, q=0.1):
    """
    Given 3 genotypes, find shared DE genes and return sliced dataframes.

    Params:
    single_muts - a list containing exactly two elements
    double_muts - a code for a double mutant
    df - a tidy dataframe. must contain columns 'target_id' and 'code'

    Output:
    x, y, xy - Dataframes containing DE genes and all relevant information
    (Betas, q-values, etc...). First dataframe corresponds to single_muts[0],
    second dataframe corresponds to genotype single_muts[1] and third is the
    double mutant information.
    """
    if type(single_muts) is not list:
        raise ValueError('single_muts must be of type list')
    if type(double_mut) is not str:
        raise ValueError('double_mut must be of type str')

    # find the overlapping gene list
    genes = find_overlap(single_muts + [double_mut], df)

    # extract the dataframes
    x = df[(df.target_id.isin(genes)) &
           (df.code == single_muts[0])]
    y = df[(df.target_id.isin(genes)) &
           (df.code == single_muts[1])]
    xy = df[(df.target_id.isin(genes)) &
            (df.code == double_mut)]

    # return the dataframes
    return x, y, xy


def f(B, x):
    """A linear function for the ODR."""
    return B*(x)


def perform_odr(add, dev, wadd, wdev):
    """A wrapper to calculate an ODR regression."""
    linear = odr.Model(f)
    # mydata = odr.Data(add, dev, wd=1./wadd, we=1./wdev)
    mydata = odr.RealData(add, dev, sx=wadd, sy=wdev)
    myodr = odr.ODR(mydata, linear, beta0=[0])
    myoutput = myodr.run()
    return myoutput


def ODR(singles, double, epistasis):
    """Find the ODR in epistasis plot between single muts and a double mut."""
    # errors:
    if len(singles) != 2:
        raise ValueError('`singles` must be a list with two dataframes!')
    if type(double) is not pd.DataFrame:
        raise ValueError('`double` must be a dataframe!')
    try:
        epistasis = epistasis.lower()
    except:
        raise ValueError('epistasis must be a string!')
    if epistasis not in ['actual', 'xy=x', 'xy=y']:
        raise ValueError('epistasis must be one of `actual`, `xy=x`, `xy=y`')

    # define the X-coordinate as the additive model of interaction
    X = singles[0].b + singles[1].b

    # fit an ODR model
    wadd = np.sqrt(singles[1].se_b**2 + singles[0].se_b**2)

    if epistasis == 'actual':
        # calculate deviation standard error:
        wdev = double.se_b**2
        for i, df in enumerate(singles):
            wdev += df.se_b**2
        wdev = np.sqrt(wdev)
        # calculate:
        output = perform_odr(X, double.b - X, wadd=wadd, wdev=wdev)
    if epistasis == 'xy=x':
        # if XY = X, then XY - X - Y = -Y
        output = perform_odr(X, -singles[1].b, wadd=wadd, wdev=singles[1].se_b)
    if epistasis == 'xy=y':
        # if XY = Y, then XY - X - Y = -X
        output = perform_odr(X, -singles[0].b, wadd=wadd, wdev=singles[0].se_b)

    return output


def plot_epistasis_regression(X, slope, **kwargs):
    """Plot the ODR line."""
    # find the xmin and xmax:
    xmin = X.min()
    xmax = X.max()
    x = np.linspace(xmin - 0.1, xmax + 0.1, 1000)
    y0 = x*slope
    # plot the models
    plt.plot(x, y0, **kwargs)


def draw_bs_sample(n):
    """Draw a bootstrap sample from a 1D data set."""
    ind = np.arange(0, n)
    return np.random.choice(ind, size=n)


def bootstrap(bframe, sebframe, epistasis='actual', nsim=1000):
    """
    Perform non-parametric bootstrapping for an epistasis ODR.

    Given a list of three numpy vectors containing betas and a separate list of
    vectors containing their standard errors, fit a model according to the
    `epistasis` parameter indicated and bootstrap it. The vectors MUST
    be provided in the order [X, Y, XY], where X is the first genotype, Y is
    the second genotype and XY is the double mutant.

    Params:
    bframe - a list of numpy vectors containing the betas for each genotype
    sebframe - a list of numpy vectors containing the se_b for each genotype
    epistasis - kind of model to simulate. One of:
                'actual', 'suppress', 'xy=x+y', 'xy=x', 'xy=y','xy=x=y'.
    nsim - number of iterations to be performed. Must be >0

    Output:
    s, se_s
    """
    nsim = int(nsim)
    # unpack
    xb, yb, xyb = bframe
    xseb, yseb, xyseb = sebframe

    s = np.zeros(nsim)

    # draw bootstrap repetitions
    for i in range(nsim):
        # sample data, keeping tuples paired:
        ind = draw_bs_sample(len(xb))

        currx = xb[ind]
        curry = yb[ind]
        currxy = xyb[ind]

        currsex = xseb[ind]
        currsey = yseb[ind]
        currsexy = xyseb[ind]

        # different bootstraps to do:
        # for the actual data, do a non-parametric bootstrap
        wadd = np.sqrt(currsex**2 + currsey**2)
        if epistasis == 'actual':
            X = currx + curry
            Y = currxy - X
            wdev = wadd**2 + currsexy**2

        elif epistasis == 'xy=x':
            X = currx + curry
            Y = -curry
            wdev = currsey

        elif epistasis == 'xy=y':
            X = currx + curry
            Y = -currx
            wdev = currsex

        # for all others, do a parametric bootstrap
        # because we know what the slope should be,
        # but we need to generate a structure to test
        # against. Non-parametric bootstrapping will
        # yield perfect lines every time.
        elif epistasis == 'xy=x+y':
            X = currx + curry
            Y = np.random.normal(0, wadd, len(X))
            wdev = wadd

        elif epistasis == 'xy=x=y':
            # flip a coin:
            coin = np.random.randint(0, 1)

            # half the time use the X data
            # half the time use the Y
            if coin == 0:
                wadd = np.sqrt(2*currsex**2)
                wdev = currsex
                X = 2*currx + np.random.normal(0, wadd, len(curry))
                Y = -currx + np.random.normal(0, wdev, len(currx))

            else:
                wadd = np.sqrt(2)*currsey
                wdev = currsey
                X = 2*curry + np.random.normal(0, wadd, len(curry))
                Y = -curry + np.random.normal(0, wdev, len(curry))

        elif epistasis == 'suppress':
            # flip a coin:
            coin = np.random.randint(0, 2)

            # half the time use the X data
            # half the time use the Y
            if coin == 0:
                wadd = np.sqrt(2)*currsex
                wdev = currsey
                X = curry + np.random.normal(0, wadd, len(curry))
                Y = -curry + np.random.normal(0, wdev, len(curry))

            else:
                wadd = np.sqrt(2)*currsex
                wdev = currsex
                X = currx + np.random.normal(0, wadd, len(currx))
                Y = -currx + np.random.normal(0, wdev, len(currx))

        # do calcs and store in vectors
        output = perform_odr(X, Y, wadd=wadd, wdev=wdev)

        # extract the slope and standard error from the output
        # and store it
        s[i] = output.beta[0]
        # se_s[i] = output.sd_beta[0]

    return s


def bootstrap_regression(singles, double, df, epistasis='actual', nsim=100):
    """
    Perform a bootstrap regression for the desired epistatic model.

    Params:
    singles - a list of 2 genotypes that make up the double mutant
    double - a string containing the ID of the double mutant.
    df - a tidy dataframe. must have columns 'target_id', 'b', 'se_b', 'qval'
         'code', and 'genotype'
    epistasis - kind of model to simulate. One of:
                'actual', 'suppress', 'xy=x+y', 'xy=x', 'xy=y','xy=x=y'.
    nsim - number of simulations to perform

    Outputs:
    s - numpy vector containing all the ODR slope values from the bootstrap
    se_s - numpy vector containing all the ODR standard error of the slope
           values from the bootstrap
    """
    nsim = int(nsim)
    x, y, xy = find_additive(singles, double, df)

    xb = x.b.values
    yb = y.b.values
    xyb = xy.b.values

    xseb = x.se_b.values
    yseb = y.se_b.values
    xyseb = xy.se_b.values

    beta = bootstrap([xb, yb, xyb],
                     [xseb, yseb, xyseb],
                     epistasis=epistasis,
                     nsim=nsim)
    return beta


def epiplot(X, Y, Y_se, **kwargs):
    """Given two arrays, X and Y, plot the points."""
    plot_unbranched = kwargs.pop('plot_unbranched', False)
    beta = kwargs.pop('beta', np.nan)
    s0 = kwargs.pop('s0', 15)

    # Calculate the point density
    points = np.vstack([X, Y])
    z = gaussian_kde(points)(points)

    # plot:
    fig, ax = plt.subplots()
    if len(X) > 50:
        ax.scatter(X, Y, c=z, s=s0/Y_se,
                   edgecolor='', cmap='viridis', alpha=0.5)
    else:
        ax.scatter(X, Y, s=s0/np.sqrt(Y_se),
                   color='#33a02c', alpha=.9)

    if plot_unbranched:
        smoothX = np.linspace(X.min() - 0.5, X.max() + 0.5, 1000)
        plt.plot(smoothX, -1/2*smoothX, color='#1f78b4', ls='--',
                 label='Unbranched Pathway')
    if beta:
        plot_epistasis_regression(X, beta, ls='-', lw=2.3,
                                  color='#33a02c', label='fit')

    plt.xlabel(r'Predicted Additive Effect')
    plt.ylabel(r'Deviation from Additive Effect')

    plt.legend()
    return ax


def make_epiplot(singles, double, df, **kwargs):
    """
    Draw an epistasis plot of the data.

    Params:
    singles - a list of 2 genotypes that make up the double mutant
    double - a string containing the ID of the double mutant.

    Output:
    x - tidy dataframe containing the DE gene data for singles[0]
    y - tidy dataframe containing the DE gene data for singles[1]
    xy - tidy dataframe containing the DE gene data for the double mutant
    ax - axis containing the plot
    """
    x, y, xy = find_additive(singles, double, df)
    actual = ODR([x, y], xy, 'actual')

    # transform coordinates:
    X = x.b + y.b
    Y = xy.b - X
    Y_se = np.sqrt(x.se_b**2 + y.se_b**2 + xy.se_b**2)

    ax = epiplot(X, Y, Y_se, plot_unbranched=True, beta=actual.beta)

    return x, y, xy, ax


def calculate_all_bootstraps(x, y, xy, df, nsim=5000):
    """
    Given two double mutants and a double find the bootstrapped epistasis coef.

    Params:
    x
    y
    xy
    df
    nsim

    Output:
    epicoef, epierr
    """
    models = epi.models
    epicoef = {}
    for model in models:
        s = bootstrap_regression([x, y], xy, df,
                                 epistasis=model, nsim=nsim)
        epicoef[model] = s
    return epicoef


def plot_bootstraps(x, y, epicoef, **kwargs):
    """Make KDE plots of the bootstrapped epistasis coefficients."""
    # make dictionaries for plotting
    colors = {'actual': '#33a02c', 'xy=x': 'blue', 'xy=y': 'k',
              'xy=x=y': '#1f78b4', 'xy=x+y': '#ff7f00', 'suppress': '#e31a1c'
              }
    labels = {'actual': 'actual', 'xy=x': label(x, y),
              'xy=y': label(y, x), 'xy=x=y': 'Unbranched',
              'xy=x+y': 'Additive', 'suppress': 'Suppression'
              }

    # checks and balances
    if type(epicoef) is not dict:
        raise ValueError('epicoef must be a dictionary')

    epistasis_choice = ['actual', 'xy=x', 'xy=y', 'xy=x=y', 'xy=x+y',
                        'suppress']

    for epistasis in epistasis_choice:
        if epistasis.lower() not in epicoef.keys():
            warning = 'epicoef must contain keys for all epistasis models'
            raise ValueError(warning)

        if len(epicoef[epistasis.lower()]) < 10:
            warning = 'too few bootstraps. Please perform >100' + \
                      'bootstraps per test'
            raise ValueError(warning)

    fig, ax = plt.subplots()
    for model, s in epicoef.items():
        try:
            sns.kdeplot(data=s, label=labels[model.lower()],
                        color=colors[model.lower()], **kwargs)
        except:
            print('{0} did not have a label'.format(model))
            next

    # plot a horizontal line wherever the actual data mean is
    plt.gca().axvline(epicoef['actual'].mean(), color='#33a02c', ls='--', lw=3)

    plt.xlabel('Epistasis Coefficient')
    plt.ylabel('Cumulative Density Function')

    return ax


def permutation_test(s):
    """Perform a permutation test on the slope the genetic data."""
    epistasis = ['xy=x', 'xy=y', 'xy=x=y', 'xy=x+y',
                 'suppress']

    diff = {}
    for epi in epistasis:
        d = [s['actual'][i] - s[epi][i] for i in range(len(s[epi]))]
        diff[epi] = d

    return diff


def message(name, pval, alpha=0.01):
    """Write a message."""
    if pval < alpha:
        return '{0} can be rejected (pval <= {1:.2g})'.format(name, pval)
    else:
        return '{0} cannot be rejected (pval = {1:.2g})'.format(name, pval)


def calculate_pval(s, diff):
    """Given `s` and  `diff`, print out the p-values for each comparison."""
    for key, array in diff.items():
        # test =
        if s[key].mean() > s['actual'].mean():
            pval = len(array[array > 0])/len(array)
        else:
            pval = len(array[array < 0])/len(array)
        if pval == 0:
            p = 1/(len(array)/10)
        else:
            p = pval
        mess = message(key, p)
        print(mess)
