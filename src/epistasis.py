"""A script that contains all functions to do RNA-seq epistasis analysis."""
# important stuff:
# import os
import pandas as pd
import numpy as np
# from IPython.core.display import HTML

# Graphics
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.odr as odr

# labeller:
import gvars

# from IPython.core.display import HTML
# # bokeh
# import bokeh.charts
# import bokeh.charts.utils
# import bokeh.io
# import bokeh.models
# import bokeh.palettes
# import bokeh.plotting
# from bokeh.plotting import figure
# from bokeh.resources import CDN
# from bokeh.embed import file_html

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


def label(code1, code2):
    """A function to make epistasis labels given two code-letters."""
    return '{0} $>$ {1}'.format(genvar.fancy_mapping[code1],
                                genvar.fancy_mapping[code2])


def find_overlap(genotypes, df, q=0.1):
    """Given a n genotypes, df and a q-value, find genes that are DE in all."""
    # find only DE genes:
    sig = df[(df.code.isin(genotypes)) & (df.qval < q)]
    grouped = sig.groupby('target_id')
    genes = []
    for target, group in grouped:
        # # make sure all q-values are significant
        # q_sig = (group.qval < q).all()
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
    mydata = odr.Data(add, dev, wd=1./wadd, we=1./wdev)
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


def bootstrap(bframe, sebframe, epistasis='actual', nsim=1000):
    """
    Perform parametric bootstrapping for an epistasis ODR.

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
    se_s = np.zeros(nsim)
    for i in range(nsim):
        # test the construction.
        indices = np.random.randint(0, len(xb), len(xb))

        # extract the array stuff:
        currx = xb[indices]
        curry = yb[indices]
        currxy = xyb[indices]

        currsex = xseb[indices]
        currsey = yseb[indices]
        currsexy = xyseb[indices]

        # different bootstraps to do:
        if epistasis == 'actual':
            X = currx + curry
            Y = currxy - X
            wadd = np.sqrt(currsex**2 + currsey**2)
            wdev = wadd**2 + currsexy**2

        elif epistasis == 'XY=X':
            X = currx + curry
            Y = -curry
            wadd = np.sqrt(currsex**2 + currsey**2)
            wdev = currsey

        elif epistasis == 'XY=Y':
            X = currx + curry
            Y = -currx
            wadd = np.sqrt(currsex**2 + currsey**2)
            wdev = currsex

        elif epistasis == 'XY=X+Y':
            X = currx + curry
            Y = 0
            wadd = np.sqrt(currsex**2 + currsey**2)
            wdev = currsex

        elif epistasis == 'XY=X=Y':
            # flip a coin:
            coin = np.random.randint(0, 1)

            # half the time use the X data
            # half the time use the Y
            if coin == 0:
                X = 2*currx
                Y = -currx
                wadd = np.sqrt(2*currsex**2)
                wdev = currsex

            else:
                X = 2*curry
                Y = -curry
                wadd = np.sqrt(2*currsey**2)
                wdev = currsex

        elif epistasis == 'suppress':
            # flip a coin:
            coin = np.random.randint(0, 1)

            # half the time use the X data
            # half the time use the Y
            if coin == 0:
                X = curry
                Y = -curry
                wadd = currsey
                wdev = currsex

            else:
                X = currx
                Y = -currx
                wadd = currsex
                wdev = currsey

        # turn X, Y into random vars:
        X = np.random.normal(X, wadd)
        Y = np.random.normal(Y, wdev)

        # do calcs and store in vectors
        output = perform_odr(X, Y, wadd=wadd, wdev=wdev)

        # extract the slope and standard error from the output
        # and store it
        s[i] = output.beta[0]
        se_s[i] = output.sd_beta[0]

    return s, se_s


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

    beta, sebeta = bootstrap([xb, yb, xyb],
                             [xseb, yseb, xyseb],
                             epistasis=epistasis,
                             nsim=nsim)
    return beta, sebeta


def epistasis_plot(singles, double, df, **kwargs):
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

    # Calculate the point density
    points = np.vstack([X, Y])
    z = gaussian_kde(points)(points)

    # plot:
    fig, ax = plt.subplots()
    if len(X) > 50:
        ax.scatter(X, Y, c=z, s=15/np.sqrt(x.se_b**2 + y.se_b**2 + xy.se_b**2),
                   edgecolor='', cmap='viridis', alpha=0.5)
    else:
        ax.scatter(X, Y, s=15/np.sqrt(x.se_b**2 + y.se_b**2 + xy.se_b**2),
                   color='#33a02c', alpha=.9)

    smoothX = np.linspace(X.min() - 0.5, X.max() + 0.5, 1000)
    plt.plot(smoothX, -1/2*smoothX, color='#1f78b4', ls='--',
             label='Unbranched Pathway')
    plot_epistasis_regression(X, actual.beta, ls='-', lw=2.3,
                              color='#33a02c', label='fit')

    plt.xlabel(r'Predicted Additive Effect')
    plt.ylabel(r'Deviation from Additive Effect')

    plt.legend()

    return x, y, xy, ax


def curves(means, errors, labels, **kwargs):
    """Will be DEPRECATED soon."""
    def normal(x, mu, sigma):
        return 1/(2*np.pi)*np.exp(-(x-mu)**2/(2*sigma**2))

    maxe = np.max(errors)
    maxm = np.max(means)
    minm = np.min(means)

    color = kwargs.pop('color', {})
    fill = kwargs.pop('fill', False)

    X = np.linspace(minm - 3*maxe, maxm+3*maxe, 1000)

    fig, ax = plt.subplots()
    for i, mean in enumerate(means):
        if labels[i] is 'fit':
            lw = 3
        else:
            lw = 2

        if len(color) > 0:
            plot = plt.plot(X, normal(X, mean, errors[i]),
                            label=labels[i], lw=lw, color=color[labels[i]])
        else:
            plot = plt.plot(X, normal(X, mean, errors[i]),
                            label=labels[i], lw=lw)
        if not fill:
            next

        if labels[i] is 'fit':
            ax.fill_between(X, 0, normal(X, mean, errors[i]), color='r',
                            alpha=0.3)

    if minm - 3*maxe < -0.5 < maxm - 3*maxe:
        plt.gca().axvline(-0.5, color='k', ls='--', label='Unbranched Pathway')
    if minm - 3*maxe < 0 < maxm - 3*maxe:
        plt.gca().axvline(0, color='k', ls='-.', label='Additive Model')

    plt.legend(loc='upper right')
    plt.title('Epistasis Coefficient Predictions vs. Observed')
    plt.xlabel('Epistasis Coefficient')
    plt.ylabel('Density')

    return plot


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
    epistasis_choice = ['actual', 'XY=X', 'XY=Y', 'XY=X=Y', 'XY=X+Y',
                        'suppress']

    epicoef, epierr = {}, {}
    for epistasis in epistasis_choice:
        s, se_b = bootstrap_regression([x, y], xy, df,
                                       epistasis=epistasis, nsim=nsim)
        epicoef[epistasis.lower()] = s
        epierr[epistasis.lower()] = se_b

    return epicoef, epierr


def plot_bootstraps(x, y, epicoef, **kwargs):
    """Make KDE plots of the bootstrapped epistasis coefficients."""
    # make dictionaries for plotting
    colors = {'actual': '#33a02c', 'xy=x': 'blue', 'xy=y': 'k',
              'xy=x=y': '#1f78b4', 'xy=x+y': '#ff7f00', 'suppressed': '#e31a1c'
              }
    labels = {'actual': 'actual', 'xy=x': label(x, y),
              'xy=y': label(y, x), 'xy=x=y': 'Unbranched',
              'xy=x+y': 'Additive', 'suppressed': '#Suppression'
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
            next

    # plot a horizontal line wherever the actual data mean is
    plt.gca().axvline(epicoef['actual'].mean(), color='#33a02c', ls='--', lw=3)

    plt.xlabel('Epistasis Coefficient')
    plt.ylabel('Cumulative Density Function')

    return ax
