"""
Functions for pretty-printing tables.
author: dangeles at caltech edu
"""


def table_print(l, space=20):
    """
    A function to pretty-print tables.

    Params:
    l - a list of strings, each entry is in a different column
    space - the space between the separators.
            Must be > than the longest entry in l

    Output:
    a print statement
    """
    # if an entry is longer than the default spacing, change the spacing.
    for i in l:
        string = str(i)
        if len(string) > space:
            space = len(string) + 5
    # assemble the message and print
    s = ''
    for i in l:
        i = str(i)
        s += i + ' '*(space-len(i))
    print(s)


def significance(pval, alpha, x, y, test_statistic):
    """
    A function to print out for statistical significance.

    Params:
    pval - pvalue
    alpha - threshold for statistical significance
    x, y - the names of the two samples tested for equality.
    test_statistic - the (plural version) of the statistic of interest

    Output:
    a print statement
    """
    if pval < alpha:
        m = '{0} and {1} have statistically significantly different {2}.'
        print(m.format(x, y, test_statistic))

    print('The p-value for this test is {0:.2g}'.format(pval))
