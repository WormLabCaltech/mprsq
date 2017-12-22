"""A script that contains all genotype variable information for mprsq."""


class genvars:
    """A class that contains important variables for the mprsq project."""

    def __init__(self):
        """Initialize the class object with all the necessary variables."""
        self.double_mapping = {'bd': 'a', 'bc': 'f'}

        # fancy mapping, for use in graphs
        self.fancy_mapping = {'a': r'\emph{egl-9;vhl-1}',
                              'f': r'\emph{egl-9 hif-1}',
                              'b': r'\emph{egl-9}',
                              'c': r'\emph{hif-1}',
                              'd': r'\emph{vhl-1}',
                              'e': r'\emph{rhy-1}',
                              'g': r'\emph{fog-2}'
                              }
        # mapping, for use in printing or networkx
        self.mapping = {'a': 'egl-9;vhl-1',
                        'b': 'egl-9',
                        'c': 'hif-1',
                        'd': 'vhl-1',
                        'e': 'rhy-1',
                        'f': 'egl-9 hif-1',
                        'g': 'fog-2'
                        }

        # mapping, for use in printing or networkx
        self.sort_muts = {'a': 6,
                          'b': 2,
                          'c': 4,
                          'd': 3,
                          'e': 1,
                          'f': 7,
                          'g': 5
                          }

        # sort pairs, for plotting pairwise correlations
        self.sort_pairs = {'eb': 1, 'be': 1,
                           'ed': 2, 'de': 2,
                           'ec': 3, 'ce': 3,
                           'eg': 4, 'ge': 4,
                           'bd': 5, 'db': 5,
                           'cb': 6, 'bc': 6,
                           'bg': 7, 'gb': 7,
                           'cd': 8, 'dc': 8,
                           'dg': 9, 'gd': 9,
                           'cg': 10, 'gc': 10
                           }

        # decode pairs for plotting pairwise correlations
        self.decode_pairs = {'eb': '\emph{rhy-1}, \emph{egl-9}',
                             'be': '\emph{rhy-1}, \emph{egl-9}',
                             'ed': '\emph{rhy-1}, \emph{vhl-1}',
                             'de': '\emph{rhy-1}, \emph{vhl-1}',
                             'ec': '\emph{rhy-1}, \emph{hif-1}',
                             'ce': '\emph{rhy-1}, \emph{hif-1}',
                             'eg': '\emph{rhy-1}, \emph{fog-2}',
                             'ge': '\emph{rhy-1}, \emph{fog-2}',
                             'bd': '\emph{egl-9}, \emph{vhl-1}',
                             'db': '\emph{egl-9}, \emph{vhl-1}',
                             'cb': '\emph{egl-9}, \emph{hif-1}',
                             'bc': '\emph{egl-9}, \emph{hif-1}',
                             'bg': '\emph{egl-9}, \emph{fog-2}',
                             'gb': '\emph{egl-9}, \emph{fog-2}',
                             'cd': '\emph{vhl-1}, \emph{hif-1}',
                             'dc': '\emph{vhl-1}, \emph{hif-1}',
                             'dg': '\emph{vhl-1}, \emph{fog-2}',
                             'gd': '\emph{vhl-1}, \emph{fog-2}',
                             'cg': '\emph{hif-1}, \emph{fog-2}',
                             'gc': '\emph{hif-1}, \emph{fog-2}'
                             }

        # plot order for all qpcr plots
        self.plot_order = {r'\emph{egl-9;vhl-1}': 4,
                           r'\emph{egl-9 hif-1}': 5,
                           r'\emph{egl-9}': 1,
                           r'\emph{hif-1}': 3,
                           r'\emph{vhl-1}': 2,
                           r'\emph{rhy-1}': 0,
                           r'\emph{fog-2}': 6
                           }

        # plot colors for all qpcr plots
        self.plot_color = {r'\emph{egl-9;vhl-1}': '#e41a1c',
                           r'\emph{egl-9 hif-1}': '#377eb8',
                           r'\emph{egl-9}': '#4daf4a',
                           r'\emph{hif-1}': '#984ea3',
                           r'\emph{vhl-1}': '#ff7f00',
                           r'\emph{rhy-1}': '#ffff33'
                           }


class epistasis:
    """
    An object that holds all possible epistasis models.

    Functions:
    init - initialize the object

    Attributes:
    models - a list of epistatic models
    """

    def __init__(self):
        """Initialize object."""
        self.models = ['actual', 'xy=x', 'xy=y', 'xy=x=y', 'xy=x+y',
                       'suppress']
