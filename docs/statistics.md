---
layout: page
title: Statistical Methods Primer
permalink: /statistics/
---

Welcome to the statistical methods section of this website. We have generated
a very brief set of tutorials for the statistical methods that we used for our
paper. Hopefully this will make the article more accessible and be helpful to
anyone who needs to replicate partially or completely some of the code we used.

## [Kallisto](http://pachterlab.github.io/kallisto/): read pseudo-alignment
We used [Kallisto](http://pachterlab.github.io/kallisto/), an excellent of
software from Lior Pachter's group, to
perform read pseudo-alignment for each mutant we analyzed. Although pseudo-alignment
is initially not as accurate as complete alignment, the speed at which the algorithm
completes means that it can be bootstrapped. The bootstrapping is what really makes
it really great (well, that and the way in which they compute k-classes). We really
love this software.

See the paper at Nature Biotechnology:
[Near-optimal probabilistic RNA-seq quantification](http://www.nature.com/nbt/journal/v34/n5/full/nbt.3519.html)

## [Sleuth](http://pachterlab.github.io/sleuth/about.html): differential expression analysis
[Sleuth](http://pachterlab.github.io/sleuth/about.html) is another great
software from Lior Pachter's group. This beautiful
library was developed to optimally accept Kallisto processed reads, although it
can work with other alignment tools. Sleuth performs the differential expression
analysis by fitting a log-linear model to explain changes in expression between
the different samples. Together, the combination of Sleuth and Kallisto are
fantastic tools for processing RNA-seq data. Equally important, they are succinct.
Kallisto is run via a single command using Terminal, and Sleuth requires about
20 lines of fairly standardized code.

See their paper at BioRxiv:
[Differential Analysis of RNA-seq incorporating quantification uncertainty](http://biorxiv.org/content/early/2016/06/10/058164)

## Bayesian robust regressions for interaction prediction

<a href="{{ site.baseurl }}/stats_tutorial/Noise Mitigation Tutorial.html"> Jupyter Notebook on Bayesian Robust Regressions</a>

## Orthogonal distance regression
Orthogonal distance regression (ODR) is a method to fit the line of best fit when you
have measurements that have errors in both the x- and the y-coordinates. Usually,
when we have measurements with no error bars, we use a method called Least Squares
to find the line of best fit. When the measurements have errors along the y-axis,
we can modify this method and use Weighted Least Squares. Weighted Least Squares
takes points with large errors and makes them less important, whereas points with
smaller errors are more important. ODR is similar in the sense that it is takes
into account errors in both x and y and weights points accordingly.

A major difference between usual Least Squares and ODR is the minimization we are performing.
For Least Squares, we are usually minimizing the vertical distance of the points
to the line. In other words, when we have minimized the equation,

$$
\sum_\mathrm{data}(y_\mathrm{data} - y__{fit}),
$$

we have found the line of best fit. However, if there are errors on both X and Y,
we need to minimize something else. In this case, it again makes sense to find the
line that also minimizes the distance between the points and the lines. However,
in this case, we will add the constraint that the ruler we use to measure this
distance must always *be at right angles to the line itself* (hence the orthogonal).

The reasoning behind this has been explained previously. For further information
on the mathematics behind regressions, we refer you to David Hogg's excellent
article on [Data Analysis recipes: Fitting a model to data](https://arxiv.org/abs/1008.4686).

In the next Jupyter notebook, we show quickly how ODR is considerably better than
least-squares when errors in both coordinates are known.

<a href="{{site.baseurl}}/stats_tutorial/ODR.html" >Jupyter Notebook on Orthogonal Regression</a>

## Bootstrapping to test a hypothesis

### Parametric Bootstrap

### Non-parametric Bootstrap


### Contact us
If there are broken links, or you have questions, you can contact us at:
[dangeles@caltech.edu](mailto:dangeles@caltech.edu)

or at

[pws@caltech.edu](mailto:pws@caltech.edu)
