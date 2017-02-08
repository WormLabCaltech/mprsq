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
<a href="{{site.baseurl}}/stats_tutorial/ODR.html" >Jupyter Notebook on Orthogonal Regression</a>

## Bootstrapping to test a hypothesis

### Parametric

### Non-parametric


### Contact us
If there are broken links, or you have questions, you can contact us at:
[dangeles@caltech.edu](mailto:dangeles@caltech.edu)

or at

[pws@caltech.edu](mailto:pws@caltech.edu)
