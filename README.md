# Reconstructing a metazoan genetic pathway with transcriptome-wide epistasis measurements
#### Authors: David Angeles-Albores, Carmie Puckett-Robinson, Brian A. Williams, Barbara Wold and Paul W. Sternberg

## Abstract
RNA-seq is commonly used to identify genetic modules that respond to perturbations. In single cells, transcriptomes have been used as phenotypes, but this concept has not been applied to whole-organism RNA-seq. Also, quantifying and interpreting epistatic effects using expression profiles remains a challenge. We developed a single coefficient to quantify transcriptome-wide epistasis that reflects the underlying interactions and which can be interpreted intuitively. To demonstrate our approach, we sequenced four single and two double mutants of *Caenorhabditis elegans*. From these mutants, we reconstructed the known hypoxia pathway. In addition, we uncovered a class of 56 genes with *hif-1*-dependent expression which have opposite changes in expression in mutants of two genes which cooperate to negatively regulate HIF-1 abundance; however, the double mutant of these genes exhibits suppression epistasis. This class violates the classical model of HIF-1 regulation, but can be explained by postulating a role of hydroxylated HIF-1 in transcriptional control.

# File Structure
```
mprsq
│   README.md
│   kallisto_commands_sh    
│
└───input -- input files used for this project. Note: due to their size, HD5
             files are gitignored, but can be found at the GEO.
└───experimental_docs - experimental documents for this project
└───sleuth_all_adjusted - contains R code for sleuth processing
└───src - contains all Jupyter notebooks with python analysis
└───output
      └─── supplementary_tables/ -- contains CSV files containing quantified results
      └─── rank_plots/ -- pairwise rank plots for all pairwise combinations
      └─── other_figs -- all other figures generated for this project
└───docs - website
```
## Contact
pws@caltech.edu

## Technical assistance
dangeles@caltech.edu

## Acknowledgements
