#!/bin/bash
# Bash commands for diff. expression analysis using Sleuth.
# Sleuth analysis command for WT_a
Rscript diff_exp_analyzer.R -d WT_a --genovar za
# Sleuth analysis command for WT_b
Rscript diff_exp_analyzer.R -d WT_b --genovar zb
# Sleuth analysis command for WT_c
Rscript diff_exp_analyzer.R -d WT_c --genovar zc
# Sleuth analysis command for WT_d
Rscript diff_exp_analyzer.R -d WT_d --genovar zd
# Sleuth analysis command for WT_e
Rscript diff_exp_analyzer.R -d WT_e --genovar ze
# Sleuth analysis command for WT_f
Rscript diff_exp_analyzer.R -d WT_f --genovar zf
