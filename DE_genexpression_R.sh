#!/bin/sh

# Write bash script that runs complete analysis of GOA gene-expression data using the full DESeq tactic.
# Approach tactic: give all non-responders day 0 as day factor.

cd "$(dirname "$0")"

# python does data preparation for DESeq2 analysis in R
python3 readcountclass.py -R

# does DESeq2 analysis using all volunteers
# non-responders are all assigned factors as if they have been measured on day 0
# Rscript deseq2_NDN.R

