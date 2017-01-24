#!/bin/sh

cd "$(dirname "$0")"

echo "Running R based pipeline"

# do data preparation for DESeq2 analysis in R
python3 read_count_class.py -R

# do DESeq2 analysis using all volunteers
# non-responders are all assigned factors as if they have been measured on day 0
Rscript deseq2_NDN.R

rm Rplots.pdf
