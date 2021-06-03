#!/bin/bash

set -x

prefix="/home/jamie/nthu_cs/research/projects/gci/project/data"

# category_prefix="statistics/run_length_wavelet_tree"

category="corpus"
# category="index"

size="1mb/*"
# size="40mb/*"
# size="raw/*"

dir=${prefix}/${category}/${size}

for entry in ${dir}
do
  # echo ${entry}
  ./main ${entry}
done

# ./main ../data/sample.txt
# ./main ${prefix}/${category}/1mb/kernel.1mb
