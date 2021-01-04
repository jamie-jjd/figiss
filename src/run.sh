#!/bin/bash

set -x

text_dir="/home/jamie/nthu_cs/research/projects/gci/project/input/text/pizza_chili_repetitive_corpus/1mb/*"
for text in ${text_dir}
do
  ./main ${text}
done
