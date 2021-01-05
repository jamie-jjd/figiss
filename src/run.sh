#!/bin/bash

set -x

text_dir="/home/jamie/nthu_cs/research/projects/gci/project/input/text/benchmark_corpus/*"
for text in ${text_dir}
do
  ./main ${text}
done

rm -f ../input/pattern/*
