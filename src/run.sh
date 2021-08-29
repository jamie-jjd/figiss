#!/usr/bin/env bash

set -x

files="../data/sample.txt"
# files="../data/1mb/*"
# files="../data/corpus/*"

for file in ${files}
do
  # echo ${file}
  ./figiss ${file}
done
