#!/bin/bash

set -x

files="../data/sample.txt"
# files="../data/corpus/*"

for file in ${files}
do
  # echo ${file}
  ./main ${file}
done
