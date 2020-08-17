#!/bin/bash

set -x

make

rm -f ../output/sample.csv
./main ../input/text/sample.txt ../output/sample.csv

# input_root_directory="../input/text/pizza_chili_repetitive_corpus/raw_corpus"
# output_root_directory="../output/csv"
#
# for directory in `ls $input_root_directory`
# do
#   for file in `ls "$input_root_directory/$directory"`
#   do
#     input_file="$input_root_directory/$directory/$file"
#     output_file="$output_root_directory/$file.csv"
#     rm -f $output_file
#     ./main $input_file $output_file
#   done
# done
