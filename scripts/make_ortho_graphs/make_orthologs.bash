#!/usr/bin/env bash

# local fs dependencies:
#
orth_bin=orthologer


if [ ! $1 ]; then 
  echo "$0 :: error: no genome dir"
  exit
fi

genome_dir=$1

out_file=$genome_dir/blast_nway.ortholog_set

if [ -e $out_file ]; then
  echo "$0 :: orthologer output already exists. skipping..."
  exit
fi


cat $genome_dir/blast_nway.filtered |\
$orth_bin >|\
$out_file

if [ ! -s $out_file ]; then
  echo "$0 :: error: orthologer run failed"
  rm -f $out_file
fi
