#!/usr/bin/env bash

# parse seqid triplets from ortholog groups and pull out
# appropriate fasta triples
#

# local fs dependencies:
#
clustal=clustalw

my_dir=$(dirname $0)
make_tuples=$my_dir/make_fasta_tuples.py
util=$my_dir/../shared/util.bash


if [ ! $1 ]; then
  echo $0 ":: error: no genome dir"
  exit
fi

genome_dir=$1

seq_dir=$genome_dir/start_seq
out_dir=$genome_dir/fasta_tuples
clustal_out_dir=$genome_dir/clustal_tuples
ortho_set=$genome_dir/blast_nway.ortholog_set
genomes=$(cd $seq_dir; for f in *; do if [ -d $f ]; then echo $f; fi; done)


. $util

all_protein() {
  for f in $genomes; do
    echo $seq_dir/$f/protein.fa.gz
  done
}

all_tag() {
  for f in $genomes; do
    echo $f
  done
}


if [ ! -e $ortho_set ] || [ ! -s $ortho_set ]; then
  echo $0 ":: error: no ortholog set"
  exit
fi


if [ -e $out_dir ]; then 
  echo $0 ":: ortholog fasta output dir already exists. skipping..."
else
  echo $0 ":: creating fasta files for orthologs"

  if [ ! -e $out_dir ]; then mkdir -p $out_dir; fi 
  (cd $out_dir/; ls | grep -e="tuple_*.{fsa,dnd}" | xargs rm -f)

  cat $ortho_set |\
  $make_tuples -fastadb-list <(all_protein) -tag-list <(all_tag) -out $out_dir/tuple_

  # test output 
  test_dir $out_dir "fasta"
fi


if [ -e $clustal_out_dir ]; then
  echo $0 ":: clustal output_dir already exists, skipping..."
  exit
fi


echo $0 ":: starting clustal alignments"


if [ ! -e $clustal_out_dir ]; then mkdir -p $clustal_out_dir; fi
(cd $clustal_out_dir/; ls | grep -e="tuple_*.{fsa,log}" | xargs rm -f)

for f in $out_dir/tuple_*.fsa; do

  outbase=$clustal_out_dir/$(basename $f .fsa)
  $clustal $f -OUTPUT=fasta -OUTFILE=$outbase.aligned.fsa >| $outbase.log 2>&1

  if [ ! -e $outbase.aligned.fsa ] || [ ! -s $outbase.aligned.fsa ]; then
    echo $0 ":: failed clustal output: " $outbase.aligned.fsa
    exit
  fi

done


# test output 
test_dir $clustal_out_dir "clustal alignment"

