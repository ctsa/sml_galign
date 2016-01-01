#!/usr/bin/env bash
#
# parse gi tuples from ortholog groups and pull out
# appropriate fasta tuples
#

# local fs dependencies:
#
my_dir=$(dirname $0)
make_rna=$my_dir/make_ncbi_pro_rna_map.py
make_rna_fsa=$my_dir/make_fsa_pro_rna_map.py
util=$my_dir/../shared/util.bash


if [ ! $1 ]; then
  echo "$0 :: no genome dir"
  exit
fi

genome_dir=$1
is_yeast_format=$2

seq_dir=$genome_dir/start_seq
clustal_out_dir=$genome_dir/clustal_tuples
rna_out_dir=$genome_dir/clustal_rna_tuples

genomes=$(cd $seq_dir; for f in *; do if [ -d $f ]; then echo $f; fi; done)


. $util


all_rna_gbk() {
  for f in $genomes; do
    echo $seq_dir/$f/rna.gbk.gz
  done
}


all_rna_fsa() {
  for f in $genomes; do
    echo $seq_dir/$f/orf_rna.fa.gz
  done
}


all_org_tag() {
  for f in $genomes; do echo $f; done
}


if [ ! -d $clustal_out_dir ]; then
  echo $0 ":: error: no clustal output"
  exit
fi

if [ -e $rna_out_dir ]; then
  echo $0 ":: rna output dir detected. skipping..."
  exit
fi

if [ ! -e $rna_out_dir ]; then mkdir -p $rna_out_dir; fi
(cd $rna_out_dir/; ls | grep -e="tuple_*.fsa" | xargs rm -f)

if [ $is_yeast_format == 0 ]; then
  $make_rna -gbk-list <(all_rna_gbk) -tag-list <(all_org_tag) \
    -indir $clustal_out_dir -outdir $rna_out_dir
else
  $make_rna_fsa -fsa-list <(all_rna_fsa) -tag-list <(all_org_tag) \
    -indir $clustal_out_dir -outdir $rna_out_dir
fi

if [ $? -ne 0 ]; then
  echo "$0 :: error in clustal_rna script. removing directory"
  rm -rf $rna_out_dir
else
  test_dir $rna_out_dir "clustal rna"
fi
