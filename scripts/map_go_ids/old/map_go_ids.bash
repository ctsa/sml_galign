#!/usr/bin/env bash
#
# map all possible go ids to each alignment, based on the annotations provided in the ncbi rna genbank file for each organism
#
#


if [ ! $1 ]; then
  echo $0 ":: error: no genome dir"
  exit
fi


genome_dir=$1

seq_dir=$genome_dir/start_seq
out_dir=$genome_dir/fasta_tuples
clustal_out_dir=$genome_dir/clustal_rna_tuples
ortho_set=$genome_dir/blast_nway.ortholog_set
genomes=$(cd $seq_dir; for f in *; do if [ -d $f ]; then echo $f; fi; done)

genome_group_label=$(echo $genome_dir | sed "s/\// /g" | awk '{print $NF}')
tmp_dir=/tmp/ctsa-$genome_group_label.$$


# local dependencies
#
my_dir=$(dirname $0)
#mid=$my_dir/map_go_ids.py
#go_assoc_dir=$HOME/data/pub/go/current/association


all_tag() {
  for f in $genomes; do
    echo $f
  done
}



if [ -e $tmp_dir ]; then rm -rf $tmp_dir; fi
mkdir -p $tmp_dir
trap 'cd; rm -rf $tmp_dir' EXIT


>| $tmp_dir/tag_file

for f in $(all_tag); do
  echo $f $seq_dir/$f/rna.gbk.gz >> $tmp_dir/tag_file
done

$mid -in $tmp_dir/tag_file -clustal-rna-dir $clustal_out_dir >| ${clustal_out_dir}_GO

