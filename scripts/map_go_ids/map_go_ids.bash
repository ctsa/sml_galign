#!/bin/bash
#
# use ncbi gbk rna files to find list of categories associated with each protein gi number
#

# local fs dependencies:
#
my_dir=$(dirname $0)
make_raw_golist=$my_dir/make_ncbi_pro_raw_go_map.py
make_neat_golist=$my_dir/make_ncbi_pro_neat_go_map.py
util=$my_dir/../shared/util.bash
sdir=/home/ctsa/data/pub/go/current/ontology


if [ ! $1 ]; then
  echo "$0 :: no genome dir"
  exit
fi

genome_dir=$1

seq_dir=$genome_dir/start_seq
rna_dir=$genome_dir/clustal_rna_tuples

genomes=$(cd $seq_dir; for f in *; do if [ -d $f ]; then echo $f; fi; done)

. $util


all_rna_gbk() {
  for f in $genomes; do
    echo $seq_dir/$f/rna.gbk.gz
  done
}


all_org_tag() {
  for f in $genomes; do echo $f; done
}


if [ ! -d $rna_dir ]; then
  echo $0 ":: no rna dir detected. skipping..."
  exit
fi


go_out=$genome_dir/go_categories
raw_go_out=$go_out.raw
if ! [ -e $raw_go_out ]; then rm $raw_go_out; fi

$make_raw_golist -gbk-list <(all_rna_gbk) -rnadir $rna_dir -out $raw_go_out

if [ $? -ne 0 ]; then
  echo "$0 :: error in go category script. removing file"
  rm $raw_go_out
  exit
fi

neat_go_out=$go_out.neat

for f in function process component; do
  $make_neat_golist -in $raw_go_out -squish $sdir/$f.ontology.8Jan2008.squish -out $neat_go_out.$f
done

