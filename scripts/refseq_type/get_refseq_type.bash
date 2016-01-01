#!/usr/bin/env bash
#

# local fs dependencies:
#
my_dir=$(cd $(dirname $0); pwd)


if [ $# -ne 1 ]; then
  echo "usage: $0 genome_dir"
  exit 1 
fi


genome_dir=$1

clustal_rna_dir=$genome_dir/clustal_rna_tuples
clustal_dir=$genome_dir/clustal_tuples


if [ ! -d $clustal_rna_dir ]; then
  echo $0 ":: error: no ortholog clustal rna alignments"
  exit 1
fi


prop_dir=$genome_dir/clustal_rna_tuples_cats/refseq_type

if [ -d $prop_dir ]; then
  echo $0 ":: $prop property directory already exists, skipping..."
  exit
else
  mkdir -p $prop_dir
fi


catlabel_f=$prop_dir/catlabel
cat << ENDE >| $catlabel_f
# divide refseq alignments into those entirely composed of curated
# sequences and all others
#
1 noncurated 
2 curated
ENDE


for f in $(ls $clustal_rna_dir | grep -e="*.fsa"); do

  ftag=$(basename $f .fsa)
  outfile=$prop_dir/$ftag.cat

  val=1
  if awk 'BEGIN {FS="|"} {if(substr($0,1,1)==">" && substr($4,1,3)!="NP_") exit 1;}' $clustal_dir/$f; then
    val=2
  fi
  echo "* $val" >| $outfile
done

