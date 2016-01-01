#!/usr/bin/env bash
#
#


# local fs dependencies:
#
my_dir=$(cd $(dirname $0); pwd)
cpgf=$my_dir/cpgf_S.py


if [ $# -ne 2 ]; then
  echo "usage: $0 genome_dir org_tag"
  exit 1 
fi


genome_dir=$1
orgid=$2

clustal_rna_dir=$genome_dir/clustal_rna_tuples

genome_group_label=$(echo $genome_dir | sed "s/\// /g" | awk '{print $NF}')
tmp_dir=/tmp/ctsa-$genome_group_label.$$


if [ ! -d $clustal_rna_dir ]; then
  echo $0 ":: error: no ortholog clustal rna alignments"
  exit 1
fi


rm -rf $tmp_dir
mkdir -p $tmp_dir
trap 'cd; rm -rf $tmp_dir' EXIT


prop_dir=$genome_dir/clustal_rna_tuples_cats/cpgisland_$orgid
#prop_dir=$my_dir/clustal_rna_tuples_cpgisland_$orgid

if [ -d $prop_dir ]; then
  echo $0 ":: $prop property directory already exists, skipping..."
  exit
else
  mkdir -p $prop_dir
fi


catlabel_f=$prop_dir/catlabel
cat << ENDE >| $catlabel_f
# max scoring segment cpg islands: +17/-1 S=50, 5' only
#
1 cpgisland_$orgid
ENDE


count=0

for f in $(ls $clustal_rna_dir | grep -e="*.fsa"); do

  ftag=$(basename $f .fsa)

  cat $clustal_rna_dir/$f |\
  awk -v orgid=$orgid '{if(substr($0,1,1)==">"){printit=0; if($2==orgid) printit=1;} if(printit) print;}' |\
  $cpgf -only5 >|\
  $tmp_dir/tmp.cat

  if [ -s $tmp_dir/tmp.cat ]; then
    mv $tmp_dir/tmp.cat $prop_dir/$ftag.cat
  fi

  ((count++))

  if ! ((count%10)); then echo -n "."; fi
  if ! ((count%800)); then echo ""; fi

done


