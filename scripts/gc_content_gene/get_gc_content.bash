#!/usr/bin/env bash
#
#


# local fs dependencies:
#
my_dir=$(cd $(dirname $0); pwd)
gc_gene=$my_dir/gc_gene.py

is_3c=0
is_3cf=0

if [ $# -gt 2 ]; then
  echo "usage: $0 [ -3c ] genome_dir"
  exit 1 
fi

if [ $2"" == "-3c" ]; then 
  is_3c=1 
elif [ $2"" == "-3cf" ]; then
  is_3cf=1
fi


genome_dir=$1

clustal_rna_dir=$genome_dir/clustal_rna_tuples

genome_group_label=$(echo $genome_dir | sed "s/\// /g" | awk '{print $NF}')
tmp_dir=/tmp/ctsa-$genome_group_label-gc_gene-$$


if [ ! -d $clustal_rna_dir ]; then
  echo $0 ":: error: no ortholog clustal rna alignments"
  exit 1
fi


rm -rf $tmp_dir
mkdir -p $tmp_dir
trap 'cd; rm -rf $tmp_dir' EXIT

gtag=""
if [ $is_3c == 1 ]; then gtag="_3c"; fi
if [ $is_3cf == 1 ]; then gtag="_3cf"; fi

prop_dir=$genome_dir/clustal_rna_tuples_cats/gc_gene$gtag

if [ -d $prop_dir ]; then
  echo $0 ":: $prop property directory already exists, skipping..."
  exit
else
  mkdir -p $prop_dir
fi

tmp_prop_file=$tmp_dir/gc_all

>| $tmp_prop_file


gflag=""
if [ $is_3c == 1 ]; then gflag="-3c"; fi
if [ $is_3cf == 1 ]; then gflag="-3cf"; fi

count=0
for f in $(ls $clustal_rna_dir | grep -e="*.fsa"); do

  ftag=$(basename $f .fsa)

  val=$($gc_gene $gflag < $clustal_rna_dir/$f)

  if [ $val"" != "" ]; then echo "$ftag $val" >> $tmp_prop_file; fi

  ((count++))

  if ! ((count%10)); then echo -n "."; fi
  if ! ((count%800)); then echo ""; fi

done


# get gc quantile thresholds:
set $(
cat << ENDE | R --vanilla --slave | awk '{if(n==1) print $2,$3,$4; n++;}'
quantile(read.table('$tmp_prop_file')[[2]])
ENDE
)

q1=$1
q2=$2
q3=$3

# write out catlabel file:
catlabel_file=$prop_dir/catlabel
cat << ENDE >| $catlabel_file
# full alignment G+C%
# divided by quartile on thresholds: $q1 $q2 $q3
#
1 lowest
2 low
3 high
4 highest
ENDE

# label gc sequences by quartile:
#
cat $tmp_prop_file |\
while read line; do
  set $line
  prop_file=$prop_dir/$1.cat

  echo $line |\
  awk -v q1=$q1 -v q2=$q2 -v q3=$q3 '{
    if     ($2<q1) a=1;
    else if($2<q2) a=2;
    else if($2<q3) a=3;
    else           a=4;
    print "*",a;
  }' >| $prop_file 

done

