#!/usr/bin/env bash
#
#


# local fs dependencies:
#
my_dir=$(cd $(dirname $0); pwd)
flen=fasta_length.py

if [ $# -ne 1 ]; then
  echo "usage: $0 genome_dir"
  exit 1 
fi


genome_dir=$1

clustal_rna_dir=$genome_dir/clustal_rna_tuples
fasta_dir=$genome_dir/fasta_tuples

genome_group_label=$(echo $genome_dir | sed "s/\// /g" | awk '{print $NF}')
tmp_dir=/tmp/ctsa-$genome_group_label.$$


if [ ! -d $clustal_rna_dir ]; then
  echo $0 ":: error: no ortholog clustal rna alignments"
  exit 1
fi
if [ ! -d $fasta_dir ]; then
  echo $0 ":: error: no fasta alignments"
  exit 1
fi


rm -rf $tmp_dir
mkdir -p $tmp_dir
trap 'cd; rm -rf $tmp_dir' EXIT


prop_dir=$genome_dir/clustal_rna_tuples_cats/glength

if [ -d $prop_dir ]; then
  echo $0 ":: $prop property directory already exists, skipping..."
  exit
else
  mkdir -p $prop_dir
fi

tmp_prop_file=$tmp_dir/glen_all

>| $tmp_prop_file


count=0
for f in $(ls $clustal_rna_dir | grep -e="*.fsa"); do

  ftag=$(basename $f .aligned.fsa)

  val=$($flen < $fasta_dir/$ftag.fsa | awk '{if($1<a||a==0)a=$1;} END {print a;}')

  if [ $val"" != "" ]; then echo "$ftag $val" >> $tmp_prop_file; fi

  ((count++))

  if ! ((count%10)); then echo -n "."; fi
  if ! ((count%800)); then echo ""; fi

done

echo ""


# get quantile thresholds:
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
# min aa length in alignment
# divided by quartile on thresholds: $q1 $q2 $q3
#
1 lowest
2 low
3 high
4 highest
ENDE

# label sequences by quartile:
#
cat $tmp_prop_file |\
while read line; do
  set $line
  prop_file=$prop_dir/$1.aligned.cat

  echo $line |\
  awk -v q1=$q1 -v q2=$q2 -v q3=$q3 '{
    if     ($2<q1) a=1;
    else if($2<q2) a=2;
    else if($2<q3) a=3;
    else           a=4;
    print "*",a;
  }' >| $prop_file 

done

