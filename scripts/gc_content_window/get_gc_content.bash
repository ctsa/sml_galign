#!/usr/bin/env bash
#
#


# local fs dependencies:
#
my_dir=$(cd $(dirname $0); pwd)
gc_win=$my_dir/gc_win.py
winsize=200


if [ $# -ne 1 ]; then
  echo "usage: $0 genome_dir"
  exit 1 
fi

genome_dir=$1

clustal_rna_dir=$genome_dir/clustal_rna_tuples

genome_group_label=$(echo $genome_dir | sed "s/\// /g" | awk '{print $NF}')
tmp_dir=/tmp/ctsa-$genome_group_label-gc_window-$$


if [ ! -d $clustal_rna_dir ]; then
  echo $0 ":: error: no ortholog clustal rna alignments"
  exit 1
fi


rm -rf $tmp_dir
mkdir -p $tmp_dir
trap 'cd; rm -rf $tmp_dir' EXIT


prop_dir=$genome_dir/clustal_rna_tuples_cats/gc_win

if [ -d $prop_dir ]; then
  echo $0 ":: $prop property directory already exists, skipping..."
  exit
else
  mkdir -p $prop_dir
fi

tmp_prop_dir=$tmp_dir/tmp_prop
mkdir -p $tmp_prop_dir


count=0
for f in $(ls $clustal_rna_dir | grep -e="*.fsa"); do

  ftag=$(basename $f .fsa)
  tmp_prop_file=$tmp_prop_dir/$ftag.gc

  cat $clustal_rna_dir/$f |\
  $gc_win -c -w $winsize >|\
  $tmp_prop_file

  ((count++))

  if ! ((count%10)); then echo -n "."; fi
  if ! ((count%800)); then echo ""; fi

done


gc_all=$tmp_dir/gc.all
cd $tmp_prop_dir
ls | grep -e="*.gc" | xargs cat >| $gc_all
cd - > /dev/null

# get gc quantile thresholds:
set $(
cat << ENDE | R --vanilla --slave | awk '{if(n==1) print $2,$3,$4; n++;}'
quantile(read.table('$gc_all')[[2]])
ENDE
)

q1=$1
q2=$2
q3=$3

# write out catlabel file:
catlabel_file=$prop_dir/catlabel
cat << ENDE >| $catlabel_file
# alignment G+C content, window: +/- $winsize bp
# divided by quartile on thresholds: $q1 $q2 $q3
#
1 lowest
2 low
3 high
4 highest
ENDE

# label gc sequences by quantile:
#
for f in $(ls $tmp_prop_dir); do
  gc_file=$tmp_prop_dir/$f
  prop_file=$prop_dir/$(basename $f .gc).cat

  cat $gc_file |\
  awk -v q1=$q1 -v q2=$q2 -v q3=$q3 '{
    if     ($2<q1) a=1;
    else if($2<q2) a=2;
    else if($2<q3) a=3;
    else           a=4;
    print $1,a;
  }' >| $prop_file 

done

