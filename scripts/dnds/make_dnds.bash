#!/usr/bin/env bash
#

# local fs dependencies:
#
bin=cq.v83

if [ ! $2 ]; then
  echo "usage: " $0 " genome_dir org_tree"
  exit
fi

genome_dir=$1
tree="$2"

rna_in_dir=$genome_dir/clustal_rna_tuples
prop_dir1=$genome_dir/clustal_rna_tuples_cats/dnds
prop_dir2=$genome_dir/clustal_rna_tuples_cats/ds

if [ ! -d $rna_in_dir ]; then
  echo $0 ":: error: no ortholog clustal rna alignments"
  exit 1
fi

is_skip_dnds=0
is_skip_ds=0

if [ -d $prop_dir1 ]; then
  echo $0 ":: $prop_dir1 property directory already exists, skipping..."
  is_skip_dnds=1
else
  mkdir -p $prop_dir1
fi

if [ -d $prop_dir2 ]; then
  echo $0 ":: $prop_dir2 property directory already exists, skipping..."
  is_skip_ds=1
else
  mkdir -p $prop_dir2
fi

if [ $is_skip_dnds == 1 ] && [ $is_skip_ds == 1 ]; then exit; fi


tmp_dir=/tmp/$USER-ds-dnds-$$

if [ -e $tmp_dir ]; then rm -rf $tmp_dir; fi
mkdir -p $tmp_dir
trap 'cd; rm -rf $tmp_dir' EXIT

site_data=$tmp_dir/data
Xds_all=$tmp_dir/Xds.all
dnds_all=$tmp_dir/dnds.all
 
>| $Xds_all
>| $dnds_all


cd $tmp_dir
for f in $rna_in_dir/*.fsa; do

  $bin -process-seq -in-seq $f -site-model codon >| $site_data 2> /dev/null

  tmr=$tmp_dir/model.report

  $bin -denovo -site-model codon -select-model single -rate-model kimura -root-model nuc-pos -lock-root -tree-file <(echo $tree) -in-data $site_data |\
  $bin -ml -in-data $site_data -in-model - 2> /dev/null |\
  $bin -report-model -in-model - >|\
  $tmr

  if [ ${PIPESTATUS[1]} != 0 ] || [ ! -e $tmr ] || [ ! -s $tmr ]; then 
    echo "WARNING:: dS & dN/dS calculation failed for $f Skipping..."
    continue
  fi

  cat $tmr |\
  awk '/^bg_time:/ {sum += $4;} END {print sum;}' >| $tmp_dir/Xds

  cat <(echo -n "$(basename $f .fsa) ") $tmp_dir/Xds >> $Xds_all
  
  cat $tmr |\
  grep "Expected dN/dS" | awk '{print $3}' >| $tmp_dir/dnds

  cat <(echo -n "$(basename $f .fsa) ") $tmp_dir/dnds >> $dnds_all
  
done
cd -



if [ $is_skip_ds == 0 ]; then

  # get ds quantile thresholds:
cat << ENDE | R --vanilla --slave | awk '{if(n==1) print $2,$3,$4; n++;}' >| $tmp_dir/Xds.quan
quantile(read.table('$Xds_all')[[2]])
ENDE

  q1=$(awk '{print $1}' $tmp_dir/Xds.quan)
  q2=$(awk '{print $2}' $tmp_dir/Xds.quan)
  q3=$(awk '{print $3}' $tmp_dir/Xds.quan)

  # label sequences by quantile:
  cat $Xds_all |\
  awk -v q1=$q1 -v q2=$q2 -v q3=$q3 '{
    if     ($2<q1) a=1;
    else if($2<q2) a=2;
    else if($2<q3) a=3;
    else           a=4;
    print $1,a;
  }' >| $tmp_dir/tmp

  mv $tmp_dir/tmp $Xds_all

  # write everything out for dnds:
  catlabel_f=$prop_dir2/catlabel
cat << ENDE >| $catlabel_f
# ML total branch length prediction by codeq v83
# divided by quartile on thresholds: $q1 $q2 $q3
#
1 lowest
2 low
3 high
4 highest
ENDE

  cat $Xds_all | while read f; do
    cat_file=$prop_dir2/$(echo $f | awk '{printf "%s.cat",$1}')

    echo "* $(echo $f | awk '{print $2}')" >| $cat_file
  done
fi



if [ $is_skip_dnds == 0 ]; then

  # get dnds quantile thresholds:
cat << ENDE | R --vanilla --slave | awk '{if(n==1) print $2,$3,$4; n++;}' >| $tmp_dir/dnds.quan
quantile(read.table('$dnds_all')[[2]])
ENDE

  q1=$(awk '{print $1}' $tmp_dir/dnds.quan)
  q2=$(awk '{print $2}' $tmp_dir/dnds.quan)
  q3=$(awk '{print $3}' $tmp_dir/dnds.quan)

  # label sequences by quantile:
  cat $dnds_all |\
  awk -v q1=$q1 -v q2=$q2 -v q3=$q3 '{
    if     ($2<q1) a=1;
    else if($2<q2) a=2;
    else if($2<q3) a=3;
    else           a=4;
    print $1,a;
  }' >| $tmp_dir/tmp

  mv $tmp_dir/tmp $dnds_all

  # write everything out for dnds:
  catlabel_f=$prop_dir1/catlabel
cat << ENDE >| $catlabel_f
# ML dN/dS prediction by codeq v83
# divided by quartile on thresholds: $q1 $q2 $q3
#
1 lowest
2 low
3 high
4 highest
ENDE

  cat $dnds_all | while read f; do
    cat_file=$prop_dir1/$(echo $f | awk '{printf "%s.cat",$1}')

    echo "* $(echo $f | awk '{print $2}')" >| $cat_file
  done
fi

