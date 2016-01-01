#!/usr/bin/env bash

# convert phobius short form to binary categories (transm or not)
#

# local fs dependencies:
#

my_dir=$(cd $(dirname $0); pwd)
usf=$my_dir/unfold_short_form.py
phob=jphobius



if [ $# -ne 2 ]; then
  echo "usage: $0 genome_dir org_id"
  exit
fi


genome_dir=$1
orgid=$2

ortho_dir=$genome_dir/fasta_tuples
clustal_dir=$genome_dir/clustal_tuples
clustal_rna_dir=$genome_dir/clustal_rna_tuples

genome_group_label=$(echo $genome_dir | sed "s/\// /g" | awk '{print $NF}')
tmp_dir=/tmp/ctsa-$genome_group_label.$$


if [ ! -d $ortho_dir ]; then
  echo $0 ":: error: no ortholog sequences"
  exit 1
fi
if [ ! -d $clustal_dir ]; then
  echo $0 ":: error: no ortholog clustal alignments"
  exit 1
fi
if [ ! -d $clustal_rna_dir ]; then
  echo $0 ":: error: no ortholog rna alignments"
  exit 1
fi


rm -rf $tmp_dir
mkdir -p $tmp_dir
trap 'cd; rm -rf $tmp_dir' EXIT


prop_dir=$genome_dir/clustal_rna_tuples_cats/transmem

if [ -d $prop_dir ]; then
  echo $0 ":: property directory already exists, skipping..."
  exit
else
  mkdir -p $prop_dir
fi

cat << ENDE >| $prop_dir/catlabel
# phobius transmembrane predictions
#
1 other
2 transmembrane
ENDE


count=0

for f in $(ls $ortho_dir | grep -e="*.fsa"); do

  fasta_file=$ortho_dir/$f
  clustal_file=$clustal_dir/$(basename $f .fsa).aligned.fsa
  clustal_rna_file=$clustal_rna_dir/$(basename $f .fsa).aligned.fsa
  prop_file=$prop_dir/$(basename $f .fsa).aligned.cat

  # check that rna alignment exists (as a quality filter)
  #
  if [ ! -e $clustal_rna_file ]; then continue; fi

  # get clustal aligned aa seq for the right organism
  gi_num=$(cat $fasta_file |\
           awk -v orgid=$orgid '{if(substr($0,1,1)==">") if($2==orgid) printit=1; else printit=0;  if(printit) print;}' |\
           sed "s/|/ /g" |\
           awk '{if(substr($0,1,1)==">") print $2;}')

  cat $clustal_file |\
  awk -v gi_num=$gi_num '{if(substr($0,1,1)==">"){printit=0; split($0,a,"|"); if(a[2]==gi_num) printit=1;} if(printit) print;}' >|\
  $tmp_dir/adjusted_seq.fsa

  # get ortho fsa's for the selected organism, with header striped down to gi number only
  #
  cat $fasta_file |\
  awk -v orgid=$orgid '{if(substr($0,1,1)==">") if($2==orgid) printit=1; else printit=0;  if(printit) print;}' |\
  sed "s/|/ /g" |\
  awk '{if(substr($0,1,1)==">") print ">",$2; else print;}' |\
  $phob -s 2> /dev/null |\
  tail +2 |\
  $usf $tmp_dir/adjusted_seq.fsa >|\
  $prop_file

  if [ ! -s $prop_file ]; then rm $prop_file; fi

  ((count++))

  if ! ((count%10)); then echo -n "."; fi
  if ! ((count%800)); then echo ""; fi
done


