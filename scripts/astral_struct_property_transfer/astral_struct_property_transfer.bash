#!/usr/bin/env bash

# external fs dependencies
astral_prop_dir=$HOME/data/local/astral_prop/data/current
astral_ss3=$astral_prop_dir/astral-dssp_ss3.fsa
astral_ss8=$astral_prop_dir/astral-dssp_ss8.fsa
astral_burial=$astral_prop_dir/astral-cn14_3.fsa
astral_relative_burial=$astral_prop_dir/astral-cn14_relative_3.fsa

my_dir=$(cd $(dirname $0); pwd)
make_prop=$my_dir/make_property_aligned_fsas.py

util=$my_dir/../shared/util.bash

#
#
if [ ! $1 ] || [ ! $3 ] || (($3 < 10)) || (($3 > 99)) || [ ! $4 ] || (($4 < 10)) || (($4 > 99)); then
  echo "usage: $0 genome_dir orgid org-genome-astral-min-seq-id org-genome-astral-max-olap"
  echo " vaid seq-id,olap in 10-99"
  exit
fi

genome_dir=$1
orgid=$2
seqid=$3
olap=$4

. $util


# if no org-id use the first lex sorted org tag:
if [ $orgid == "" ]; then
  orgid=$(cd $seq_dir; for f in *; do if [ -d $f ]; then echo $f; exit; fi; done)
fi
astral_genome_alignment=$genome_dir/${orgid}_to_astral_align_seqid${seqid}_olap${olap}.fsa

if [ ! -e $astral_genome_alignment ] || [ ! -s $astral_genome_alignment ]; then
  echo $0 ":: error: no astral genome alignment"
  exit
fi

rna_dir=$genome_dir/clustal_rna_tuples

if [ ! -d $rna_dir ]; then
  echo $0 ":: error: no rna alignments"
  exit
fi


for prop in burial relative_burial ss3 ss8; do

  if [ $prop == "ss3" ]; then
    astral_prop=$astral_ss3
  elif [ $prop == "ss8" ]; then
    astral_prop=$astral_ss8
  elif [ $prop == "burial" ]; then
    astral_prop=$astral_burial
  elif [ $prop == "relative_burial" ]; then
    astral_prop=$astral_relative_burial
  fi

  if [ ! -e $astral_prop ]; then
    echo $0 ":: error: can't find astral property file, $astral_prop "
    exit
  fi

  prop_dir=$genome_dir/clustal_rna_tuples_cats/$prop

  if [ -d $prop_dir ]; then
    echo $0 ":: $prop property directory already exists, skipping..."
    continue
  else
    echo $0 ":: aligning astral domain property: $prop"
    mkdir -p $prop_dir
  fi


  $make_prop \
    -in $rna_dir \
    -out $prop_dir \
    -astral-property-file $astral_prop \
    -org-genome-astral-file $astral_genome_alignment \
    -orgid $orgid

  test_dir $prop_dir "astral $prop alignment"

done


