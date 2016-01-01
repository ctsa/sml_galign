#!/usr/bin/env bash

# align pdb to sister1 genome
#

# local fs dependencies:
#
# in path, original source on genome is $HOME/proj/subs_ml/disorder/holladay_disorder/disorder_prediction
pdis=predict_disorder.exe

my_dir=$(cd $(dirname $0); pwd)
mac=$my_dir/make_aligned_cat.py
disorder_scorefile=$HOME/proj/subs_ml/cat_bio_studies/disorder/holladay_disorder/disorder_prediction/PredictorParameterFiles/sw35_st30_5_train_avg_params_700.xml


#
#
if [ ! $1 ] || [ ! $2 ]; then
  echo "usage: $0 genome_dir org-id"
  exit
fi

genome_dir=$1
orgid=$2
thresh=2.0

# if no org-id use the first lex sorted org tag:
if [ $orgid == "" ]; then
  orgid=$(cd $seq_dir; for f in *; do if [ -d $f ]; then echo $f; exit; fi; done)
fi
ortho_dir=$genome_dir/fasta_tuples
clustal_dir=$genome_dir/clustal_tuples
clustal_rna_dir=$genome_dir/clustal_rna_tuples

genome_group_label=$(echo $genome_dir | sed "s/\// /g" | awk '{print $NF}')
tmp_dir=/tmp/ctsa-$genome_group_label.$$


get_org_seqs() {
  cd $ortho_dir
  ls |\
  grep -e="*.fsa" |\
  xargs cat |\
  awk -v orgid=$orgid '{if(substr($0,1,1)==">"){printit=0; if($2==orgid) printit=1;} if(printit) print;}'
}


if [ ! -d $ortho_dir ]; then
  echo $0 ":: error: no ortholog sequences"
  exit
fi
if [ ! -d $clustal_dir ]; then
  echo $0 ":: error: no ortholog clustal alignments"
  exit
fi
if [ ! -d $clustal_rna_dir ]; then
  echo $0 ":: error: no rna alignments"
  exit
fi



rm -rf $tmp_dir
mkdir -p $tmp_dir
trap 'cd; rm -rf $tmp_dir' EXIT



prop_dir=$genome_dir/clustal_rna_tuples_cats/disorder

if [ -d $prop_dir ]; then
   echo $0 ":: $prop property directory already exists, skipping..."
   exit
else
   mkdir -p $prop_dir
fi



catlabel_f=$prop_dir/catlabel
cat << ENDE >| $catlabel_f
# holladay & grishin sequence based disorder prediction
# with option: -tails
# disordered threshold >= 2.0
#
1 ordered
2 disordered
ENDE



(
cd $pdis_path
for f in $(ls $ortho_dir | grep -e="*.fsa"); do

  # check that rna alignment exists (as a quality filter)
  clustal_rna_f=$clustal_rna_dir/$(basename $f .fsa).aligned.fsa
  if [ ! -e $clustal_rna_f ]; then continue; fi

  cat $ortho_dir/$f |\
  awk -v orgid=$orgid '{if(substr($0,1,1)==">"){printit=0; if($2==orgid) printit=1;} if(printit) print;}' |\
  $pdis -tails -params $disorder_scorefile |\
  awk -v thresh=$thresh '{if(substr($0,1,1)!=">") {if ($3>=thresh) print $1,$2,2; else print $1,$2,1;}}' >|\
  $tmp_dir/disorder_tails 

  if [ ${PIPESTATUS[2]} != 0 ]; then 
    echo "WARNING:: disorder prediction failed for $f Skipping..."
    continue
  fi

  orgtag=$(cat $ortho_dir/$f | awk -v orgid=$orgid '{if(substr($0,1,1)==">" && $2==orgid) print $1;}')

  clustal_f=$clustal_dir/$(basename $f .fsa).aligned.fsa

  cat $clustal_f |\
  awk -v orgtag=$orgtag '{if(substr($0,1,1)==">"){printit=0; if($1==orgtag) printit=1;} if(printit) print;}' >|\
  $tmp_dir/adjusted_seq.fsa

  prop_f=$prop_dir/$(basename $f .fsa).aligned.cat

  $mac $tmp_dir/adjusted_seq.fsa $tmp_dir/disorder_tails >|\
  $prop_f
done
)

