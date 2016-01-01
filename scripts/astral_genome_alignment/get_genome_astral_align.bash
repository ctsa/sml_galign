#!/usr/bin/env bash

# align pdb to sister1 genome
#

# local fs dependencies:
#
data_dir=$HOME/data/pub
astral_bdb=$data_dir/astral/current/fasta/astral.fsa
pdbres_file=$data_dir/pdb/derived_data/resolu.idx
ncbi_data_dir=$HOME/opt/noarch/blast_data

psiblast=blastpgp
fasta_length_filter=fasta_length_filter.py
fasta_split=fasta_split.py
formatdb=formatdb

my_dir=$(cd $(dirname $0); pwd)
pick_best=$my_dir/pick_best_pdb.py


#
#
if [ ! $1 ] || [ ! $3 ] || (($3 < 10)) || (($3 > 99)) || [ ! $4 ] || (($4 < 10)) || (($4 > 99)); then
  echo "usage: $0 genome_dir org-id org-genome-astral-min-seq-id org-genome-astral-max-olap"
  echo " vaid seq-id,olap in 10-99"
  exit
fi

genome_dir=$1
orgid=$2
seqid=$3
olap=$4

# if no org-id use the first lex sorted org tag:
if [ $orgid == "" ]; then
  orgid=$(cd $seq_dir; for f in *; do if [ -d $f ]; then echo $f; exit; fi; done)
fi
ortho_dir=$genome_dir/fasta_tuples
blastoptions="-J T -I T -m 7 -v 10 -b 10 -s T"
output_tag=${orgid}_to_astral_align_seqid${seqid}_olap$olap
genome_group_label=$(echo $genome_dir | sed "s/\// /g" | awk '{print $NF}')
tmp_dir=/tmp/ctsa-$genome_group_label.$$
split_dir=$tmp_dir/split


get_org_seqs() {
  cd $ortho_dir
  ls |\
  grep -e="*.fsa" |\
  xargs cat |\
  awk -v orgid=$orgid '{if(substr($0,1,1)==">"){printit=0; if($2==orgid) printit=1;} if(printit) print;}'
}

if [ ! -e $astral_bdb.psq ]; then
  $formatdb -i $astral_bdb
fi

if [ ! -d $ortho_dir ]; then
  echo $0 ":: error: no ortholog alignments"
  exit
fi

if [ -e $genome_dir/$output_tag.fsa ] && [ -s $genome_dir/$output_tag.fsa ]; then
  echo $0 ":: pdb align output detected, skipping..."
  exit
fi


rm -rf $tmp_dir
mkdir -p $tmp_dir
trap 'cd; rm -rf $tmp_dir' EXIT

cd $tmp_dir
echo "[NCBI]" >| .ncbirc
echo "Data=$ncbi_data_dir" >> .ncbirc

# split up input fasta into sane chunks AND move
# onto local disk:
#
if [ -e $split_dir ]; then rm -rf $split_dir; fi

get_org_seqs |\
$fasta_length_filter -min 40 |\
$fasta_split -out $split_dir -n 1000



(
  >| $genome_dir/$output_tag.fsa
  for splitfile in $split_dir/*.fsa; do
    cat $splitfile |\
    $psiblast $blastoptions -d $astral_bdb |\
    $pick_best -seqid 0.$seqid -length 40 -eval 0.01 -overlap 0.$olap >>\
    $genome_dir/$output_tag.fsa
  done
) 2>| $genome_dir/$output_tag.log

# -pdbres 3.9 -allow-nmr -resfile $pdbres_file


