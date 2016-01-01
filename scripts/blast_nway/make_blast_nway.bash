#!/usr/bin/env bash

# run all proteins from 2 sister genomes and an outgroup against each other to
# find ortholog cliques
#

# local fs dependencies:
#
ncbi_data_dir=$HOME/opt/noarch/blast_data
blast=blastall
formatdb=formatdb
fasta_split=fasta_split.py


#
#
if [ ! $1 ]; then
  echo "no genome dir specified"
  exit
fi 


genome_dir=$1

seq_dir=$genome_dir/start_seq
blastdb_dir=$genome_dir/blastdb
out_dir=$genome_dir/blast_nway

genome_group_label=$(echo $genome_dir | sed "s/\// /g" | awk '{print $NF}')
tmp_dir=/tmp/ctsa-$genome_group_label.$$
split_dir=$tmp_dir/split


# check that genome group has been set up
#
if [ ! -e $seq_dir ]; then
  echo "Cannot find input seq data directory:" $seq_dir
  exit
fi

# get genome codes from the genome code list
# \todo - just read this from the folders
#
genomes=$(cd $seq_dir; for f in *; do if [ -d $f ]; then echo $f; fi; done)


# first make sure blast databases have been built for
# each genome
#
if [ ! -e $blastdb_dir ]; then mkdir -p $blastdb_dir; fi

for f in $genomes; do
  if [ -e $blastdb_dir/$f ] && [ -e $blastdb_dir/$f/protein.fa.psq ]; then continue; fi
  if [ ! -e $blastdb_dir/$f ]; then mkdir $blastdb_dir/$f; fi
  (
    echo "Building blastdb for group: $f"
    cd $blastdb_dir/$f
    mkfifo protein.fa
    gzip -dc $seq_dir/$f/protein.fa.gz >| protein.fa &
    $formatdb -i protein.fa -o T
    rm protein.fa
  )
done



if [ -e $tmp_dir ]; then rm -rf $tmp_dir; fi
mkdir -p $tmp_dir
trap 'cd; rm -rf $tmp_dir' EXIT

if [ ! -e $out_dir ]; then mkdir -p $out_dir; fi

cd $tmp_dir
echo "[NCBI]" >| .ncbirc
echo "Data=$ncbi_data_dir" >> .ncbirc


for f in $genomes; do
  nfs_infsa=$seq_dir/$f/protein.fa.gz

  # split up input fasta into sane chunks AND move
  # onto local disk:
  #
  if [ -e $split_dir ]; then rm -rf $split_dir; fi

  gzip -dc $nfs_infsa |\
  $fasta_split -out $split_dir -n 1000

  for g in $genomes; do

    if [ $f == $g ]; then continue; fi

    echo "Starting search from->to group: $f -> $g"

    outtag=$out_dir/${f}_${g}.blast
    if [ -e $outtag.out.gz ]; then
      echo "$outtag.out.gz exists. skipping... "
      continue
    fi

    nfs_blastdb=$blastdb_dir/$g/protein.fa
    blastdb=$tmp_dir/$g/protein.fa

    # get blast database onto local disk
    mkdir -p $tmp_dir/$g
    cp ${nfs_blastdb}.p* $tmp_dir/$g

    for splitfile in $split_dir/*.fsa; do

      splitbase=$(basename $splitfile .fsa)

      splitouttag=$outtag.$splitbase

      if [ -e $splitouttag.out ]; then
        echo "$splitouttag.out exists. skipping... "
        continue
      fi

      nice $blast \
        -p blastp \
        -i $splitfile \
        -d $blastdb \
        -a 2 \
        -e 1 -v 5 -b 0 -I T -J T \
        -o $splitouttag.out 2>&1 |\
      tee $splitouttag.log

    done

    cat $outtag.*.out | gzip -c >| $outtag.out.gz
    rm -f $outtag.*.out

    cat $outtag.*.log >| $outtag.log
    rm -f $outtag.*.log

  done

done

