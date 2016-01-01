#!/usr/bin/env bash
#
# cat together the raw output from the 3way genome blast output
# into a good format for the ortholog clustering prog
#
# final format is:
# org1_tag seq1_tag org2_tag seq2_tag [ eval ]
#

# local fs dependencies:
#
my_dir=$(dirname $0)
top_hit_parser=$my_dir/blast_top_hit.py

#
#
if [ ! $1 ]; then
  echo "$0 :: error: no genome dir"
  exit
fi

genome_dir=$1


seq_dir=$genome_dir/start_seq
blastout_dir=$genome_dir/blast_nway
filterout=$genome_dir/blast_nway.filtered

org_tags=$(cd $seq_dir; for f in *; do if [ -d $f ]; then echo $f; fi; done)
output_tag=".blast.out.gz"



if [ ! -e $blastout_dir ]; then
  echo "$0 :: error: no blast results in genome dir"
  exit
fi


# !!non-clobber version!!
#
if [ -e $filterout ] && [ -s $filterout ]; then
  echo "$0 :: output already exists. skipping..."
  exit
fi


>| $filterout

for f in $blastout_dir/*$output_tag; do

  # parse file name to get organism tags
  #
  basef=$(basename $f $output_tag)

  org1_tag=""
  org2_tag=""

  # parse org1,org2 from filename
  #
  for g in $org_tags; do
    frontmatch=$(expr match "$basef" "\($g\)")
    backmatch=$(expr match "$basef" ".*\($g\)$")
    if [ $g == "$frontmatch" ];   then org1_tag=$g; fi
    if [ $g == "$backmatch" ]; then org2_tag=$g; fi
  done
  
  if [ ! $org1_tag ] || [ ! $org2_tag ] || [ $org1_tag == $org2_tag ]; then
    echo $0 "can't parse org tags: " $org1_tag $org2_tag
    exit
  fi

  cat $f |\
  gzip -dc |\
  $top_hit_parser |\
  awk -v o1t=$org1_tag -v o2t=$org2_tag '{print o1t,o1t "_" $1,o2t, o2t "_" $2,$3}' >>\
  $filterout

done
