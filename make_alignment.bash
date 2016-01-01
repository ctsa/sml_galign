#!/usr/bin/env bash

#
# usage:
# -g genome_dir
# -a astral_target_org
# -d disorder org [ disordered protein regions ]
# -c cpg island org
# -t transmem org [ phobius transmem helix prediction ]
# -s "org newick tree" [ dN/dS calculation ]
# -y [ yeast genome format ]
#


genome_dir=""
astral_target_org=""
is_yeast_format=0
disorder_target_org=""
cpgi_target_org=""
transmem_org=""
dnds_tree=""


while getopts ":g:a:d:c:t:s:yw" Option
do
case $Option in
  g ) genome_dir=$OPTARG;;
  a ) astral_target_org=$OPTARG;;
  d ) disorder_target_org=$OPTARG;;
  c ) cpgi_target_org=$OPTARG;;
  t ) transmem_org=$OPTARG;;
  s ) dnds_tree="$OPTARG";;
  y ) is_yeast_format=1;;
  * ) echo "Unimplemented option chosen."
      exit 1;;
esac
done

shift $(($OPTIND - 1))


if [ ! $genome_dir ]; then
  echo "no genome dir"
  exit
fi



my_dir=$(dirname $0)

scripts_dir=$my_dir/scripts
blast_nway=$scripts_dir/blast_nway/make_blast_nway.bash
filter_blast_nway=$scripts_dir/filter_blast_nway/blast_nway_to_ortholog_graph.bash
make_ortho_groups=$scripts_dir/make_ortho_graphs/make_orthologs.bash
make_clustal=$scripts_dir/make_clustalw_align/make_clustal.bash
make_rna=$scripts_dir/ncbi_reverse_transcribe/make_rna.bash
astral_genome_alignment=$scripts_dir/astral_genome_alignment/get_genome_astral_align.bash
astral_prop_transfer=$scripts_dir/astral_struct_property_transfer/astral_struct_property_transfer.bash
make_dnds=$scripts_dir/dnds/make_dnds.bash
map_go=$scripts_dir/map_go_ids/map_go_ids.bash
get_disorder=$scripts_dir/disorder/get_disorder.bash
get_cpgi=$scripts_dir/cpg_island/get_start_cpg.bash
make_transmem=$scripts_dir/transmem/make_phobius_cats.bash
gc_win=$scripts_dir/gc_content_window/get_gc_content.bash
gc_gene=$scripts_dir/gc_content_gene/get_gc_content.bash
glen=$scripts_dir/gene_length/get_gene_length_cats.bash



echo
echo "########## step 1 -- prepare all against all blast alignments #############"
echo

$blast_nway $genome_dir 


echo
echo "########## step 2 -- format blast results into a simple column flatfile ############"
echo

$filter_blast_nway $genome_dir


echo
echo "########## step 3 -- make orthologous groups from blast results ############"
echo

$make_ortho_groups $genome_dir


echo
echo "########## step 4 -- do clustal alignments for orthologous groups ############"
echo

$make_clustal $genome_dir


echo
echo "########## step 5 -- \"reverse transcribe\" the clustal alignment via db lookup ############"
echo

$make_rna $genome_dir $is_yeast_format


if [ "$cpgi_target_org" != "" ]; then
  echo
  echo "########## optional step 1 -- calculate cpg islands ############"
  echo

  $get_cpgi $genome_dir "$cpgi_target_org"
fi


if [ "$astral_target_org" != "" ]; then

  echo
  echo "########## optional step 2, part 1 -- align representative genome to structure ############"
  echo

  min_pdb_seqid=30
  max_overlap=40

  $astral_genome_alignment $genome_dir "$astral_target_org" $min_pdb_seqid $max_overlap


  echo
  echo "########## optional step 2, part 2 -- label genome to structure alignments with structure properties ############"
  echo

  $astral_prop_transfer $genome_dir "$astral_target_org" $min_pdb_seqid $max_overlap

fi


if [ "$disorder_target_org" != "" ]; then

  echo
  echo "########## optional step 3 -- find disordered regions ############"
  echo

  $get_disorder $genome_dir "$disorder_target_org"
fi


if [ "$dnds_tree" != "" ]; then
  echo
  echo "########## optional step 4 -- calculate dnds ############"
  echo

  $make_dnds $genome_dir "$dnds_tree"
fi


if [ "$transmem_org" != "" ]; then
  echo
  echo "########## optional step 5 -- find transmem helices ############"
  echo

  $make_transmem $genome_dir "$transmem_org"
fi


echo
echo "########## optional step 6a -- get windowed G+C% regions ############"
echo

$gc_win $genome_dir


echo
echo "########## optional step 6b -- get gene G+C% ############"
echo

$gc_gene $genome_dir

echo
echo "########## optional step 6c -- get gene_3c G+C% ############"
echo

$gc_gene $genome_dir -3c

echo
echo "########## optional step 6d -- get gene_3cf G+C% ############"
echo

$gc_gene $genome_dir -3cf

echo
echo "########## optional step 6e -- get gene length ############"
echo

$glen $genome_dir



#echo
#echo "########### optional step -- assign go categories ###########"
#echo

#$map_go $genome_dir


