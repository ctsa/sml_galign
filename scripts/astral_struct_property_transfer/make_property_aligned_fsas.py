#!/usr/bin/env python
"""given astral-fsa alignments and astral-property alignments,
produce matches to original fsa's

usage: $0 -in fastadir -out propdir 
required:
  -astral-property-file file
  -org-genome-astral-file file   (required for ss or burial)
  -orgid (tag for organism aligned to astral)
"""

# add auto-kn splitting into script
#


from bioseq_util import *



class codon_prop_t :
  "catch-all per-codon property struct"
  pass




def init_outfp(outfilename,symbol_map) :
  outfp = open(outfilename,"w")
  outfp.write("# category symbols:\n")
  for i in range(len(symbol_map)) :
    outfp.write("# %i:%s\n" % (i,symbol_map[i]))
  outfp.write("#\n")
  return outfp




def count_subs(file,target_orgid,struct_prop_map,outfilename) :
  fp = open(file,"r")

  gapchar="-"

  target_seq=""
  target_seqid=""

  # slurp all three in and separate the outgroup
  #
  fastano = 0
  seqid = orgid = ""
  for line in fp :
    if line[0] == '>' :
      fastano += 1
      
      word=line.strip().split()

      # my (ctsa) input files have the organism annotated
      # publishers (kondrashov) files are in rn-mm-hs order
      #
      if len(word) == 3 :
        orgid=word[1]
        seqid=word[2]

        if target_orgid == orgid :
          target_seqid = seqid
      else :
        raise Exception("Unexpected fasta header: %s" % line)

    else :
      if orgid == target_orgid :
        target_seq += line.strip()

  ncodon=len(target_seq)/3

  target_codon_no = 0
  
  outfp = 0
  
  for i in range(ncodon) :
    base_pos=i*3
    target_codon=target_seq[base_pos:(base_pos)+3]

    # record target codon number map to look up aligned structural properties:
    tprop = 0
    if target_codon != "---" :
      target_codon_no += 1
 
      if target_seqid != "" :
        tprop = struct_prop_map.get_prop(target_seqid,target_codon_no)

    if tprop != 0 : 
      if outfp == 0 : outfp = init_outfp(outfilename,struct_prop_map.get_symbols())
      outfp.write("%i %i\n" % (i,tprop) )




def main() :

  import os
  import sys

  indir = ""
  outdir = ""
  suffix = "cat"
  astral_property_file = ""
  org_genome_astral_file = ""
  orgid = ""

  arg = sys.argv[1:]

  i = 0
  
  while i < len(arg)  :
    if arg[i] == "-in" :
      i += 1
      indir = arg[i]
    elif arg[i] == "-out" :
      i += 1
      outdir = arg[i]
    elif arg[i] == "-astral-property-file" :
      i += 1
      astral_property_file = arg[i]
    elif arg[i] == "-org-genome-astral-file" :
      i += 1
      org_genome_astral_file = arg[i]
    elif arg[i] == "-orgid" :
      i += 1
      orgid = arg[i]
    else :
      print __doc__
      return      
    i += 1

  if indir == "" or outdir == "" or suffix == "" :
    print __doc__
    return

  # check for files required for any structural calculation
  if org_genome_astral_file == "" or astral_property_file == "" :
    print __doc__
    return

  # crank up the big fat lookup db
  import fat_pig_db
  struct_prop_map = fat_pig_db.genome_struct_property_map(org_genome_astral_file,
                                                          astral_property_file)

  symbol_map = struct_prop_map.get_symbols()

  # create output directory
  #
  if not os.path.isdir(outdir) : os.mkdir(outdir)

  # create catmap
  #
  outfilename = outdir+"/catlabel"
  outfp = open(outfilename,"w")
  for i in range(1,len(symbol_map)) :
    outfp.write("%i %s\n" % (i,symbol_map[i]))
  outfp.close()

  # over all aligned fasta files:
  #
  for infile in os.listdir(indir) :

    # file - suffix
    outfile = infile
    r = infile.rfind(".")
    if r != -1 : outfile = infile[:r]
    
    outfilename = outdir+"/"+outfile+"."+suffix
    count_subs(indir+"/"+infile,orgid,struct_prop_map,outfilename)


if __name__ == '__main__':
  main()
