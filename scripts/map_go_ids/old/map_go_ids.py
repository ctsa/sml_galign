#!/usr/bin/env python
"""
  -in tag_file
  -clustal-rna-dir dir
  -out out_file

tag_file format:
---
tag1 protein-fsa-file1 go-annotation-file1
tag2...
---
"""



def build_gi_rs_map(gmap,pfile) :

  import gzip

  pfp = gzip.open(pfile)

  for line in pfp :
    if line[0] != ">" : continue
    word = line.strip().split("|")
    gi=int(word[1])
    rs=word[3].split(".")[0]
    gmap[gi] = rs



def build_rs_go_map(rmap,gfile) :

  gfp = open(gfile)

  for line in gfp :

    word = line.strip().split('\t')

    if word[0] != "RefSeq" : continue

    rs = word[1]
    go = word[4][3:]

    if not rmap.has_key(rs) : rmap[rs] = []

    rmap[rs].append(go)



def main() :

  import os
  import sys

  tag_file = ""
  out_file = ""
  clustal_rna_dir = ""

  arg = sys.argv[1:]

  i = 0
  
  while i < len(arg)  :
    if arg[i] == "-in" :
      i += 1
      tag_file = arg[i]
    elif arg[i] == "-out" :
      i += 1
      out_file = arg[i]
    elif arg[i] == "-clustal-rna-dir" :
      i += 1
      clustal_rna_dir = arg[i]
    else :
      print __doc__
      return
    i += 1

  if tag_file == "" or clustal_rna_dir == "" :
    print __doc__
    return


  tags = {}

  # parse tag file:
  tagfp = open(tag_file)
  for line in tagfp :
    word = line.strip().split()
    tag = word[0]
    pro_fsa_file = word[1]
    go_file = word[2]
    tags[tag] = (pro_fsa_file,go_file)

  tagfp.close()


  # make gi to refseq map:
  gi_rs_map = {}
  for t in tags.keys() :
    build_gi_rs_map(gi_rs_map,tags[t][0])


  # make refseq to go map:
  rs_go_map = {}
  for t in tags.keys() :
    build_rs_go_map(rs_go_map,tags[t][1])


  # over all aligned fasta files:
  #
  outfp = sys.stdout

  for infile in os.listdir(clustal_rna_dir) :

    outtag = infile
    r = infile.rfind(".")
    if r != -1 : outtag = infile[:r]
    
    outfp.write("%s" % (outtag))

    goset = {}

    ifp = open(clustal_rna_dir+"/"+infile)
    for line in ifp :
      if line[0] != '>' : continue

      gi = int(line.strip().split()[1])
      rs = gi_rs_map[gi]
      
      if rs_go_map.has_key(rs) :
        for g in rs_go_map[gi_rs_map[gi]] : goset[g] = 1
        
    for g in goset.keys() :
      outfp.write(" %s" % g)

    outfp.write("\n")



if __name__ == '__main__':
  main()
