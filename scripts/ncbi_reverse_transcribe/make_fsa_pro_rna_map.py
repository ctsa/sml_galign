#!/usr/bin/env python
"""
reverse transcriptor

approximatly total brute force solution -- scan a set of ncbi mRNA
.gbk files for relevent information and store it *ALL*, then read and
convert one folder of fastas into another folder of fastas

usage: -indir -outdir -gbk-list [gbk files] -tag-list [org id files]
"""


import sys

import os

from rt_fsa_file import rt_fsa_file




def loaddb(cdsseq,org_id,tag,db) :
  if cdsseq != "" :
    db[tag+"_"+org_id] = (org_id,cdsseq)




def read_fsa_db(db,fp,org_id) :
  """
  read from the fsa-formated stream fp, and add cds information to db
  """

  cdsseq=""
  tag=""
  for line in fp :
    if line[0] == '>' :
      loaddb(cdsseq,org_id,tag,db)

      tag = line[1:].strip().split()[0]
      tag=tag.replace("ORFN:","ORFP_")
      cdsseq = ""
    else :
      cdsseq += line.strip()
 
  loaddb(cdsseq,org_id,tag,db)




def get_gnum(line) :
  """
  get ncbi number
  """

  words = line.strip().split("|")

  if words[0][-2:] == "gi" :
    return words[1]
  else :
    #raise Exception("unexpected blast output")
    return line[1:].strip().split()[0]





def main() :

  indir = outdir = fsalist_file = ""

  arg = sys.argv[1:]

  i = 0
  while i < len(arg) :
    if arg[i] == "-indir" :
      i += 1
      indir = arg[i]
    elif arg[i] == "-outdir" :
      i += 1
      outdir = arg[i]
    elif arg[i] == "-fsa-list" :
      i += 1
      fsalist_file = arg[i]
    elif arg[i] == "-tag-list" :
      i += 1
      taglist_file = arg[i]
    else :
      print __doc__
      return
    i += 1

  # check param:
  #
  if indir == "" or outdir == "" or fsalist_file == "" :
    print __doc__
    return

  # suck in lists:
  fsalist = []
  for line in open(fsalist_file) :
    fsalist.append(line.strip())

  taglist = []
  for line in open(taglist_file) :
    taglist.append(line.strip())

  if len(taglist) != len(fsalist) :
    print "ERROR:: different size for tag and fsa lists!"
    sys.exit(1)


  # make gigantic fat-piggy database of rna:
  #
  fsa_db = {}

  import gzip

  for i in range(len(fsalist)) :
    ffile = fsalist[i]
    org_id = taglist[i]
    sys.stderr.write("starting fsa db file: %s\n" % ffile)
    ffp = gzip.open(ffile)
    read_fsa_db(fsa_db,ffp,org_id) 
    sys.stderr.write("finished fsa db file: %s\n" % ffile)
    ffp.close()


  # start reading through fasta files in
  # indir
  for file in os.listdir(indir) :
    if file[-4:] == ".fsa" :
      rt_fsa_file(indir+"/"+file,outdir,fsa_db,get_gnum)





if __name__ == '__main__':
    main()
