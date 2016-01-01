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




def get_cdsseq(cds,seq) :
  """
  process cds tag to produce rna sequence
  """

  is_complement = False
  
  if len(cds) >= 10 and cds[:10] == "complement" :
    is_complement = True
    
    # strip off complement tag and continue with regular processing
    #
    cds = cds[11:-1]
  
  if len(cds) >= 4 and cds[:4] == "join" :
    cds_ranges = cds[5:-1].split(",")
  else :
    cds_ranges = [ cds ]

  regions = []

  is_complete = True
  for r in cds_ranges :
    d = r.find("..")

    start = r[:d]
    if start[0] == "<" :
      start = start[1:]
      is_complete = False
    start = int(start)

    stop = r[d+2:]
    if stop[0] == ">" :
      stop = stop[1:]
      is_complete = False
    stop = int(stop)

    regions.append( (start,stop) )

  cdsseq = ""
  for start,stop in regions :
    cdsseq += seq[start-1:stop]

  cdsseq = cdsseq.upper()

  import bioseq_util

  if is_complement :
    string_reverse = lambda s: ''.join([s[i] for i in xrange(len(s)-1, -1, -1)])
    cdsseq = string_reverse(cdsseq)
    cdsseq = bioseq_util.base_comp_str(cdsseq)

  return cdsseq,is_complete








class cds_info :
  pass




def process_cds(cds_buffer) :
  """
  parse a single cds buffer into the cds_info struct and return
  """

  cds_line_continue = False

  this_cds = cds_info()

  is_gi_found = False

  for line in cds_buffer :

    if cds_line_continue :
      this_cds.cds += line.strip().lstrip()
      
      # check to see if the cds wraps
      if this_cds.cds.count("(") != this_cds.cds.count(")") :
        cds_line_continue = True
      else :
        cds_line_continue = False

    if len(line) >= 8 and line[:8] == "     CDS" :
      word = line.strip().split()
      this_cds.cds = word[1]

      # check to see if the cds wraps
      if this_cds.cds.count("(") != this_cds.cds.count(")") :
        cds_line_continue = True

    if line.find("/protein_id=") != -1 :
      word = line.strip().split('"')
      pro_id = word[1]

    if line.find("/db_xref=\"GI:") != -1 :
      word = line.strip().split(":")[1][:-1]
      this_cds.pro_gi = word
      is_gi_found = True

    # look out for special cases that this script is unable to handle
    #
    if line.find("/exception=\"artificial frameshift\"") != -1 or \
       line.find("/exception=\"unclassified translation discrepancy\"") != -1 :
      return None

    # ..and biology we aren't interested in:
    #
    if line.find("/note=\"pseudogene\"") != -1 :
      return None

    if line.find("/transl_table=") != -1 :
      tablecode=int(line.split('=')[1])
      if tablecode != 1 and tablecode != 11: # euk & bacteria
	return None

    if line.find("/codon_start=") != -1 :
      if int(line.split('=')[1]) != 1 :
	return None

  if is_gi_found :
    return this_cds
  else :
    return None






def read_gbk_db(db,fp,org_id) :
  """
  read from the gbk-formated stream fp, and add cds information to db
  """

  cds_list = []
  suckseq = False
  seqarray = []
  cds_scan = False
  cds_buffer = []
  
  for line in fp :
    if len(line) >= 5 and line[:5] == "LOCUS" :   
      # start of gbk record
      #
      word = line.strip().split()
      # don't do anything with this anymore:
      rna_id = word[1]

    if cds_scan and line[:6] != "      " :
      cds_scan = False
      cds_stuff = process_cds(cds_buffer)
      if cds_stuff != None : cds_list.append(cds_stuff)

    if not cds_scan and len(line) >= 8 and line[:8] == "     CDS" :
      cds_scan = True
      cds_buffer = []

    if cds_scan :
      cds_buffer.append(line)

      
    if len(line) >= 2 and line[:2] == "//" :
      # end of gbk record
      #
      seq = "".join(seqarray)
      for cds in cds_list :
        if cds.pro_gi == "" : continue
        cdsseq,is_complete = get_cdsseq(cds.cds,seq)
	if not is_complete : continue
	if len(cdsseq)%3 != 0 :
	  sys.stderr.write("WARNING::cds sequence not divisible to codons.  pro id: %s seq_size: %i\n" % (cds.pro_gi,len(cdsseq)))
	  continue

        db[cds.pro_gi] = (org_id,cdsseq)

      # reset for next record
      cds_list = []
      suckseq = False
      seqarray = []

    # read in rna/genome sequence:
    #
    if suckseq :
      word = line.strip().split()
      for w in word[1:] :
        seqarray.append(w)

    if len(line) >= 6 and line[:6] == "ORIGIN" :
      suckseq = True   






def get_gnum(line) :
  """
  get ncbi number
  """

  words = line.strip().split("|")

  if words[0][-2:] == "gi" :
    return words[1]
  else :
    # shouldn't happen for current app
    raise Exception("unexpected blast output")
    return words[1]






def main() :

  indir = outdir = gbklist_file = ""

  arg = sys.argv[1:]

  i = 0
  while i < len(arg)  :
    if arg[i] == "-indir" :
      i += 1
      indir = arg[i]
    elif arg[i] == "-outdir" :
      i += 1
      outdir = arg[i]
    elif arg[i] == "-gbk-list" :
      i += 1
      gbklist_file = arg[i]
    elif arg[i] == "-tag-list" :
      i += 1
      taglist_file = arg[i]
    else :
      print __doc__
      return
    i += 1

  # check param
  #
  if indir == "" or outdir == "" or gbklist_file == "" :
    print __doc__
    return

  # suck in lists:
  gbklist = []
  for line in open(gbklist_file) :
    gbklist.append(line.strip())

  taglist = []
  for line in open(taglist_file) :
    taglist.append(line.strip())

  if len(taglist) != len(gbklist) :
    print "ERROR:: different size for tag and gbk lists!"
    sys.exit(1)


  # make gigantic fat-piggy database of rna:
  #
  gbk_db = {}

  import gzip

  for i in range(len(gbklist)) :
    gfile = gbklist[i]
    org_id = taglist[i]
    sys.stderr.write("starting gbk file: %s\n" % gfile)
    gbkfp = gzip.open(gfile)
    read_gbk_db(gbk_db,gbkfp,org_id) 
    sys.stderr.write("finished gbk file: %s\n" % gfile)
    gbkfp.close()


  # start reading through fasta files in
  # indir
  for file in os.listdir(indir) :
    if file[-4:] == ".fsa" :
      rt_fsa_file(indir+"/"+file,outdir,gbk_db,get_gnum)





if __name__ == '__main__':
    main()
