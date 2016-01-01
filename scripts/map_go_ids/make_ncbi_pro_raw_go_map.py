#!/usr/bin/env python
"""
get list of associated gene categories for each rna alignment (annotated with
ncbi protein GI numbers, from ncbi rna gbk files for each organism in the alignment)

usage: -rnadir d -out f -gbk-list [gbk files] -tag-list [org id files]
"""

import sys
import os



def get_go_from_note(n,s) :
  """
  parse a stripped version of the refseq /note section for goid info
  """
  wa = n.split(";")
  for w in wa :
    g=w.find("[goid")
    if g != -1 :
      wg=w[g+6:]
      wge=wg.find("]")
      if wge == -1 : raise Exception("unexpected wge")
      s.add(wg[:wge])



class cds_info :
  pass



def process_cds(cds_buffer) :
  """
  parse a single cds buffer into the cds_info struct and return
  """

  is_note_lines = False
  is_gi_found = False

  this_cds = cds_info()
  this_cds.goid=set()

  for line in cds_buffer :
    if is_note_lines :
      note_buffer += line.strip()
      nls += 1
      if nls > 300 : raise Exception("long note segment")
      if note_buffer[-1] == "\"" :
        get_go_from_note(note_buffer,this_cds.goid)
	is_note_lines = False

    if line.find("/note=") != -1 :
      note_buffer = line.strip()
      nls = 1
      if note_buffer[-1] == "\"" :
        get_go_from_note(note_buffer,this_cds.goid)
      else :
        is_note_lines = True

    if line.find("/db_xref=\"GI:") != -1 :
      word = line.strip().split(":")[1][:-1]
      this_cds.pro_gi = int(word)
      is_gi_found = True

  if is_gi_found : return this_cds
  else :           return None



def read_gbk_db(db,fp) :
  """
  read from the gbk-formated stream fp, and add go information to db
  """

  cds_list = []
  cds_scan = False
  cds_buffer = []
  
  for line in fp :
    if cds_scan and line[:6] != "      " :
      cds_scan = False
      cds_tmp = process_cds(cds_buffer)
      if cds_tmp != None : cds_list.append(cds_tmp)

    if not cds_scan and len(line) >= 8 and line[:8] == "     CDS" :
      cds_scan = True
      cds_buffer = []

    if cds_scan :
      cds_buffer.append(line)

      
    if len(line) >= 2 and line[:2] == "//" :
      # end of gbk record
      #
      for cds in cds_list :
        if cds.pro_gi != "" : db[cds.pro_gi] = cds.goid
      cds_list = []



def get_gi_list(file) :

  l=[]
  fp=open(file)
  for line in fp :
    if line[0] == ">" :
      w=line.strip().split()
      l.append(int(w[2]))
  return l



def main() :

  rnadir = outfile = gbklist_file = ""

  arg = sys.argv[1:]

  i = 0
  while i < len(arg)  :
    if arg[i] == "-rnadir" :
      i += 1
      rnadir = arg[i]
    elif arg[i] == "-out" :
      i += 1
      outfile = arg[i]
    elif arg[i] == "-gbk-list" :
      i += 1
      gbklist_file = arg[i]
    else :
      print __doc__
      return
    i += 1

  # check param
  #
  if rnadir == "" or outfile == "" or gbklist_file == "" :
    print __doc__
    return

  # suck in lists:
  gbklist = []
  for line in open(gbklist_file) :
    gbklist.append(line.strip())

  # make gigantic fat-piggy database of pro-gi -> go relationships:
  #
  gbk_db = {}

  import gzip
  import os.path

  for i in range(len(gbklist)) :
    gfile = gbklist[i]
    sys.stderr.write("starting gbk file: %s\n" % gfile)
    gbkfp = gzip.open(gfile)
    read_gbk_db(gbk_db,gbkfp) 
    sys.stderr.write("finished gbk file: %s\n" % gfile)
    gbkfp.close()

  outfp=open(outfile,"w")

  for file in os.listdir(rnadir) :
    if file[-4:] != ".fsa" : continue
    idend=file.find(".aligned.fsa")
    id=file[:idend]
    s=set()
    pro_gi_list=get_gi_list(os.path.join(rnadir,file))
    for p in pro_gi_list : 
      if gbk_db.has_key(p) : s.update(gbk_db[p])
    outfp.write(id)
    for i in s : outfp.write(" %s" % i)
    outfp.write("\n")

           


if __name__ == '__main__':
    main()
