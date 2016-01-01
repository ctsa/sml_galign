#!/usr/bin/env python
"""
get list of associated gene categories for each rna alignment (annotated with
ncbi protein GI numbers, from ncbi rna gbk files for each organism in the alignment)

usage: -infile f -outfile f -squish f
"""


import sys


def main() :

  infile = outfile = squish_file = ""

  arg = sys.argv[1:]

  i = 0
  while i < len(arg)  :
    if arg[i] == "-in" :
      i += 1
      infile = arg[i]
    elif arg[i] == "-out" :
      i += 1
      outfile = arg[i]
    elif arg[i] == "-squish" :
      i += 1
      squish_file = arg[i]
    else :
      print __doc__
      return
    i += 1

  # check param
  #
  if infile == "" or outfile == "" or squish_file == "" :
    print __doc__
    return

  # make gigantic fat-piggy database of go -> golabel map
  #
  go_db = {}
  go_label = []

  for line in open(squish_file) :
    if line[0] == '>' :
      label=line[1:].strip().split()[0]
      go_label.append(label)
      lid=len(go_label)-1
    else :
      id=line.strip().split()[0]
      go_db[id]=lid

  outfp=open(outfile,"w")

  for line in open(infile) :
    w=line.strip().split()
    tag=w[0]
    s=set()
    for wa in w[1:] :
      if go_db.has_key(wa) : s.add(go_db[wa])
    ss=list(s)
    ss.sort()
      
    outfp.write("%s" % tag)
    for sss in ss :
      outfp.write(" %s" % go_label[sss])
    outfp.write("\n")


if __name__ == '__main__':
    main()
