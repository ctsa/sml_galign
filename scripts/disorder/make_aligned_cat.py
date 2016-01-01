#!/usr/bin/env python

import sys

if len(sys.argv) < 3 :
  print "usage: ",sys.argv[0]," fsa_file disorder_file"


outfp=sys.stdout

fsafile=sys.argv[1]
disfile=sys.argv[2]


seq=""


fsafp=open(fsafile)
for line in fsafp :
  if line[0] == ">" : continue
  seq += line.strip()
fsafp.close()


dis_seq = []

disfp=open(disfile)
for line in disfp :
  if line[0] == "R" : continue
  word=line.strip().split()
  dis_seq.append( (word[1],word[2]) )
disfp.close()
  
i=0
disi=0
for s in seq :
  if s != "-" and s != "X" :
    if s != dis_seq[disi][0] :
      sys.stderr.write("AA mismatch s,di,disi,seq: %s %s %i %s\n" % (s,dis_seq[disi][0],disi,seq))
      fsafp=open(fsafile)
      for line in fsafp :
	sys.stderr.write(line)
      sys.exit(1)

    outfp.write("%i %s\n" % (i,dis_seq[disi][1]))
    disi+=1

  i+=1


