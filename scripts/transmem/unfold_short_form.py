#!/usr/bin/env python
"""
parse phobius prediction line (short form) and corresponding
aligned fsa sequence to make codeq cat file
"""

import re
import sys


if len(sys.argv) < 2 :
  print "usage: ",sys.argv[0]," fsa_file"
  sys.exit(1)


outfp=sys.stdout

fsafile=sys.argv[1]



#outfp.write(head)


seq=""


fsafp=open(fsafile)
for line in fsafp :
  if line[0] == ">" : continue
  seq += line.strip()
fsafp.close()


line=sys.stdin.readline()

word=line.strip().split()

tcount=0
if len(word) == 4 : tcount=int(word[1])

is_transmem=True
prop_seq = []

if tcount == 0 :
  is_transmem=False
else :
  nosp=word[3].split('/')
  if len(nosp) == 2 : nosp=nosp[1]
  else              : nosp=nosp[0]

  tmseq=re.split('[io]',nosp)

  for t in tmseq :
    if t.find('-') != -1 :
      begin,end=t.split('-')
      begin=int(begin)-1
      end=int(end)
      for i in range(begin-len(prop_seq)) :
	prop_seq.append(1)
      for i in range(end-begin) :
	prop_seq.append(2)


i=0
propi=0
for s in seq :
  if s == "-" or s == "X" :
    if is_transmem :
      pass
#      outfp.write("%i 0\n" % i)
    else :
      outfp.write("%i 1\n" % i)
  else :
    if not is_transmem or propi >= len(prop_seq) :
      outfp.write("%i %s\n" % (i,1))
    else :
      outfp.write("%i %s\n" % (i,prop_seq[propi]))
      propi+=1

  i+=1



if propi < len(prop_seq) :
  sys.stderr.write("prop_seq len mismatch propi,len,line: %i %i %s\n" % (propi,len(prop_seq),line) )
  fsafp=open(fsafile)
  for line in fsafp :
    sys.stderr.write(line)
    sys.exit(1)


