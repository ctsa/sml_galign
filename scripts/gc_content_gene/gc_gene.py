#!/usr/bin/env python
"""
cheap alignment G+C% 
usage: $0 [ -3c | -3cf ] < multi-seq fsa > alignment G+C% 
"""

# for use with -3cf option only...
filter_threshold=0.08

import sys

if len(sys.argv) > 2 :
  print __doc__
  sys.exit(1)


is_3c=False
is_3cf=False

if len(sys.argv) > 1 :
  if sys.argv[1] == "-3c" : 
    is_3c=True
  if sys.argv[1] == "-3cf" : 
    is_3c=True
    is_3cf=True

infp = sys.stdin
outfp = sys.stdout

seq = []
for line in infp :
  if line[0] == '>' :
    seq.append("")
  else :
    seq[-1] += line.strip()

ls=len(seq[0])

count=gc_count=0

scount=[0]*len(seq)
sgc_count=[0]*len(seq)

start=0
stop=ls

for j in range(start,stop) :
  if is_3c and j%3 != 2 : continue
  for k in range(len(seq)) :
    s=seq[k]
    if s[j] == '-' or s[j] == 'N' :
      pass
    else :
      count += 1
      scount[k] += 1
      if s[j] == 'C' or s[j] == 'G' :
        gc_count += 1
	sgc_count[k] += 1

gc=float(gc_count)/float(count)
is_skip=False
  
if is_3cf :
  sgc=[0]*len(seq)
  for k in range(len(seq)) :
    sgc[k]=float(sgc_count[k])/float(scount[k])

  for k1 in range(len(seq)) :
    for k2 in range(k1+1,len(seq)) :
      if abs(sgc[k1]-sgc[k2]) > filter_threshold :
	is_skip=True
      

if not is_skip : outfp.write("%f" % (gc) )
outfp.write("\n")

