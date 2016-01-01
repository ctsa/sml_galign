#!/usr/bin/env python
"""
parse blast output for the top hit for each query -- regardless of e-val
"""



def get_gnum(line) :
  """get ncbi gi number"""

  words = line.strip().split("|")

  if words[0][-2:] == "gi" :
    return words[1]
  else :
    # if no ginum, get the $1 tag instead
    words = line.strip().split()
    return words[0]


#raise Exception("unexpected blast output")





import sys

infp = sys.stdin
outfp = sys.stdout

bestline = 0


for line in infp :

  # get the current query
  #
  if len(line) >= 6 and line[0:6] == "Query=" :
    query_gnum = get_gnum(line[6:])
      
  # now get the best hit, if there is one  
  #
  if line.find("Sequences producing significant alignments:") != -1 :
    bestline = 3
    
  if bestline == 1 :
    besthit_gnum = get_gnum(line)

    # add the e-val, just for kicks
    besthit_eval = line.strip().split()[-1]

    outfp.write("%s %s %s\n" % (query_gnum,besthit_gnum, besthit_eval))

  if bestline != 0 : bestline -= 1
