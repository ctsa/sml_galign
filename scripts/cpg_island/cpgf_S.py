#!/usr/bin/env python
"""
cheap cpg-island finder
usage: $0 [ -S val -only5 ] < fsa > codon_cat_labels
"""

# fixed segment scoring constants:
#
p_cpg=17
p_else=-1
S=50
is_only5=False
#k=50
k=(-p_cpg/p_else)


import sys

arg = sys.argv[1:]

i = 0
while i < len(arg)  :
  if arg[i] == "-S" :
    i += 1
    S = float(arg[i])
  elif arg[i] == "-only5" :
    is_only5 = True 
  else :
    print __doc__
    sys.exit()
  i += 1


infp=sys.stdin
outfp=sys.stdout


seq = ""
for line in infp :
  if line[0] == '>' : continue
  seq += line.strip()



ls=len(seq)
gapless_index = []
for j in range(ls) :
  if seq[j] != '-' : gapless_index.append(j)


lg=len(gapless_index)

cover=0
cum = max = 0
start = end = j = 1

while j < lg :
  if seq[gapless_index[j-1]] == 'C' and seq[gapless_index[j]] == 'G' :
    cum += p_cpg
  else :
    cum += p_else

  if cum >= max :
    max = cum
    end = j

  if cum < 0 or j == lg-1 :
    if max >= S :
      if not is_only5 or (start-1) < k:
        cover += end-start+1
        # actual segment is from start-1 to end.
        # but, output categories for codons:
        cstart=(gapless_index[start-1])/3
        if is_only5 : cstart=0

        cend=(gapless_index[end])/3
        for j2 in range(cstart,cend+1) :
          outfp.write("%i 1\n" % (j2))
        j = end
    max = 0
    cum = 0
    start = j+1
    end = j+1
  j += 1



#sys.stderr.write("ungapped len/cov: %i %i\n" % (lg,cover))


