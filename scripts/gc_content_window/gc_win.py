#!/usr/bin/env python
"""
cheap windowed gc-content
usage: $0 [ -c -w winsize ] < multi-seq fsa > base/codon gc%
"""




def get_gc_content(winsize,is_codon,infp,outfp) :

  seq = []
  for line in infp :
    if line[0] == '>' :
      seq.append("")
    else :
      seq[-1] += line.strip()

  ls=len(seq[0])

  i = 0
  while i<ls :
    count=gc_count=0

    start=max(i-winsize,0)
    stop=min(i+winsize,ls)
    if is_codon : stop=min(stop+3,ls)

    for j in range(start,stop) :
      for s in seq :
        if s[j] == '-' or s[j] == 'N' :
          pass
        else :
          count += 1
          if s[j] == 'C' or s[j] == 'G' :
            gc_count += 1

    gc=float(gc_count)/float(count)

    index=i
    if is_codon : index /= 3
    outfp.write("%i %f\n" % (index,gc) )

    i += 1
    if is_codon : i += 2




def main():

  import sys

  is_codon=False
  winsize=200

  arg = sys.argv[1:]

  i = 0
  while i < len(arg)  :
    if arg[i] == "-c" :
      is_codon=True
    elif arg[i] == "-w" :
      i += 1
      winsize = int(arg[i])
    else :
      print __doc__
      return
    i += 1

  infp = sys.stdin
  outfp = sys.stdout

  get_gc_content(winsize,is_codon,infp,outfp)


if __name__ == '__main__':
  main()


