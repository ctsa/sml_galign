#!/usr/bin/env python
"""
for blast xml on stdin, pick best hit from the pdb given:

1) a minimum seq-identity
2) a max pdb resolution (if resolution file is supplied
3) e-val?


usage: $0 -seqid x [options]

options:
  -seqid x      ; min seqid to hit
  -coverage x   ; min hit coverage of query
  -length x     ; min hit length
  -eval x       ; max blast e-val
  -pdbres x     ; max pdb resolution
  -allow-nmr    ; consider nmr structures (as last resort)
  -resfile file ; pdb resolution index, as formated from the rcsb
"""

def wrapwrite(fp,str,wrapn) :
  "trivial print utility"

  for i,s in enumerate(str) :
    if i != 0 and i%wrapn == 0 :
      fp.write("\n")
    fp.write(s)

  fp.write("\n")  




def compress_overlaps(ranges) :
  "combine overlapping ranges from a range list"  

  r = range(len(ranges))
  r.reverse()
  for n,i in enumerate(r) :
    ff,tt = ranges[i]

    for j in r[(n+1):] :
      f,t = ranges[j]

      if f < tt and t > ff :
        if f > ff : f = ff
        if t < tt : t = tt
        ranges[j] = (f,t)
        del ranges[i]
        break
      
      






def main():

  import sys

  resfile = ""

  min_seqid = -1.
  min_coverage = -1.
  min_length = -1.
  max_eval = -1.
  max_pdbres = -1.
  max_overlap = -1.
  
  use_seqid_filter = False
  use_coverage_filter = False
  use_length_filter = False
  use_eval_filter = False
  use_res_filter = False
  use_overlap_filter = False

  allow_nmr = False
  
  arg = sys.argv[1:]

  i = 0
  while i < len(arg)  :
    if arg[i] == "-seqid" :
      i += 1
      min_seqid = float(arg[i])
      use_seqid_filter = True
    elif arg[i] == "-coverage" :
      i += 1
      min_coverage = float(arg[i])
      use_coverage_filter = True
    elif arg[i] == "-length" :
      i += 1
      min_length = int(arg[i])
      use_length_filter = True
    elif arg[i] == "-eval" :
      i += 1
      max_eval = float(arg[i])
      use_eval_filter = True
    elif arg[i] == "-overlap" :
      i += 1
      max_overlap = float(arg[i])
      use_overlap_filter = True
    elif arg[i] == "-pdbres" :
      i += 1
      max_pdbres = float(arg[i])
      use_res_filter = True
    elif arg[i] == "-allow-nmr" :
      allow_nmr = True
      use_res_filter = True
    elif arg[i] == "-resfile" :
      i += 1
      resfile = arg[i]
    else :
      print __doc__
      return
    i += 1

  if use_seqid_filter :
    if min_seqid > 1. or min_seqid < 0. :
      print __doc__
      return

  if use_coverage_filter :
    if min_coverage > 1. or min_coverage < 0. :
      print __doc__
      return

  if use_res_filter :
    if resfile == "" or max_pdbres < 0. :
      print __doc__
      return

  if use_overlap_filter :
    if max_overlap > 1. or max_overlap < 0. :
      print __doc__
      return



  # load up the structure resolution hash if needed...
  #
  resdb = {}

  if use_res_filter :
    resfp = open(resfile,"r")

    for line in resfp :
      word = line.strip().split()
      if len(word) == 3 and word[1] == ";" :
        resdb[word[0]] = float(word[2])

    resfp.close()


  infp = sys.stdin
  outfp = sys.stdout
  logfp = sys.stderr

  import blast_xml_util

  # parse blast
  #
  for query,queryid,querylen in blast_xml_util.blast_xml_querys(infp) :

    class hitprop :
      def __init__(self) :
        self.overlap = 0.
    
    hits_passed = []

    # keep track of hit overlap
    q_covered = []

    for hit in blast_xml_util.blast_xml_hits(query) :

      this_hitprop = hitprop()

      hitprop.desc = hit.hitdef
      
      # filter on structure resolution
      #
      if use_res_filter :
        # get pdb-id(s)
        #
        pdbids = []
        words = (hit.hitid+hit.hitdef).split("|")

        for i,word in enumerate(words) :
          if word == "pdb" :
            # get pdbid and chainid
            #
            pdbid = words[i+1].upper()
            chainid = words[i+2][0].upper()
            if chainid == " " or chainid == "0" :
              chainid = "_"
            
            pdbids.append((pdbid,chainid))

        pdbselected = None
        selectedres = None
        for pdbid,chainid in pdbids :
          if not resdb.has_key(pdbid) :
            sys.stderr.write("Unknown pdb: %s\n" % pdbid)
            continue

          # store the highest res pdb
          # 
          thisres = resdb[pdbid]
          if thisres < 0. and allow_nmr :
            thisres = max_pdbres
          if thisres <= max_pdbres and thisres > 0. :
            if selectedres == None or thisres < selectedres :
              pdbselected = pdbid,chainid
              selectedres = thisres
              this_hitprop.res = thisres
              this_hitprop.pdb = pdbid,chainid

        if pdbselected == None :
          # reject this hit
          #
          continue

      #filter on eval
      #
      if use_eval_filter :
        if hit.eval > max_eval :
          continue
        else :
          this_hitprop.eval = hit.eval

      #filter on seqid
      #
      if use_seqid_filter :
        seqid = float(hit.idcount) / float(hit.lencount)
        if seqid < min_seqid :
          continue
        else :
          this_hitprop.seqid = seqid


      if use_coverage_filter :
        cov = float(hit.lencount) / float(querylen)
        if cov < min_coverage :
          continue
        else :
          this_hitprop.cov = cov
      else :
	this_hitprop.cov = 1.

      if use_length_filter :
        if hit.lencount < min_length :
          continue
        else :
          this_hitprop.len = hit.lencount

      if use_overlap_filter :
        qlen = hit.q_to - hit.q_from + 1
        qlap = 0

        for f,t in q_covered :
          if f < hit.q_to and t > hit.q_from :
            if f < hit.q_from : f = hit.q_from
            if t > hit.q_to   : t = hit.q_to
            qlap += t - f + 1

        overlap = float(qlap)/float(qlen)

        if overlap > max_overlap :
          continue
        else :
          this_hitprop.overlap = overlap

          
      ######### made it!
      #
      #  write out hit alignment
      #

      # put current hit into coverage map so as to combine
      # any overlaps
      #
      q_covered.append((hit.q_from,hit.q_to))
      compress_overlaps(q_covered)
      
      
      outfp.write(">%s %i" % (queryid,hit.q_from))
      
      extra_fsa = False
      if extra_fsa :
        outfp.write("seqid: %.2f cov: %.2f len: %i eval: %.2f olap: %.2f" %\
                    (this_hitprop.seqid,
                     this_hitprop.cov,
                     this_hitprop.len,
                     this_hitprop.eval,
                     this_hitprop.overlap))
      outfp.write("\n")
      
      wrapwrite(outfp,hit.qseq,50)
      outfp.write(">%s %i\n" % (hit.hitdef.split()[0],hit.h_from))
      wrapwrite(outfp,hit.hseq,50)

      hits_passed.append(this_hitprop)
      
      continue
    
    if len(hits_passed) > 0 :
      for hp in hits_passed :
        logfp.write("%s seqid: %.2f cov: %.2f len: %i eval: %.2f olap: %.2f\n" %\
                    (queryid,hp.seqid,hp.cov,hp.len,hp.eval,hp.overlap))
    else :
      logfp.write("%s miss\n" % (queryid))


if __name__ == '__main__':
  main()

