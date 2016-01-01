#!/usr/bin/env python
"""
make fasta tuples for ortholog tuple sets specified by gi number --

this is a brute force approach of sucking all proteomes in ram and
spitting them out as individually numbered files, with enough
organisms something more efficient will be required

usage:

$0 -fastadb-list fasta.db.filelist -tag-list organism.tag.list -out fasta_prefix < tuple_list
"""


def get_gnum(line) :
  """
  get ncbi gi number from ncbi fasta header
  """

  words = line.strip().split("|")

  if words[0][-2:] == "gi" :
    return words[1]
  else :
    words = line.strip().split()
    return words[0]

    # raise Exception("unexpected fasta header")






def add_fasta_to_db(fastafile,org_id,fastahash) :
  """
  add file of fasta sequences to global fasta db
  """

  import gzip

  dbfp = gzip.open(fastafile,"r")

  gnum = thisfasta = thisheader = ""

  for line in dbfp :
    if line[0] == '>' :
      if gnum != "" :
        fastahash[org_id+"_"+gnum] = (thisheader,org_id,thisfasta)
      thisfasta = ""
      gnum = get_gnum(line[1:])
      thisheader = line.strip().split()[0]

    else :
      thisfasta += line

  if gnum != "" :
    fastahash[org_id+"_"+gnum] = (thisheader,org_id,thisfasta)

  dbfp.close()






def main():

  import sys

  fastadb_list = tag_list = outtag = ""

  arg = sys.argv[1:]

  i = 0
  while i < len(arg)  :
    if arg[i] == "-fastadb-list" :
      i = i + 1
      fastadb_list = arg[i]
    elif arg[i] == "-tag-list" :
      i = i + 1
      tag_list = arg[i]
    elif arg[i] == "-out" :
      i = i + 1
      outtag = arg[i]
    else :
      print __doc__
      return
    i = i + 1

  if fastadb_list == "" or tag_list == "" or outtag == "" :
    print __doc__
    sys.exit()


  infp = sys.stdin

  fastahash = {}

  fastadb = []

  for line in open(fastadb_list) :
    fastadb.append(line.strip())

  tagdb = []

  for line in open(tag_list) :
    tagdb.append(line.strip())

  if len(tagdb) != len(fastadb) :
    print "fastadb and tag lists are different lengths"
    exit(1)

  for i in range(len(fastadb)) :
    add_fasta_to_db(fastadb[i],tagdb[i],fastahash)

  outputcount = 0

  for line in infp :
    outputcount += 1
    outfile = "%s%05i.fsa" % (outtag,outputcount)
    outfp = open(outfile,"w")

    word = line.strip().split()

    if len(word) < 2 :
      raise Exception("unexpected tuple format; less than 2 seqids: %s" % line)

    for i in range(len(word)) :
      if not fastahash.has_key(word[i]) :
        raise Exception("seqid does not exist in database: %s" % word[i])

    for i in range(len(word)):
      header,org_id,seq = fastahash[word[i]]

      # org_id is fused to the seq_id in word1 of the fasta header to
      # guarantee the uniqueness of all sequence labels in each fasta
      # alignment. So far, this has only been a problem for the fungi
      # alignments, b/c of their short seq_id's. For ncbi sourced
      # sequences, the seq_id's alone (with gi number embedded) should
      # be unique w/o the added org_id.
      #
      seq_id=header[1:].strip()
      outfp.write(">%s_%s %s\n" % (seq_id,org_id,org_id) )
      outfp.write(seq)

    outfp.close()



if __name__ == '__main__':
    main()
