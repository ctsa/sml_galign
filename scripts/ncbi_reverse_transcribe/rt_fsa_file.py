
import sys
import os
import os.path

from bioseq_util import *


def rt_fsa_file(fsafile,outdir,db,get_gnum) :
  """
  make the reverse transcript fasta file from the input aa
  fasta file according to database of mRNA cds's
  """

  gapchar = "-"

  fsafp = open(fsafile)

  outfile = outdir+"/"+os.path.basename(fsafile)
  outfp = open(outfile,"w")
  
  for line in fsafp :
    if line[0] == '>' :
      pro_gi = get_gnum(line)
      if not db.has_key(pro_gi) :
	sys.stderr.write("WARNING:: protein id missing in mRNA database: %s\n" % (pro_gi))
        outfp.close()
        os.remove(outfile)
        #raise Exception("protein id missing in mRNA database")
        return

      org_id,cdsseq = db[pro_gi]
      codon_no=0

      outfp.write("> %s %s\n" % (org_id,pro_gi))

    else :
      for aa in line.strip() :
        if aa == gapchar :
          outfp.write(gapchar*3)
        else :
          base_no=codon_no*3
          codon=cdsseq[base_no:base_no+3]
          if len(codon) != 3 :
            raise Exception("sequence not divisible to codons.  pro id: %s, seq: %s\n" % (pro_gi,cdsseq))
	    #sys.stderr.write("WARNING::sequence not divisible to codons.  pro id: %s seq_size: %i\n" % (pro_gi,len(cdsseq)))
            #outfp.close()
      	    #os.remove(outfile)
	    #return
          elif aa != codon_trans(codon) :
            #raise Exception("bad codon in mRNA for codon/aa: %s/%s" % (codon,aa))
            sys.stderr.write("WARNING::bad codon in mRNA for codon/aa: %s/%s, pro id: %s\n" % (codon,aa,pro_gi))
            outfp.close()
      	    os.remove(outfile)
      	    return

          outfp.write(codon)
          codon_no += 1

      outfp.write("\n")
