


def get_prop_fsa_db(file) :
  "read in astral<->struct property alignment database"

  db = {}
  astralid = ""
  for line in open(file) :
    if line[0] == '>' :
      if astralid != "" : db[astralid] = seq
      astralid = line.strip().split()[0][2:]
      seq = ""
    else : seq += line.strip()

  if astralid != "" : db[astralid] = seq

  return db





class genome_struct_property_map :
  """
  class to parse genome-astral / astral-property files
  and respond to queries on the property{ss,burial} by
  (genome sequence id + codon number)
  """

  # primary lockup db:
  # accepts [protein_gi][codonno]
  prop_db = {}
  symbol_db = []


  def __init__(self,genome_astral_file,astral_property_file) :

    class seqinfo :
      def __init__(self) :
        self.seq = ""

  
    def process_alignment(genome_sinfo,astral_sinfo,outdb,astral_prop_db,symbol_db) :
      """
      given genome<->astral alignment and astral property alignment info, build 
      data into the property database 
      """

      genomecodon = genome_sinfo.startpos-1
      astralcodon = astral_sinfo.startpos-1

      # no prop seq is listed when prop is undetermined ('X') for entire sequence
      #
      is_prop_seq=False
      if astral_prop_db.has_key(astral_sinfo.idx) :
	is_prop_seq=True
        prop_seq = astral_prop_db[astral_sinfo.idx]

      if not outdb.has_key(genome_sinfo.pro_gi) :
        outdb[genome_sinfo.pro_gi] = {}

      for i,gs in enumerate(genome_sinfo.seq) :
        aseq = astral_sinfo.seq[i]
        if gs != "-" : genomecodon += 1
        if aseq != "-" : astralcodon += 1

        if gs == "-" or aseq == "-" : continue

        # don't double-write on overlaps, unless you can erase an 0:
        if outdb[genome_sinfo.pro_gi].has_key(genomecodon) and \
           outdb[genome_sinfo.pro_gi][genomecodon] != 0 :
          continue

        if is_prop_seq : symbol = prop_seq[astralcodon-astral_sinfo.startpos]
	else           : symbol = 'X'

        if not symbol in symbol_db : symbol_db.append(symbol)
        
        outdb[genome_sinfo.pro_gi][genomecodon] = symbol_db.index(symbol)




    self.symbol_db.append('X')

    tmp_prop_db = get_prop_fsa_db(astral_property_file) 


    import copy

    sinfo = seqinfo()
    genome_sinfo = None
    
    fastacount = 0
    is_genome = False


    # flip between genome sequences and their aligned astral
    # counterparts, use is_genome to track which one we're on
    #
    for line in open(genome_astral_file) :
      if line[0] == '>' :
        if fastacount != 0 :
          if is_genome :
            genome_sinfo = copy.deepcopy(sinfo)
          else :
            process_alignment(genome_sinfo,sinfo,self.prop_db,tmp_prop_db,self.symbol_db)
     
        fastacount += 1
        is_genome = bool(fastacount%2)
        
        words = line.strip().split()
        sinfo.idx = words[0][1:]
        sinfo.startpos = int(words[1])
        sinfo.seq = ""

	if is_genome :
	  if sinfo.idx[:3] != "gi|" and sinfo.idx[:4] != "lcl|" :
            raise Exception("GEN_STRUCT_PROP_MAP:: Unexpected fasta format in genome-astral alignment")

        if is_genome :
          sinfo.pro_gi = sinfo.idx.split("|")[1]
        else :
          # astral id's shortened by 1             
          sinfo.idx = sinfo.idx[1:]
        
      else :
        sinfo.seq += line.strip()

    # process the last alignment:
    if fastacount != 0 and not is_genome :
      process_alignment(genome_sinfo,sinfo,self.prop_db,tmp_prop_db,self.symbol_db)



  def get_symbols(self) :
    "list symbols corresponding to category numbers"

    return self.symbol_db


  def get_prop(self,gi,codon) :
    "get category corresponding to protein gi num and (1-indexed) codon number"
    
    if self.prop_db.has_key(gi) :
      if self.prop_db[gi].has_key(codon) :
        return self.prop_db[gi][codon]
    return 0
