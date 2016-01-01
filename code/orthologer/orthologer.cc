// -*- mode: c++; indent-tabs-mode: nil; -*-
// 


// cheap cog finder -- finds orthologous clusters given the set of
// best hit edges between organisms for a set of genes (a la Koonin).
// 
// 1) looks for mutual best hit triangles
// 2) joins triangles with a common edge into larger graphs (COGS)
// 3) returns all COGS with at least three orgs
//
// current implementation is not designed to scale well, if you want
// to work with 10+ organisms, start rewriting: the triangle joiner is
// a great place to start
//
//
// file on stdin expected in format:
//   org1_label seq1_id org2_label seq2_id [ e-val ]
//   NOTE: all seqid must be unique within and between orgs
//
//
// output clusters:
//   for 1..n: seqn_id
//


#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <exception>
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

using namespace std;



template <typename T>
struct sorting_triplet {

  sorting_triplet(const T& t0,
                  const T& t1,
                  const T& t2){
    t[0] = t0;
    t[1] = t1;
    t[2] = t2;
    if(t[0]>t[1]) swap(t[0],t[1]);
    if(t[1]>t[2]) swap(t[1],t[2]);
    if(t[0]>t[1]) swap(t[0],t[1]);
  }

  bool 
  operator<(const sorting_triplet& rhs) const {
    for(unsigned i(0);i<3;++i) if(t[i] != rhs.t[i]) return t[i] < rhs.t[i];
    return false;
  }

  T t[3];
};





// id string to id num
//
struct name_lup{

  name_lup() : idmap(), strmap() {}

  unsigned getid(const string& s){
    const map<string,unsigned>::const_iterator a(idmap.find(s));
    if(a==idmap.end()){ 
      const unsigned id(strmap.size());
      idmap[s]=id;
      strmap.push_back(s);
      return id;
    }
    return a->second;
  }

  const string& getstr(unsigned i) const{
    if(i>=strmap.size()) throw exception();
    return strmap[i];
  }

private:
  map<string,unsigned> idmap;
  vector<string> strmap;
};



struct seq_vertex;
typedef sorting_triplet<seq_vertex*> triangle;



struct seq_vertex {

  seq_vertex(unsigned char _orgid = 0,
             unsigned _seqid = 0)
    : orgid(_orgid), seqid(_seqid), outedge(0), nedge(0) {}

  ~seq_vertex() { if(outedge) free(outedge); }

  void 
  add_edge(seq_vertex* to) {
    nedge++;
    outedge = (seq_vertex**) realloc((void*)outedge,nedge*sizeof(seq_vertex*));
    outedge[nedge-1] = to;
  }

  bool
  is_connected_to(const seq_vertex* s) const {
    for(unsigned i(0);i<nedge;++i){
      if(outedge[i]==s) return true;
    }
    return false;
  }

  void
  get_triangles(set<triangle>& s) {
    for(unsigned i(0);i<nedge;++i){
      if(! outedge[i]->is_connected_to(this)) continue;
      s.insert(triangle(this,outedge[i],0));

      for(unsigned j(i+1);j<nedge;++j){
        if(outedge[j]->is_connected_to(this) &&
           outedge[j]->is_connected_to(outedge[i]) &&
           outedge[i]->is_connected_to(outedge[j]) ) {
          s.insert(triangle(this,outedge[i],outedge[j]));
        }
      }
    }
  }

  unsigned get_seqid() const { return seqid; }

  //  unsigned char get_orgid() const { return orgid; }

private:
  unsigned char orgid;
  unsigned seqid;
  seq_vertex** outedge;
  unsigned char nedge;
};




static
bool
triangle_share_edge(const triangle& a,
                    const triangle& b){

  bool is_match(false);
  for(int i(0);i<3;++i){
    if(i==2 && ! is_match) break;
    for(int j(0);j<3;++j){
      if(a.t[i] && a.t[i]==b.t[j]) {
        if(is_match) return true;
        is_match=true;
        break;
      }
    }
  }
  return false;
}





static
void
get_cogs_from_triangles(const set<triangle>& tri_set,
                        vector<vector<triangle> >& cog_set){

  list<triangle> work_set;
  copy(tri_set.begin(),tri_set.end(),back_inserter(work_set));

  list<triangle>::iterator i,j,i_end;

  i=work_set.begin();
  i_end=work_set.end();
  for(;i!=i_end;++i){
    vector<triangle> cog;
    cog.push_back(*i);
    //    std::cerr << "new cog: " << i->t[0] << " " << i->t[1] << " " << i->t[2] << "\n";
    j=i;
    ++j;
    while(j!=i_end){
      bool is_erase(false);
      for(unsigned c(0);c<cog.size();++c){
        if(triangle_share_edge(cog[c],*j)){
          cog.push_back(*j);
          //std::cerr << "cog: " << j->t[0] << " " << j->t[1] << " " << j->t[2] << "\n";
          is_erase=true;
          break;
        }
      }
      if(is_erase) {
        work_set.erase(j);
        j=i;
      }
      ++j;
    }
    cog_set.push_back(cog);
  }
}




struct ortho_graph {
  typedef map<unsigned,seq_vertex*> vlup_t;


  ~ortho_graph() {
    vlup_t::iterator i=_vlup.begin(),i_end=_vlup.end();
    for(;i!=i_end;++i) delete i->second;
  }

  void
  add_blast_edge(unsigned char orgid1, unsigned seqid1,
                 unsigned char orgid2, unsigned seqid2){

    seq_vertex* v0(get_vertex(orgid1,seqid1));
    seq_vertex* v1(get_vertex(orgid2,seqid2));

    v0->add_edge(v1);
  }

  void
  get_cogs(vector<vector<unsigned> >& cog_seqid_set){
    set<triangle> triangle_set;
    
    vlup_t::const_iterator i=_vlup.begin(),i_end=_vlup.end();
    for(;i!=i_end;++i) i->second->get_triangles(triangle_set);

    vector<vector<triangle> > cog_set;

    get_cogs_from_triangles(triangle_set,cog_set);

    for(unsigned j(0);j<cog_set.size();++j){
      set<unsigned> seqid_set;
      for(unsigned k(0);k<cog_set[j].size();++k){
        for(unsigned t(0);t<3;++t){
	  if(cog_set[j][k].t[t]){
            seqid_set.insert(cog_set[j][k].t[t]->get_seqid());
          }
	}
      }
      cog_seqid_set.push_back(vector<unsigned>());
      copy(seqid_set.begin(),seqid_set.end(),back_inserter(cog_seqid_set.back()));
    }
  }

private:

  seq_vertex*
  get_vertex(unsigned char orgid,
             unsigned seqid){

    const vlup_t::const_iterator a(_vlup.find(seqid));
    if( a == _vlup.end() ){
      seq_vertex* v = new seq_vertex(orgid,seqid);
      _vlup[seqid] = v;
      return v;
    }
    return a->second;
  }

  vlup_t _vlup; // map to find vertices by seqid
};





int main() {

  // setup io
  std::istream& blastedge_infp(std::cin);
  std::ostream& outfp(std::cout);

  // graph:
  ortho_graph og;

  // organism name ids:
  name_lup olup;

  // seq_label ids:
  name_lup slup;

  // read best blast edges, add to graph:
  char linebuffer[256];
  char org1label[128];
  char seq1label[128];
  char org2label[128];
  char seq2label[128];
  while( blastedge_infp.getline(linebuffer,256,'\n') ){

    int rval = sscanf(linebuffer,"%s %s %s %s",org1label,seq1label,org2label,seq2label);
    if(rval < 4){
      std::cerr << "incorrect input format\n";
      abort();
    }
    const unsigned oid1(olup.getid(org1label));
    const unsigned oid2(olup.getid(org2label));

    const unsigned sid1(slup.getid(seq1label));
    const unsigned sid2(slup.getid(seq2label));
    og.add_blast_edge(oid1,sid1,oid2,sid2);
  }

  vector<vector<unsigned> > cog_seqid_set;
  og.get_cogs(cog_seqid_set);

  for(unsigned i(0);i<cog_seqid_set.size();++i){
    for(unsigned j(0);j<cog_seqid_set[i].size();++j){
      if(j) outfp << " ";
      outfp << slup.getstr(cog_seqid_set[i][j]);
    }
    outfp << "\n";
  }
}
