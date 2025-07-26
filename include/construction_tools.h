#pragma once

#include <stack>
#include <algorithm>

#include "robin_hood.h"
#include "lemon/smart_graph.h"
#include "lemon/maps.h"
#include "lemon/connectivity.h"
#include "lemon/dfs.h"
#include "lemon/adaptors.h"
#include "utils.h"


template< typename GR >
class SerializingWriteMap
{
public:
  typedef typename GR::Node  Key;
  typedef bool  Value;
  typedef std::list<Key>  OrderedList;

  OrderedList  order;

  SerializingWriteMap(const GR&) {}

  void set(const Key& k, const Value& v)
  {
    if( v ) order.push_front(k);
  }

  Value operator[] (const Key&)
  {
    return false;
  }
};

constexpr std::pair<size_t, size_t> OPT_QUANT(0,1);



using nfa_t = lemon::SmartDigraph;
using node_t = nfa_t::Node;
using arc_t = nfa_t::Arc;
using lmap_t = nfa_t::NodeMap<int>;
using gmap_t = nfa_t::NodeMap<size_t>;
using amap_t = robin_hood::unordered_map<int, std::pair<node_t, node_t>>;
// using node_pair_t = std::pair<node_t, node_t>;
// using node_pair_t = std::tuple<node_t, node_t, size_t>;

struct Subgraph
{
  node_t start;
  node_t end;
  size_t split_run_count{};
  size_t raw_split_count{};
  uint8_t origin;
};

using nfa_stack_t = std::stack<Subgraph>;
using wmap_t = SerializingWriteMap<nfa_t>;

struct Catsite
{
  node_t cleavage_site_; // Node before a high-complexity subgraph begins
  node_t cleavage_start_; // The entry node for a high-complexity subgraph
  node_t cleavage_end_; // Exit node for a high-complexity subgraph
  arc_t arc_; // I only need this to print the graph. Maybe I can get rid of it later
  node_t downstream_; // Connection point for the cleavage site node

  size_t cleavage_site_id_;
  size_t cleavage_start_id_;
  size_t cleavage_end_id_;
  size_t downstream_id_;

  void addIDs(const nfa_t& nfa)
  {
    cleavage_site_id_ = nfa.id(cleavage_site_);
    cleavage_start_id_ = nfa.id(cleavage_start_);
    cleavage_end_id_ = nfa.id(cleavage_end_);
  }

  void complete(const nfa_t& nfa, const amap_t& arcs)
  {
    downstream_ = arcs.find(nfa.id(cleavage_end_))->second.first;
    downstream_id_ = nfa.id(downstream_);
  }

  void dumpInfo(const std::vector<int>& ranks)
  {
    seqan3::debug_stream << "FUSE: " << ranks[cleavage_site_id_] << " CLEAVE START: " << ranks[cleavage_start_id_]<< " CLEAVE END: " << ranks[cleavage_end_id_] << " FUSE END: " << ranks[downstream_id_] << std::endl; 
  }

};

using catsites_t = std::vector<Catsite>;



enum NodeTypes
{
    Match = 256,
    Ghost = 257,
    Split = 258,
    Gap = 259
};

enum Procedures
{
  Default = 0,
  Concat = 1,
  Union = 2,
  Optional = 3,
  Kleene = 4,
  Plus = 5
};


using buffer_t = std::stack<int>;

void print_graph(nfa_t &NFA, lmap_t &nmap, const catsites_t& cats, const bool& augment);

std::string generate_kmer_seq(uint64_t &kmer, uint8_t &k);

arc_t update_arc_map(nfa_t &NFA, lmap_t &node_map, amap_t &arc_map, const node_t &source, const node_t &target);

void print_node_pointers(const amap_t &arc_map, nfa_t &nfa);

void print_kgraph_arcs(const nfa_t &NFA);

void print_node_ids(nfa_t &NFA, lmap_t &nmap);

void copy_subgraph(const Subgraph &subgraph, nfa_t &NFA, lmap_t &node_map, Subgraph &subgraph_copy, amap_t &arc_map);