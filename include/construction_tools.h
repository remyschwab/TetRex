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
using node_pair_t = std::tuple<node_t, node_t, size_t>;
using nfa_stack_t = std::stack<node_pair_t>;
using wmap_t = SerializingWriteMap<nfa_t>;
// using catsites_t = std::vector<arc_t>;
using catsites_t = std::vector<std::tuple<arc_t, node_t>>;

enum
{
    Match = 256,
    Ghost = 257,
    Split = 258,
    Gap = 259
};


using buffer_t = std::stack<int>;

void print_graph(nfa_t &NFA, lmap_t &nmap, const catsites_t& cats);

std::string generate_kmer_seq(uint64_t &kmer, uint8_t &k);

arc_t update_arc_map(nfa_t &NFA, lmap_t &node_map, amap_t &arc_map, const node_t &source, const node_t &target);

void print_node_pointers(const amap_t &arc_map, nfa_t &nfa);

void print_kgraph_arcs(const nfa_t &NFA);

void print_node_ids(nfa_t &NFA, lmap_t &nmap);

void export_nfa_img(nfa_t &nfa, std::string &title);

size_t copy_subgraph(node_pair_t &subgraph, nfa_t &NFA, lmap_t &node_map, node_pair_t &subgraph_copy, amap_t &arc_map);
