#pragma once

#include <stack>
#include <algorithm>

#include "robin_hood.h"
#include "lemon/smart_graph.h"
#include "lemon/maps.h"
#include "lemon/connectivity.h"
#include "lemon/dfs.h"
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

using nfa_t = lemon::SmartDigraph;
using node_t = nfa_t::Node;
using arc_t = nfa_t::Arc;
using lmap_t = nfa_t::NodeMap<int>;
using amap_t = robin_hood::unordered_map<int, std::pair<node_t, node_t>>;
using node_pair_t = std::pair<node_t, node_t>;
using nfa_stack_t = std::stack<node_pair_t>;
using wmap_t = SerializingWriteMap<nfa_t>;

enum
{
    Match = 256,
    Ghost = 257,
    Split = 258
};


using buffer_t = std::stack<int>;

void print_graph(const nfa_t &NFA, const amap_t &arc_map, lmap_t &nmap);

std::string generate_kmer_seq(uint64_t &kmer, uint8_t &k);

void update_arc_map(nfa_t &NFA, lmap_t &node_map, amap_t &arc_map, node_t &source, node_t &target);

void print_node_pointers(const amap_t &arc_map, nfa_t &nfa);

void print_kgraph_arcs(const nfa_t &NFA);

void print_node_ids(nfa_t &NFA, lmap_t &nmap);

void export_nfa_img(nfa_t &nfa, std::string &title);

void copy_subgraph(node_pair_t &subgraph, nfa_t &NFA, lmap_t &node_map, node_pair_t &subgraph_copy, amap_t &arc_map);

std::vector<int> run_top_sort(nfa_t &NFA);
