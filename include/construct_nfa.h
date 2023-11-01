#pragma once

#include <stack>

#include "robin_hood.h"
#include "lemon/smart_graph.h"
#include "lemon/maps.h"
#include "lemon/connectivity.h"


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
using amap_t = robin_hood::unordered_map<int, std::pair<const node_t*, const node_t*>>;
using node_pair_t = std::pair<node_t, node_t>;
using nfa_stack_t = std::stack<node_pair_t>;
using wmap_t = SerializingWriteMap<nfa_t>;


enum
{
    Match = 256,
    Ghost = 257,
    SplitU = 258, // Union Split
    SplitP = 259, // + Split
    SplitK = 260, // Kleene Split
};

void update_arc_map(nfa_t &NFA, lmap_t &node_map, amap_t &arc_map, const node_t &source, const node_t &target);

void print_kgraph_arcs(const nfa_t &NFA);

void export_nfa_img(nfa_t &nfa, std::string &title);

void copy_subgraph(node_pair_t &node_pair);

std::vector<int> run_top_sort(nfa_t &NFA);

void default_procedure(nfa_t &nfa, nfa_stack_t &stack, lmap_t &node_map, const int &symbol);

void concat_procedure(nfa_t &nfa, lmap_t &node_map, nfa_stack_t &stack, amap_t &arc_map);

void union_procedure(nfa_t &nfa, nfa_stack_t &stack, lmap_t &node_map, amap_t &arc_map);

void optional_procedure(nfa_t &nfa, nfa_stack_t &stack, lmap_t &node_map, amap_t &arc_map);

void kleene_procedure(nfa_t &nfa, nfa_stack_t &stack, lmap_t &node_map, const uint8_t &k, amap_t &arc_map);

void plus_procedure(nfa_t &nfa, nfa_stack_t &stack, lmap_t &node_map, const uint8_t &k, amap_t &arc_map);

void construct_kgraph(const std::string &postfix, nfa_t &nfa, lmap_t &node_map, amap_t &arc_map, const uint8_t &k);
