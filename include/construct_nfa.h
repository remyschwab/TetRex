#pragma once

#include <stack>

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

  SerializingWriteMap( const GR& ) {}

  void  set( const Key& k, const Value& v )
  {
    if( v ) order.push_front(k);
  }

  Value  operator[] ( const Key& )
  {
    return false;
  }
};

using nfa_t = lemon::SmartDigraph;
using node_t = lemon::SmartDigraph::Node;
using lmap_t = lemon::SmartDigraph::NodeMap<int>;
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

void export_nfa_img(nfa_t &nfa, std::string &title);

void copy_subgraph(node_pair_t &node_pair);

wmap_t run_top_sort(nfa_t &NFA, lmap_t &node_map);

void default_procedure(nfa_t &nfa, nfa_stack_t &stack, lmap_t &node_map, const int &symbol);

void concat_procedure(nfa_t &nfa, nfa_stack_t &stack);

void union_procedure(nfa_t &nfa, nfa_stack_t &stack, lmap_t &node_map);

void optional_procedure(nfa_t &nfa, nfa_stack_t &stack, lmap_t &node_map);

void kleene_procedure(nfa_t &nfa, nfa_stack_t &stack, lmap_t &node_map, const uint8_t &k);

void plus_procedure(nfa_t &nfa, nfa_stack_t &stack, lmap_t &node_map, const uint8_t &k);

void construct_kgraph(const std::string &postfix, nfa_t &nfa, lmap_t &node_map, const uint8_t &k);
