#pragma once

#include <stack>

#include "lemon/smart_graph.h"
#include "lemon/graph_to_eps.h"


using nfa_t = lemon::SmartDigraph;
using node_t = lemon::SmartDigraph::Node;
using lmap_t = lemon::SmartDigraph::NodeMap<int>;
using node_pair_t = std::pair<node_t, node_t>;
using nfa_stack_t = std::stack<node_pair_t>;

enum
{
    Match = 256,
    Ghost = 257,
    SplitU = 258, // Union Split
    SplitP = 259, // + Split
    SplitK = 260, // Kleene Split
};

void export_nfa_img(nfa_t &nfa);

node_pair_t copy_subgraph(node_pair_t &node_pair);

void default_procedure(nfa_t &nfa, nfa_stack_t &stack, lmap_t &node_map, const int &symbol);

void concat_procedure(nfa_t &nfa, nfa_stack_t &stack);

void union_procedure(nfa_t &nfa, nfa_stack_t &stack, lmap_t &node_map);

void optional_procedure(nfa_t &nfa, nfa_stack_t &stack, lmap_t &node_map);

void kleene_procedure(nfa_t &nfa, nfa_stack_t &stack, lmap_t &node_map, const uint8_t &k);

void plus_procedure(nfa_t &nfa, nfa_stack_t &stack, lmap_t &node_map, const uint8_t &k);

void construct_graph(const std::string &postfix, nfa_t &nfa, lmap_t &node_map, const uint8_t &k);
