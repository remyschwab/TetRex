#pragma once

#include "construction_tools.h"


bool redundancy_test(buffer_t &buffer);

bool twin_test(Subgraph &node_pair);

void twin_procedure(Subgraph &node_pair, buffer_t &buffer, nfa_t &nfa, nfa_stack_t &stack, lmap_t &node_map);

Subgraph add_node(nfa_t &nfa, lmap_t &node_map, const int &symbol);

void default_procedure(buffer_t &buffer, const int symbol, nfa_stack_t &stack);

void concat_procedure(nfa_t &nfa, lmap_t &node_map, nfa_stack_t &stack, amap_t &arc_map, buffer_t &buffer, catsites_t& cats);

void union_procedure(nfa_t &nfa, nfa_stack_t &stack, lmap_t &node_map, amap_t &arc_map, buffer_t &buffer);

void optional_procedure(nfa_t &nfa, nfa_stack_t &stack, lmap_t &node_map, amap_t &arc_map, buffer_t &buffer);

void kleene_procedure(nfa_t &nfa, nfa_stack_t &stack, lmap_t &node_map, const uint8_t &k, amap_t &arc_map, buffer_t &buffer);

void plus_procedure(nfa_t &nfa, nfa_stack_t &stack, lmap_t &node_map, const uint8_t &k, amap_t &arc_map, buffer_t &buffer);

void quant_procedure(nfa_t &nfa, nfa_stack_t &stack, lmap_t &node_map, const uint8_t &k, amap_t &arc_map, const size_t min, const size_t max, catsites_t& cats, buffer_t& buffer);

catsites_t construct_reduced_kgraph(const std::string &postfix, nfa_t &nfa, lmap_t &node_map, amap_t &arc_map, const uint8_t &k);
