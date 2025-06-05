#pragma once

#include "construction_tools.h"


void default_procedure(nfa_t &nfa, nfa_stack_t &stack, lmap_t &node_map, const int &symbol);

void concat_procedure(nfa_t &nfa, lmap_t &node_map, nfa_stack_t &stack, amap_t &arc_map, catsites_t& cats);

void union_procedure(nfa_t &nfa, nfa_stack_t &stack, lmap_t &node_map, amap_t &arc_map);

void optional_procedure(nfa_t &nfa, nfa_stack_t &stack, lmap_t &node_map, amap_t &arc_map);

void kleene_procedure(nfa_t &nfa, nfa_stack_t &stack, lmap_t &node_map, const uint8_t &k, amap_t &arc_map);

void plus_procedure(nfa_t &nfa, nfa_stack_t &stack, lmap_t &node_map, const uint8_t &k, amap_t &arc_map);

std::pair<size_t, size_t> parse_quant(const std::string& postfix, size_t quant_start);

void quant_procedure(nfa_t &nfa, nfa_stack_t &stack, lmap_t &node_map, const uint8_t &k, amap_t &arc_map, catsites_t& cats);

catsites_t construct_kgraph(const std::string &postfix, nfa_t &nfa, lmap_t &node_map, amap_t &arc_map, const uint8_t &k);
