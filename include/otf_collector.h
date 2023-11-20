#pragma once

#include <iostream>
#include <vector>
#include "custom_queue.h"
#include "index.h"
#include "construct_nfa.h"


void update_kmer(const int &symbol, kmer_t &kmer, IndexStructure &ibf);

void update_path(auto &current_state, int &symbol, auto &agent, IndexStructure &ibf, cache_t &cache);

bool all_bits_zero(bitvector const & bitvector) noexcept;

void split_procedure(const amap_t &arc_map, int &id, auto &top, CustomQueue &minheap, nfa_t &NFA);

bitvector collect_Top(nfa_t &NFA, IndexStructure &ibf, lmap_t &nfa_map, const std::vector<int> &rank_map, const amap_t &arc_map);
