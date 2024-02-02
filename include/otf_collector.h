#pragma once

#include <iostream>
#include <vector>
#include "custom_queue.h"
#include "index_base.h"
#include "construct_nfa.h"


void print_in_order(size_t &node_count, const std::vector<int> &ranks);

template<index_structure::is_valid flavor, molecules::is_molecule mol_t>
void update_kmer(const int &symbol, kmer_t &kmer, TetrexIndex<flavor, mol_t> &ibf);

template<index_structure::is_valid flavor, molecules::is_molecule mol_t>
void update_path(auto &current_state, int &symbol, auto &agent, TetrexIndex<flavor, mol_t> &ibf, cache_t &cache);

bool all_bits_zero(bitvector const & bitvector) noexcept;

void split_procedure(const amap_t &arc_map, int &id, auto &top, CustomQueue &minheap, nfa_t &NFA);

template<index_structure::is_valid flavor, molecules::is_molecule mol_t>
bitvector collect_Top(nfa_t &NFA, TetrexIndex<flavor, mol_t> &ibf, lmap_t &nfa_map, const std::vector<int> &rank_map, const amap_t &arc_map);
