#pragma once

#include <iostream>
#include <vector>
#include <queue>
#include "index.h"
#include "construct_nfa.h"



// Define custom lambda comparator for the min-heap that uses the topological rankings
const auto customComparator = [](const int& a, const int& b, const std::vector<int> &rankings) {
    return rankings[a] > rankings[b];
};

using kmer_t = uint64_t;
using path_t = bitvector;
using cache_t = robin_hood::unordered_map<uint64_t, bitvector>;
using comp_table_t = std::vector<std::vector<std::pair<kmer_t, bitvector>>>;
using minheap_t = std::priority_queue<int, std::vector<int>, std::function<bool(int, int)>>;

struct CollectionItem
{
    node_t node;
    kmer_t kmer_;
    path_t path_;
    uint8_t shift_count_;
};

void update_kmer(const int &symbol, kmer_t &kmer, IndexStructure &ibf);

void update_path(kmer_t &kmer, uint8_t &shift_count, int &symbol, auto &agent, bitvector &path, IndexStructure &ibf, cache_t &cache);

bool all_bits_zero(bitvector const & bitvector) noexcept;

void condense_queue(std::queue<CollectionItem> &queue);

bitvector collect_BFS(nfa_t &NFA, IndexStructure &ibf, lmap_t &nfa_map);

bitvector collect_TOP(nfa_t &NFA, IndexStructure &ibf, lmap_t &nfa_map, const std::vector<int> &rank_map, const amap_t &arc_map);
