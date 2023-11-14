#pragma once

#include <iostream>
#include <vector>
#include <queue>
#include "index.h"
#include "construct_nfa.h"


using kmer_t = uint64_t;
using path_t = bitvector;
using cache_t = robin_hood::unordered_map<uint64_t, bitvector>;

struct CollectorsItem
{
    node_t node;
    int id_;
    uint8_t shift_count_;
    kmer_t kmer_;
    path_t path_;
};

class CustomCompare
{
    private:
        std::vector<int> ranks_;

    public:
        bool operator() (CollectorsItem &item1, CollectorsItem &item2)
        {
            return ranks_[item1.id_] > ranks_[item2.id_];
        }

    CustomCompare() = default;

    explicit CustomCompare(std::vector<int> ranks) : ranks_{ranks} {}
};

using minheap_t = std::priority_queue<CollectorsItem, std::vector<CollectorsItem>, CustomCompare>;

void update_kmer(const int &symbol, kmer_t &kmer, IndexStructure &ibf);

void update_path(auto &current_state, int &symbol, auto &agent, IndexStructure &ibf, cache_t &cache);

bool all_bits_zero(bitvector const & bitvector) noexcept;

void split_procedure(const amap_t &arc_map, int &id, auto &top, minheap_t &minheap, CollectorsItem &item, nfa_t &NFA);

void condense_queue(std::queue<CollectorsItem> &queue);

bitvector collect_BFS(nfa_t &NFA, IndexStructure &ibf, lmap_t &nfa_map);

bitvector collect_Top(nfa_t &NFA, IndexStructure &ibf, lmap_t &nfa_map, const std::vector<int> &rank_map, const amap_t &arc_map);
