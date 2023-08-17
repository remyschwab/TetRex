#include <iostream>
#include <vector>
#include <stack>
#include "robin_hood.h"
#include "nfa_pointer.h"


using kmer_t = std::string;
using path_t = robin_hood::unordered_set<kmer_t>;
using cache_t = robin_hood::unordered_map<uint32_t, robin_hood::unordered_set<kmer_t>>;

struct CollectionItem
{
    State *nfa_state_;
    kmer_t kmer_;
    path_t path_;
    uint8_t cycles_;
};

void update_path(kmer_t &kmer, path_t &path_ref, uint64_t &threshold, int &symbol);

void collect_kNFA(State *NFA, uint8_t &k);
