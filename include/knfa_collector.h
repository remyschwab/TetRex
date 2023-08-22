#include <iostream>
#include <vector>
#include <stack>
#include "utils.h"
#include "robin_hood.h"
#include "nfa_pointer.h"


using kmer_t = uint64_t;
using path_t = bitvector;
using cache_t = robin_hood::unordered_map<uint32_t, bitvector>;

struct CollectionItem
{
    State *nfa_state_;
    kmer_t kmer_;
    path_t path_;
    uint8_t cycles_;
    uint8_t shift_count_;
};

void update_kmer(const int &symbol, kmer_t &kmer, const uint64_t &selector_mask, const uint64_t &threshold);

void update_path(kmer_t &kmer, uint8_t &shift_count, uint64_t &threshold, int &symbol, auto &agent, bitvector &path, const uint64_t &bitmask);

bitvector collect_kNFA(State *NFA, uint8_t &k, IndexStructure &ibf);
