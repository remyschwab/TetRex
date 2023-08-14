#include <iostream>
#include <vector>
#include <stack>
#include "robin_hood.h"
#include "nfa_pointer.h"


using kmer_t = std::string;
using path_t = robin_hood::unordered_set<kmer_t>;

struct CollectionItem
{
    State *nfa_state_;
    kmer_t kmer_;
    path_t path_;
};

kmer_t update_kmer(kmer_t kmer, uint64_t &threshold, int &symbol);

void update_path(kmer_t &kmer, path_t &path_ref, uint64_t &threshold, int &symbol);

void collect_kNFA(State *NFA, uint8_t &k);
