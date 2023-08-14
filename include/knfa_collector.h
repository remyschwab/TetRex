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

void collect_kNFA(State *NFA, uint8_t &k);
