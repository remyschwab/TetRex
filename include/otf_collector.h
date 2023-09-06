#include <iostream>
#include <vector>
#include <stack>
#include <queue>
#include "utils.h"
#include "robin_hood.h"
#include "nfa_pointer.h"


using kmer_t = uint64_t;
using path_t = bitvector;
using cache_t = robin_hood::unordered_map<uint64_t, bitvector>;

struct CollectionItem
{
    State *nfa_state_;
    kmer_t kmer_;
    path_t path_;
    uint8_t cycles_;
    uint8_t shift_count_;
};

void update_kmer(const int &symbol, kmer_t &kmer, IndexStructure &ibf);

void update_path(kmer_t &kmer, uint8_t &shift_count, int &symbol, auto &agent, bitvector &path, IndexStructure &ibf, cache_t &cache);

bool all_bits_zero(bitvector const & bitvector) noexcept;

void condense_queue(std::queue<CollectionItem> &queue);

bitvector collect_DFS(State *NFA, IndexStructure &ibf);

bitvector collect_BFS(State *NFA, IndexStructure &ibf);
