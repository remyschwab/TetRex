#include "otf_collector.h"


void update_kmer(const int &symbol, kmer_t &kmer, IndexStructure &ibf)
{
    if(ibf.molecule_ == "na")
    {
        uint64_t fb = (symbol>>1)&3;
        uint64_t forward;
        uint64_t reverse;
        forward = ((kmer<<2)&ibf.selection_mask_) | fb;
        reverse = revComplement(forward, ibf.k_);
        kmer = forward <= reverse ? forward : reverse;
        return;
    }
}


void update_path(auto &current_state, int &symbol, auto &agent, IndexStructure &ibf, cache_t &cache)
{
    bitvector hits;
    if(current_state.shift_count_ < (ibf.k_-1)) // Corresponds to a new kmer < threshold size (--A)
    {
        update_kmer(symbol, current_state.kmer_, ibf);
        current_state.shift_count_++;
    }
    else if(current_state.shift_count_ == (ibf.k_-1)) // If the kmer just needs to be updated one more time to be a valid kmer (-AC)
    {
        update_kmer(symbol, current_state.kmer_, ibf);
        if(cache.find(current_state.kmer_) == cache.end())
        {
            hits = agent.bulk_contains(current_state.kmer_);
            cache[current_state.kmer_] = hits;
        }
        // hits = cache[current_state.kmer_];
        current_state.path_.raw_data() &= hits.raw_data();
        current_state.shift_count_++;
    }
    else if(current_state.shift_count_ == (ibf.k_)) // (ACG) Next kmer will be valid
    {
        update_kmer(symbol, current_state.kmer_, ibf);
        if(cache.find(current_state.kmer_) == cache.end())
        {
            hits = agent.bulk_contains(current_state.kmer_);
            cache[current_state.kmer_] = hits;
        }
        hits = cache[current_state.kmer_];
        current_state.path_.raw_data() &= hits.raw_data();
    }
}


bool all_bits_zero(bitvector const & bitvector) noexcept
{
    uint64_t const * const ptr = bitvector.raw_data().data();
    size_t const number_of_words{(bitvector.size() + 63u) >> 6};
    bool result{false};

    for (size_t i{}; !result && i < number_of_words; ++i)
        result |= ptr[i];

    return !result;
}

void split_procedure(const amap_t &arc_map, int &id, auto &top, CustomQueue &minheap, nfa_t &NFA)
{
    node_t n1 = arc_map.at(id).first;
    CollectorsItem item1 = {n1, NFA.id(n1), top.shift_count_, top.kmer_, top.path_};
    minheap.push(item1);
    node_t n2 = arc_map.at(id).second;
    CollectorsItem item2 = {n2, NFA.id(n2), top.shift_count_, top.kmer_, top.path_};
    minheap.push(item2);
}


bitvector collect_Top(nfa_t &NFA, IndexStructure &ibf, lmap_t &nfa_map, const std::vector<int> &rank_map, const amap_t &arc_map)
{
    bitvector path_matrix{ibf.getBinCount()};
    CustomQueue minheap(rank_map, NFA, ibf.k_);

    cache_t kmer_cache;
    auto && ibf_ref = ibf.getIBF();
    auto agent = ibf_ref.membership_agent();
    bitvector hit_vector{ibf.getBinCount()};
    std::fill(hit_vector.begin(), hit_vector.end(), true);
    
    int id = 0;
    uint64_t kmer_init = 0;
    node_t graph_head = NFA.nodeFromId(id);
    CollectorsItem item = {graph_head, id, 0, kmer_init, hit_vector};
    minheap.push(item);
    
    node_t next1;

    // size_t loop_count = 0;
    while(!minheap.empty())
    {
        // loop_count++;
        auto top = minheap.top();
        minheap.pop();
        id = top.id_;
        int symbol = nfa_map[top.node];
        // seqan3::debug_stream << symbol << std::endl;
        switch(symbol)
        {
            case Match:
                path_matrix.raw_data() |= top.path_.raw_data();
                break;
            case Ghost:
                next1 = arc_map.at(id).first;
                item = {next1, NFA.id(next1), top.shift_count_, top.kmer_, top.path_};
                minheap.push(item);
                break;
            case Split:
                split_procedure(arc_map, id, top, minheap, NFA);
                break;
            default:
                update_path(top, symbol, agent, ibf, kmer_cache);
                if(all_bits_zero(top.path_)) break; // Immediately get rid of deadend paths
                next1 = arc_map.at(id).first;
                item = {next1, NFA.id(next1), top.shift_count_, top.kmer_, top.path_};
                minheap.push(item);
                break;
        }
    }
    // seqan3::debug_stream << path_matrix << std::endl;
    // seqan3::debug_stream << loop_count << std::endl;
    return path_matrix;
}
