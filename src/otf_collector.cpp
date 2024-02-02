#include "otf_collector.h"


void print_in_order(size_t &node_count, const std::vector<int> &ranks)
{
    CustomCompare ranker(ranks);
    minheap_t minheap(ranker);
    for(size_t i = 0; i < node_count; ++i)
    {
        std::cout << i << " ";
        minheap.push(i);   
    }
    std::cout << std::endl;
    while(!minheap.empty())
    {
        std::cout << minheap.top() << " ";
        minheap.pop();
    }
    std::cout << std::endl;
}


// template<index_structure::is_valid flavor, molecules::is_molecule mol_t>
// void update_kmer(const int &symbol, kmer_t &kmer, TetrexIndex<flavor, mol_t> &ibf)
// {
//     if(ibf.molecule_ == "na")
//     {
//         uint64_t fb = (symbol>>1)&3;
//         uint64_t forward = ((kmer<<2)&ibf.selection_mask_) | fb;
//         kmer = forward;
//         return;
//     }
// }

template<index_structure::is_valid flavor, molecules::is_molecule mol_t>
void update_path(auto &current_state, int &symbol, auto &agent, TetrexIndex<flavor, mol_t> &ibf, cache_t &cache)
{
    bitvector hits;
    uint64_t canonical_kmer = 0;
    uint64_t reverse = 0;
    if(current_state.shift_count_ < (ibf.k_-1)) // Corresponds to a new kmer < threshold size (--A)
    {
        ibf.update_kmer(symbol, current_state.kmer_, ibf);
        current_state.shift_count_++;
    }
    else if(current_state.shift_count_ == (ibf.k_-1)) // If the kmer just needs to be updated one more time to be a valid kmer (-AC)
    {
        ibf.update_kmer(symbol, current_state.kmer_, ibf);
        reverse = revComplement(current_state.kmer_, ibf.k_);
        canonical_kmer = current_state.kmer_ <= reverse ? current_state.kmer_ : reverse;
        if(cache.find(current_state.kmer_) == cache.end())
        {
            cache[current_state.kmer_] = agent.bulk_contains(canonical_kmer);
        }
        hits = cache[current_state.kmer_];
        current_state.path_.data() &= hits.data();
        current_state.shift_count_++;
    }
    else if(current_state.shift_count_ == (ibf.k_)) // (ACG) Next kmer will be valid
    {
        ibf.update_kmer(symbol, current_state.kmer_, ibf);
        reverse = revComplement(current_state.kmer_, ibf.k_);
        canonical_kmer = current_state.kmer_ <= reverse ? current_state.kmer_ : reverse;
        if(cache.find(current_state.kmer_) == cache.end())
        {
            cache[current_state.kmer_] = agent.bulk_contains(canonical_kmer);
        }
        hits = cache[current_state.kmer_];
        current_state.path_.data() &= hits.data();
    }
}


// bool all_bits_zero(bitvector const & bitvector) noexcept
// {
//     uint64_t const * const ptr = bitvector.raw_data().data();
//     size_t const number_of_words{(bitvector.size() + 63u) >> 6};
//     bool result{false};

//     for (size_t i{}; !result && i < number_of_words; ++i)
//         result |= ptr[i];

//     return !result;
// }


void split_procedure(const amap_t &arc_map, int &id, auto &top, CustomQueue &minheap, nfa_t &NFA)
{
    node_t n1 = arc_map.at(id).first;
    CollectorsItem item1 = {n1, NFA.id(n1), top.shift_count_, top.kmer_, top.path_};
    minheap.push(item1);
    node_t n2 = arc_map.at(id).second;
    CollectorsItem item2 = {n2, NFA.id(n2), top.shift_count_, top.kmer_, top.path_};
    minheap.push(item2);
}


template<index_structure::is_valid flavor, molecules::is_molecule mol_t>
bitvector collect_Top(nfa_t &NFA, TetrexIndex<flavor, mol_t> &ibf, lmap_t &nfa_map, const std::vector<int> &rank_map, const amap_t &arc_map)
{
    bitvector path_matrix{ibf.getBinCount()};
    CustomQueue minheap(rank_map, NFA, ibf.k_);

    cache_t kmer_cache;
    auto ibf_ref = ibf.getIBF();
    auto agent = ibf_ref.membership_agent();
    bitvector hit_vector(ibf.getBinCount(), true);
    std::fill(hit_vector.begin(), hit_vector.end(), true);
    
    int id = 0;
    uint64_t kmer_init = 0;
    node_t next = NFA.nodeFromId(id);
    CollectorsItem item = {next, id, 0, kmer_init, hit_vector};
    minheap.push(item);

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
                path_matrix.data() |= top.path_.data();
                break;
            case Ghost:
                next = arc_map.at(id).first;
                item = {next, NFA.id(next), top.shift_count_, top.kmer_, top.path_};
                minheap.push(item);
                break;
            case Split:
                split_procedure(arc_map, id, top, minheap, NFA);
                break;
            default:
                update_path(top, symbol, agent, ibf, kmer_cache);
                if(top.path_.none()) break; // Immediately get rid of deadend paths
                next = arc_map.at(id).first;
                item = {next, NFA.id(next), top.shift_count_, top.kmer_, top.path_};
                minheap.push(item);
                break;
        }
    }
    // seqan3::debug_stream << path_matrix << std::endl;
    // seqan3::debug_stream << loop_count << std::endl;
    return path_matrix;
}
