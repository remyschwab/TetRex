#pragma once

#include <iostream>
#include <vector>
#include "custom_queue.h"
#include "index_base.h"
#include "construct_nfa.h"


inline void print_in_order(size_t &node_count, const std::vector<int> &ranks)
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


template<index_structure::is_valid flavor, molecules::is_molecule mol_t>
void update_path(auto &current_state, int &symbol, TetrexIndex<flavor, mol_t> &ibf, cache_t &cache)
{
    bitvector hits;
    uint64_t canonical_kmer = 0;
    if(current_state.shift_count_ < (ibf.k_-1)) // Corresponds to a new kmer < threshold size (--A)
    {
        ibf.update_kmer(symbol, current_state.kmer_);
        current_state.shift_count_++;
    }
    else if(current_state.shift_count_ == (ibf.k_-1)) // If the kmer just needs to be updated one more time to be a valid kmer (-AC)
    {
        // The canonical kmer is returned but the reference kmer is also updated!!!
        canonical_kmer = ibf.update_kmer(symbol, current_state.kmer_);
        if(cache.find(current_state.kmer_) == cache.end())
        {
            cache[current_state.kmer_] = ibf.query(canonical_kmer);
        }
        hits = cache[current_state.kmer_];
        current_state.path_ &= hits;
        current_state.shift_count_++;
    }
    else if(current_state.shift_count_ == (ibf.k_)) // (ACG) Next kmer will be valid
    {
        canonical_kmer = ibf.update_kmer(symbol, current_state.kmer_);
        if(cache.find(current_state.kmer_) == cache.end())
        {
            cache[current_state.kmer_] = ibf.query(canonical_kmer);
        }
        hits = cache[current_state.kmer_];
        current_state.path_ &= hits;
    }
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


template<index_structure::is_valid flavor, molecules::is_molecule mol_t>
bitvector collect_Top(nfa_t &NFA, TetrexIndex<flavor, mol_t> &ibf, lmap_t &nfa_map, const std::vector<int> &rank_map, const amap_t &arc_map)
{
    bitvector path_matrix(ibf.getBinCount());
    CustomQueue minheap(rank_map, NFA, ibf.k_);
    seqan3::debug_stream << path_matrix << std::endl;
    cache_t kmer_cache;
    bitvector hit_vector(ibf.getBinCount(), true);
    // std::fill(hit_vector.begin(), hit_vector.end(), true);
    
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
                path_matrix |= top.path_;
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
                update_path(top, symbol, ibf, kmer_cache);
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
