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


// bitvector collect_BFS(nfa_t &NFA, IndexStructure &ibf, lmap_t &nfa_map)
// {
//     bitvector path_matrix{ibf.getBinCount()};
//     std::fill(path_matrix.begin(), path_matrix.end(), false);
    
//     std::queue<CollectionItem> queue;
//     cache_t kmer_cache;

//     bitvector hit_vector{ibf.getBinCount()};
//     std::fill(hit_vector.begin(), hit_vector.end(), true);
//     // Spawn IBF membership agent in this scope because it is expensive
//     auto && ibf_ref = ibf.getIBF();
//     auto agent = ibf_ref.membership_agent();

//     CollectionItem queue_init(NFA, 0, hit_vector, 0, 0);
//     queue.emplace(queue_init);

//     State *current_state;
//     kmer_t kmer;
//     path_t current_path;
//     uint8_t shift_counts;

//     CollectionItem item1;
//     CollectionItem item2;

//     while(!queue.empty())
//     {
//         auto front = std::move(queue.front());
//         queue.pop();
//         current_state = front.nfa_state_;
//         kmer = front.kmer_;
//         current_path = front.path_;
//         shift_counts = front.shift_count_;
//         queue.pop();
        
//         int symbol = current_state->c_;
//         switch(symbol)
//         {
//             case Match:
//                 path_matrix.raw_data() |= current_path.raw_data();
//                 break;
//             case SplitU:
//                 item1 = {current_state->out1_, kmer, current_path, split_hits, shift_counts};
//                 item2 = {current_state->out2_, kmer, current_path, split_hits, shift_counts};
//                 if((current_state->out1_->c_ > SplitU)) item1.cycles_++;
//                 if((current_state->out2_->c_ > SplitU)) item2.cycles_++;
//                 queue.push(item1);
//                 queue.push(item2);
//                 break;
//             case SplitP:
//                 if(split_hits == (ibf.k_-1)) // The k-1 removes redundant homopolymer kmers
//                 {
//                     item1 = {current_state->out1_, kmer, current_path, 0, shift_counts};
//                     queue.push(item1);
//                     break;
//                 }
//                 item1 = {current_state->out1_, kmer, current_path, 0, shift_counts};
//                 item2 = {current_state->out2_, kmer, current_path, split_hits, shift_counts};
//                 queue.push(item1);
//                 queue.push(item2);
//                 break;
//             case SplitK:
//                 if(split_hits == ibf.k_)
//                 {
//                     item1 = {current_state->out1_, kmer, current_path, 0, shift_counts};
//                     queue.push(item1);
//                     break;
//                 }
//                 item1 = {current_state->out1_, kmer, current_path, 0, shift_counts};
//                 item2 = {current_state->out2_, kmer, current_path, split_hits, shift_counts};
//                 queue.push(item1);
//                 queue.push(item2);
//                 break;
//             default:
//                 if(current_state->out1_->c_ > SplitU) split_hits++;
//                 update_path(kmer, shift_counts, symbol, agent, current_path, ibf, kmer_cache);
//                 if(all_bits_zero(current_path)) break; // Immediately get rid of deadend paths
//                 item1 = {current_state->out1_, kmer, current_path, split_hits, shift_counts};
//                 queue.push(item1);
//                 break;
//         }
//     }
//     seqan3::debug_stream << path_matrix << std::endl;
//     return path_matrix;
// }


bitvector collect_Top(nfa_t &NFA, IndexStructure &ibf, lmap_t &nfa_map, const std::vector<int> &rank_map, const amap_t &arc_map)
{
    bitvector path_matrix{ibf.getBinCount()};
    CustomCompare ranker(rank_map);
    minheap_t minheap(ranker);

    cache_t kmer_cache;
    auto && ibf_ref = ibf.getIBF();
    auto agent = ibf_ref.membership_agent();
    bitvector hit_vector{ibf.getBinCount()};
    std::fill(hit_vector.begin(), hit_vector.end(), true);
    
    int id = 0;
    uint64_t kmer_init = 0;
    node_t graph_head = NFA.nodeFromId(id);
    CollectorsItem item = {graph_head, id, 0, kmer_init, hit_vector};
    CollectorsItem item2; // For Splits...
    minheap.push(item);
    
    node_t next;
    
    while(!minheap.empty())
    {
        auto top = minheap.top();
        minheap.pop();
        id = top.id_;
        std::cout << id << std::endl;
        int symbol = nfa_map[top.node];
        seqan3::debug_stream << id << " " << symbol << std::endl;
        switch(symbol)
        {
            case Match:
                path_matrix.raw_data() |= top.path_.raw_data();
                break;
            case Ghost:
                next = *(arc_map.at(id).first);
                id = NFA.id(next);
                item = {next, id, top.shift_count_, top.kmer_, top.path_};
                minheap.push(item);
                break;
            case SplitU:
                next = *(arc_map.at(id).first);
                id = NFA.id(next);
                item = {next, id, top.shift_count_, top.kmer_, top.path_};
                minheap.push(item);
                next = *(arc_map.at(id).second);
                id = NFA.id(next);
                item = {next, id, top.shift_count_, top.kmer_, top.path_};
                minheap.push(item);
                break;
            case SplitP:
                next = *(arc_map.at(id).first);
                id = NFA.id(next);
                item = {next, id, top.shift_count_, top.kmer_, top.path_};
                minheap.push(item);
                next = *(arc_map.at(id).second);
                id = NFA.id(next);
                item = {next, id, top.shift_count_, top.kmer_, top.path_};
                minheap.push(item);
                break;
            case SplitK:
                next = *(arc_map.at(id).first);
                id = NFA.id(next);
                item = {next, id, top.shift_count_, top.kmer_, top.path_};
                minheap.push(item);
                next = *(arc_map.at(id).second);
                id = NFA.id(next);
                item = {next, id, top.shift_count_, top.kmer_, top.path_};
                minheap.push(item);
                break;
            default:
                update_path(top, symbol, agent, ibf, kmer_cache);
                if(all_bits_zero(top.path_)) break; // Immediately get rid of deadend paths
                next = *(arc_map.at(id).first);
                id = NFA.id(next);
                item = {next, id, top.shift_count_, top.kmer_, top.path_};
                minheap.push(item);
                break;
        }
    }
    seqan3::debug_stream << path_matrix << std::endl;
    return path_matrix;
}
