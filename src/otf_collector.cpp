// #include "otf_collector.h"


// void update_kmer(const int &symbol, kmer_t &kmer, IndexStructure &ibf)
// {
//     if(ibf.molecule_ == "na")
//     {
//         uint64_t fb = (symbol>>1)&3;
//         uint64_t forward;
//         uint64_t reverse;
//         forward = ((kmer<<2)&ibf.selection_mask_) | fb;
//         reverse = revComplement(forward, ibf.k_);
//         kmer = forward <= reverse ? forward : reverse;
//         return;
//     }
// }


// void update_path(kmer_t &kmer, uint8_t &shift_count, int &symbol, auto &agent, bitvector &path, IndexStructure &ibf, cache_t &cache)
// {
//     bitvector hits;
//     if(shift_count < (ibf.k_-1)) // Corresponds to a new kmer < threshold size
//     {
//         update_kmer(symbol, kmer, ibf);
//         shift_count++;
//     }
//     else if(shift_count == (ibf.k_-1))
//     {
//         update_kmer(symbol, kmer, ibf);
//         if(cache.find(kmer) == cache.end())
//         {
//             hits = agent.bulk_contains(kmer);
//             cache[kmer] = hits;
//         }
//         hits = cache[kmer];
//         path.raw_data() &= hits.raw_data();
//         shift_count++;
//     }
//     else if(shift_count == (ibf.k_))
//     {
//         update_kmer(symbol, kmer, ibf);
//         if(cache.find(kmer) == cache.end())
//         {
//             hits = agent.bulk_contains(kmer);
//             cache[kmer] = hits;
//         }
//         hits = cache[kmer];
//         path.raw_data() &= hits.raw_data();
//     }
// }


// bool all_bits_zero(bitvector const & bitvector) noexcept
// {
//     uint64_t const * const ptr = bitvector.raw_data().data();
//     size_t const number_of_words{(bitvector.size() + 63u) >> 6};
//     bool result{false};

//     for (size_t i{}; !result && i < number_of_words; ++i)
//         result |= ptr[i];

//     return !result;
// }


// bitvector collect_BFS(State *NFA, IndexStructure &ibf)
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
//     uint8_t split_hits;
//     uint8_t shift_counts;

//     CollectionItem item1;
//     CollectionItem item2;

//     while(!queue.empty())
//     {
//         current_state = queue.front().nfa_state_;
//         kmer = queue.front().kmer_;
//         current_path = queue.front().path_;
//         split_hits = queue.front().cycles_;
//         shift_counts = queue.front().shift_count_;
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