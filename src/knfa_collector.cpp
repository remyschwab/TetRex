#include "knfa_collector.h"


void update_kmer(const int &symbol, kmer_t &kmer, const uint64_t &selector_mask, const uint64_t &threshold)
{
    uint64_t fb = (symbol>>1)&3;
    uint64_t forward;
    uint64_t reverse;
    forward = ((kmer<<2)&selector_mask) | fb;
    reverse = revComplement(forward, threshold);
    kmer = forward <= reverse ? forward : reverse;
}


void update_path(kmer_t &kmer, uint8_t &shift_count, uint64_t &threshold, int &symbol, auto &agent, bitvector &path, const uint64_t &bitmask)
{
    bitvector hits;
    if(shift_count < (threshold-1)) // Corresponds to a new kmer < threshold size
    {
        update_kmer(symbol, kmer, bitmask, threshold);
        shift_count++;
    }
    else if(shift_count == (threshold-1))
    {
        update_kmer(symbol, kmer, bitmask, threshold);
        hits = agent.bulk_contains(kmer);
        path.raw_data() &= hits.raw_data();
        shift_count++;
    }
    else if(shift_count == (threshold))
    {
        update_kmer(symbol, kmer, bitmask, threshold);
        hits = agent.bulk_contains(kmer);
        path.raw_data() &= hits.raw_data();
    }
}


bitvector collect_kNFA(State *NFA, uint8_t &k, IndexStructure &ibf)
{
    bitvector path_matrix{ibf.getBinCount()};
    std::fill(path_matrix.begin(), path_matrix.end(), false);
    std::stack<CollectionItem> search_stack;

    bitvector hit_vector{ibf.getBinCount()};
    std::fill(hit_vector.begin(), hit_vector.end(), true);
    uint64_t threshold = k;
    // Spawn IBF membership agent in this scope because it is expensive
    auto && ibf_ref = ibf.getIBF();
    auto agent = ibf_ref.membership_agent();

    CollectionItem stack_init(NFA, 0, hit_vector, 0, 0);
    search_stack.push(stack_init);

    State *current_state;
    kmer_t kmer;
    path_t current_path;
    uint8_t split_hits;
    uint8_t shift_counts;

    CollectionItem item1;
    CollectionItem item2;
    while(!search_stack.empty())
    {
        current_state = search_stack.top().nfa_state_;
        kmer = search_stack.top().kmer_;
        current_path = search_stack.top().path_;
        split_hits = search_stack.top().cycles_;
        shift_counts = search_stack.top().shift_count_;
        search_stack.pop();
        
        int symbol = current_state->c_;
        switch(symbol)
        {
            case Match:
                path_matrix.raw_data() |= current_path.raw_data();
                break;
            case SplitU:
                item1 = {current_state->out1_, kmer, current_path, split_hits, shift_counts};
                item2 = {current_state->out2_, kmer, current_path, split_hits, shift_counts};
                if((current_state->out1_->c_ > SplitU)) item1.cycles_++;
                if((current_state->out2_->c_ > SplitU)) item2.cycles_++;
                search_stack.push(item2);
                search_stack.push(item1);
                break;
            case SplitP:
                if(split_hits == k)
                {
                    item1 = {current_state->out1_, kmer, current_path, 0, shift_counts};
                    search_stack.push(item1);
                    break;
                }
                item1 = {current_state->out1_, kmer, current_path, 0, shift_counts};
                item2 = {current_state->out2_, kmer, current_path, split_hits, shift_counts};
                search_stack.push(item2);
                search_stack.push(item1);
                break;
            case SplitK:
                if(split_hits == (k+1))
                {
                    item1 = {current_state->out1_, kmer, current_path, 0, shift_counts};
                    search_stack.push(item1);
                    break;
                }
                item1 = {current_state->out1_, kmer, current_path, 0, shift_counts};
                item2 = {current_state->out2_, kmer, current_path, split_hits, shift_counts};
                search_stack.push(item2);
                search_stack.push(item1);
                break;
            default:
                if(current_state->out1_->c_ > SplitU) split_hits++;
                update_path(kmer, shift_counts, threshold, symbol, agent, current_path, ibf.selection_mask_);
                item1 = {current_state->out1_, kmer, current_path, split_hits, shift_counts};
                search_stack.push(item1);
                break;
        }
    }
    return path_matrix;
}
