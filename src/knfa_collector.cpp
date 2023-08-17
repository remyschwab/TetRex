#include "knfa_collector.h"


void update_path(kmer_t &kmer, path_t &path_ref, uint64_t &threshold, int &symbol)
{
    if(kmer.length() < (threshold-1))
    {
        kmer += (char)symbol;
    }
    else if(kmer.length() == (threshold-1))
    {
        kmer += (char)symbol;
        path_ref.emplace(kmer);
    }
    else if(kmer.length() == threshold)
    {
        kmer = kmer.substr(1) + (char)symbol;
        path_ref.emplace(kmer);
    }
}


void collect_kNFA(State *NFA, uint8_t &k)
{
    std::vector<path_t> path_matrix;
    std::stack<CollectionItem> search_stack;
    CollectionItem stack_init(NFA, "", {}, 0);
    search_stack.push(stack_init);

    uint64_t threshold = k;

    // uint64_t stack_max = 0;

    State *current_state;
    kmer_t kmer;
    path_t current_path;
    uint8_t split_hits;

    CollectionItem item1;
    CollectionItem item2;
    while(!search_stack.empty())
    {
        current_state = search_stack.top().nfa_state_;
        kmer = search_stack.top().kmer_;
        current_path = search_stack.top().path_;
        split_hits = search_stack.top().cycles_;
        search_stack.pop();
        
        int symbol = current_state->c_;
        switch(symbol)
        {
            case Match:
                path_matrix.push_back(current_path);
                break;
            case SplitU:
                item1.nfa_state_ = current_state->out1_;
                item1.kmer_ = kmer;
                item1.path_ = current_path;
                item1.cycles_ = split_hits;
                item2.nfa_state_ = current_state->out2_;
                item2.kmer_ = kmer;
                item2.path_ = current_path;
                item2.cycles_ = split_hits;
                search_stack.push(item2);
                search_stack.push(item1);
                break;
            case SplitP:
                item1.nfa_state_ = current_state->out1_;
                item1.kmer_ = kmer;
                item1.path_ = current_path;
                if(split_hits == k)
                {
                    item1.cycles_ = 0;
                    search_stack.push(item1);
                    break;
                }
                item1.cycles_ = split_hits;
                item2.nfa_state_ = current_state->out2_;
                item2.kmer_ = kmer;
                item2.path_ = current_path;
                item2.cycles_ = split_hits;
                search_stack.push(item2);
                search_stack.push(item1);
                break;
            case SplitK:
                item1.nfa_state_ = current_state->out1_;
                item1.kmer_ = kmer;
                item1.path_ = current_path;
                if(split_hits == (k+1))
                {
                    item1.cycles_ = 0;
                    search_stack.push(item1);
                    break;
                }
                item1.cycles_ = split_hits;
                item2.nfa_state_ = current_state->out2_;
                item2.kmer_ = kmer;
                item2.path_ = current_path;
                item2.cycles_ = split_hits;
                search_stack.push(item2);
                search_stack.push(item1);
                break;
            default:
                if((current_state->out1_->c_ == SplitP) || (current_state->out1_->c_ == SplitK)) split_hits++;
                update_path(kmer, current_path, threshold, symbol);
                item1 = {current_state->out1_, kmer, current_path, split_hits};
                search_stack.push(item1);
                break;
        }
        // stack_max = (stack_max < search_stack.size()) ? search_stack.size() : stack_max;
    }
    seqan3::debug_stream << path_matrix << std::endl;
}
