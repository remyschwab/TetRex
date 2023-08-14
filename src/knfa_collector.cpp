#include "knfa_collector.h"



kmer_t update_kmer(kmer_t kmer, uint64_t &threshold, int &symbol)
{
    kmer_t new_kmer;
    if(kmer.length() < (threshold-1))
    {
        new_kmer += (char)symbol;
    }
    else if(kmer.length() == (threshold-1))
    {
        new_kmer += (char)symbol;
    }
    else if(kmer.length() == threshold)
    {
        new_kmer = kmer.substr(1) + (char)symbol;
    }
    return new_kmer;
}


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
    CollectionItem stack_init(NFA, "", {});
    search_stack.push(stack_init);

    uint64_t threshold = k;

    bool kmer_repeat;
    bool split_ptr;

    CollectionItem item1;
    CollectionItem item2;
    while(!search_stack.empty())
    {
        State *current_state = search_stack.top().nfa_state_;
        kmer_t kmer = search_stack.top().kmer_;
        path_t current_path = search_stack.top().path_;
        search_stack.pop();
        
        int symbol = current_state->c_;
        switch(symbol)
        {
            case Match:
                path_matrix.push_back(current_path);
                break;
            case Split:
                item1.nfa_state_ = current_state->out1_;
                item1.kmer_ = kmer;
                item1.path_ = current_path;
                search_stack.push(item1);
                item2.nfa_state_ = current_state->out2_;
                item2.kmer_ = kmer;
                item2.path_ = current_path;
                search_stack.push(item2);
                break;
            default:
                kmer_repeat = current_path.find(update_kmer(kmer, threshold, symbol)) != current_path.end();
                split_ptr = current_state->out1_->c_ == Split;
                if(kmer_repeat && split_ptr)
                {
                    break;
                }
                update_path(kmer, current_path, threshold, symbol);
                item1 = {current_state->out1_, kmer, current_path};
                search_stack.push(item1);
                break;
        }
    }
    seqan3::debug_stream << path_matrix << std::endl;
}
