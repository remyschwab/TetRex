#include <queue>
#include "robin_hood.h"
#include "construct_nfa.h"


using kmer_t = uint64_t;
using path_t = bitvector;
using cache_t = robin_hood::unordered_map<uint64_t, bitvector>;


struct CollectorsItem
{
    node_t node;
    int id_;
    uint8_t shift_count_;
    kmer_t kmer_;
    path_t path_;
};


using comp_table_t = std::vector<robin_hood::unordered_map<uint64_t, CollectorsItem>>;


class CustomCompare
{
    private:
        std::vector<int> ranks_;

    public:
        bool operator() (int &item1, int &item2) const
        {
            return ranks_[item1] > ranks_[item2];
        }

    CustomCompare() = default;

    explicit CustomCompare(std::vector<int> ranks) : ranks_{ranks} {}

    void print_ranks()
    {
        for(size_t i = 0; i < ranks_.size(); ++i)
        {
            std::cout << i << " " << ranks_[i] << std::endl;
        }
    }
};

using minheap_t = std::priority_queue<int, std::vector<int>, CustomCompare>;


class CustomQueue
{
    private:
        CustomCompare ranker_;
        minheap_t minheap_;
        uint64_t submask_{};
        uint8_t k_;

    public:
        comp_table_t comp_table_;

    CustomQueue() = default;
    ~CustomQueue() = default;

    explicit CustomQueue(std::vector<int> ranks, nfa_t &NFA, uint8_t &k)
    {
        ranker_ = CustomCompare(ranks);
        minheap_ = minheap_t(ranker_);
        comp_table_.resize(NFA.nodeNum());
        k_ = k;
        create_selection_bitmask(k);
    }

    void push(CollectorsItem &item)
    {
        uint64_t subhash = extract_hash(item);
        // seqan3::debug_stream << "Subhash: " << subhash << std::endl;
        if(comp_table_[item.id_].find(subhash) == comp_table_[item.id_].end())
        {
            comp_table_[item.id_][subhash] = item;
            minheap_.push(item.id_);
            seqan3::debug_stream << "+SLOT " << item.id_  << "\t" << generate_kmer_seq(item.kmer_, k_) << "/" << subhash << std::endl;
        }
        else
        {
            seqan3::debug_stream << "ABSORBING " << generate_kmer_seq(item.kmer_, k_) << " AT SLOT " << item.id_ << " WITH SUBHASH " << subhash << "\t";
            absorb(subhash, item);
        }
        // seqan3::debug_stream << minheap_.size() << std::endl;
    }

    void absorb(uint64_t &subhash, CollectorsItem &item)
    {
        seqan3::debug_stream << "[-->" << item.id_ << "<--]" << std::endl;
        comp_table_[item.id_][subhash].path_.data() |= item.path_.data();
    }

    void pop()
    {
        int topid = minheap_.top();
        auto it = comp_table_[topid].begin(); // This creates some randomization
        seqan3::debug_stream << "\t-SLOT " << topid  << "\t" << it->second.kmer_ << "/" << it->first << std::endl;
        // seqan3::debug_stream << std::endl;
        comp_table_[topid].erase(it);
        minheap_.pop();
    }

    CollectorsItem top()
    {
        int heap_top = minheap_.top();
        return comp_table_[heap_top].begin()->second;
    }

    bool empty()
    {
        return minheap_.empty();
    }

    uint64_t extract_hash(CollectorsItem &item)
    {
        uint64_t subhash = (item.kmer_ & submask_);
        return subhash;
    }

    void create_selection_bitmask(uint8_t &k)
    {
        /*
        Example with k=4 creates a bitmask like 
        0b00000000-00000000-00000000-00000000-00000000-00000000-00000000-00111111
        to detect collisions for absorption
        */
        size_t countdown = (k-1);
        while(countdown > 0)
        {
            submask_ = (submask_<<2) | 0b11;
            countdown--;
        }
    }

    void print_ranks()
    {
        ranker_.print_ranks();
    }
};
