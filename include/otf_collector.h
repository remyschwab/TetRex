#pragma once

#include <iostream>
#include <vector>
#include "robin_hood.h"
#include "index_base.h"
#include "construct_nfa.h"
#include "lemon/core.h"


namespace CollectionUtils
{
    using kmer_t = uint64_t;
    using path_t = bitvector;
    using cache_t = robin_hood::unordered_map<uint64_t, bitvector>;
    using rank_t = std::vector<int>;
    struct CollectorsItem
    {
        node_t node;
        int id_;
        uint8_t shift_count_;
        kmer_t kmer_;
        path_t path_;
    };
    using comp_table_t = std::vector<robin_hood::unordered_map<uint64_t, CollectorsItem>>;
}

template<index_structure::is_valid flavor, molecules::is_molecule mol_t>
class OTFCollector
{
    private:
        std::unique_ptr<nfa_t> NFA_{};
        std::unique_ptr<lmap_t> nfa_map_{};
        size_t node_count_{};
        TetrexIndex<flavor, mol_t> *ibf_{};
        amap_t arc_map_{};
        CollectionUtils::comp_table_t comp_table_{};
        CollectionUtils::rank_t rank_map_{};
        uint64_t submask_{};
        CollectionUtils::cache_t kmer_cache_{};
    
    public:
        OTFCollector() = default;

        explicit OTFCollector(std::unique_ptr<nfa_t> nfa,
                            std::unique_ptr<lmap_t> nfa_map,
                            TetrexIndex<flavor, mol_t> &ibf,
                            CollectionUtils::rank_t &&rank_map,
                            amap_t const &&arc_map) :
                    NFA_(std::move(nfa)),
                    nfa_map_(std::move(nfa_map)),
                    node_count_{NFA_->nodeNum()},
                    ibf_{&ibf},
                    arc_map_{std::move(arc_map)},
                    comp_table_{node_count_},
                    rank_map_{std::move(rank_map)}
        {
            create_selection_bitmask();
            ibf_->spawn_agent(); // Not done by the IBFIndex constructor during deserialization
        }

    std::string kmer2string(uint64_t kmer, uint8_t ksize)
    {
        robin_hood::unordered_map<uint64_t, char> nucleotide_map{{0x0ULL,'A'}, {0x1ULL,'C'}, {0x3ULL,'G'}, {0x2ULL,'T'}};
        std::string kmer_seq = "";
        for(size_t i = 0; i < ksize; ++i)
        {
            uint64_t base = kmer&3;
            kmer_seq += nucleotide_map[base];
            kmer = kmer>>2;
        }
        std::reverse(kmer_seq.begin(), kmer_seq.end());
        return kmer_seq;
    }

    uint64_t extract_hash(CollectionUtils::CollectorsItem &item)
    {
        uint64_t subhash = (item.kmer_ & submask_);
        return subhash;
    }

    void create_selection_bitmask()
    {
        /*
        Example with k=4 creates a bitmask like 
        0b00000000-00000000-00000000-00000000-00000000-00000000-00000000-00111111
        to detect collisions for absorption
        */
        size_t countdown = ibf_->k_-1;
        while(countdown > 0)
        {
            submask_ <<= ibf_->decomposer_.lshift_;
            submask_ |= ibf_->decomposer_.rmask_;
            --countdown;
        }
        // seqan3::debug_stream << submask_ << std::endl;
    }

    void print_ranks() const
    {
        seqan3::debug_stream << rank_map_ << std::endl;
    }

    void push(CollectionUtils::CollectorsItem &item)
    {
        uint64_t subhash = extract_hash(item);
        int idx = rank_map_[item.id_];
        if(comp_table_[idx].find(subhash) == comp_table_[idx].end())
        {
            comp_table_[idx][subhash] = item;
        }
        else
        {
            absorb(subhash, item, idx);
        }
    }

    void absorb(uint64_t &subhash, CollectionUtils::CollectorsItem &item, int idx)
    {
        // seqan3::debug_stream << "ABSORBING" << std::endl;
        comp_table_[idx][subhash].path_ |= item.path_;
    }

    void scrub(int const topid)
    {
        auto it = comp_table_[topid].begin();
        comp_table_[topid].erase(it);
    }

    void update_path(auto &current_state, int &symbol)
    {
        bitvector hits;
        CollectionUtils::kmer_t canonical_kmer = 0;
        if(current_state.shift_count_ < (ibf_->k_-1)) // Corresponds to a new kmer < threshold size (--A)
        {
            ibf_->update_kmer(symbol, current_state.kmer_);
            current_state.shift_count_++;
        }
        else if(current_state.shift_count_ == (ibf_->k_-1)) // If the kmer just needs to be updated one more time to be a valid kmer (-AC)
        {
            // The canonical kmer is returned but the reference kmer is also updated!!!
            canonical_kmer = ibf_->update_kmer(symbol, current_state.kmer_);
            if(kmer_cache_.find(current_state.kmer_) == kmer_cache_.end())
            {
                kmer_cache_[current_state.kmer_] = ibf_->query(canonical_kmer);
            }
            hits = kmer_cache_[current_state.kmer_];
            current_state.path_ &= hits;
            current_state.shift_count_++;
        }
        else if(current_state.shift_count_ == (ibf_->k_)) // (ACG) Next kmer will be valid
        {
            canonical_kmer = ibf_->update_kmer(symbol, current_state.kmer_);
            if(kmer_cache_.find(current_state.kmer_) == kmer_cache_.end())
            {
                kmer_cache_[current_state.kmer_] = ibf_->query(canonical_kmer);
            }
            hits = kmer_cache_[current_state.kmer_];
            current_state.path_ &= hits;
        }
    }

    void split_procedure(int &id, auto &top)
    {
        node_t n1 = arc_map_.at(id).first;
        CollectionUtils::CollectorsItem item1 = {n1, NFA_->id(n1), top.shift_count_, top.kmer_, top.path_};
        push(item1);
        node_t n2 = arc_map_.at(id).second;
        CollectionUtils::CollectorsItem item2 = {n2, NFA_->id(n2), top.shift_count_, top.kmer_, top.path_};
        push(item2);
    }

    static int sumBitVector(bitvector const &bits)
    {
        int sum = 0;
        for(bool bit : bits)
        {
            sum += bit ? 1 : 0;
        }
        return sum;
    }

    bitvector collect()
    {
        bitvector path_matrix(ibf_->getBinCount());
        bitvector hit_vector(ibf_->getBinCount(), true);

        int id = 0;
        CollectionUtils::kmer_t kmer_init = 0;
        node_t next = NFA_->nodeFromId(id);
        CollectionUtils::CollectorsItem item = {next, id, 0, kmer_init, hit_vector};
        push(item);
        for(size_t i = 0; i < node_count_; ++i)
        {
            while(!comp_table_[i].empty()) // [hash, itm]
            {
                auto top = comp_table_[i].begin()->second;
                scrub(i);
                id = top.id_;
                int symbol = (*nfa_map_)[top.node];
                switch(symbol)
                {
                    case Match:
                        path_matrix |= top.path_;
                        break;
                    case Ghost:
                        next = arc_map_.at(id).first;
                        item = {next, NFA_->id(next), top.shift_count_, top.kmer_, top.path_};
                        push(item);
                        break;
                    case Split:
                        split_procedure(id, top);
                        break;
                    default:
                        update_path(top, symbol);
                        if(top.path_.none()) break; // Immediately get rid of deadend paths
                        next = arc_map_.at(id).first;
                        item = {next, NFA_->id(next), top.shift_count_, top.kmer_, top.path_};
                        push(item);
                        break;
                }
            }
        }
        // seqan3::debug_stream << path_matrix << std::endl;
        return path_matrix;
    }
};
