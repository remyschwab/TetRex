#pragma once

#include <iostream>
#include <vector>
#include "robin_hood.h"
#include "index_base.h"
#include "construct_nfa.h"
#include "lemon/core.h"
#include "lemon/list_graph.h"
#include "dGramIndex.h"


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
        bool gapped_{false};
    };
    using comp_table_t = std::vector<robin_hood::unordered_map<uint64_t, CollectorsItem>>;
    using cmplx_t = std::vector<std::vector<int>>;
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
        std::unique_ptr<gmap_t> gap_map_{};
        bool has_dibf_{};
        DGramIndex *dibf_{};
        CollectionUtils::cmplx_t cmplx_mtrx_{};
        std::vector<size_t> rank_to_id_{};
    
    public:
        OTFCollector() = default;

        explicit OTFCollector(std::unique_ptr<nfa_t> nfa,
                            std::unique_ptr<lmap_t> nfa_map,
                            TetrexIndex<flavor, mol_t> &ibf,
                            amap_t const &&arc_map, std::unique_ptr<gmap_t>gap_map,
                            const bool& has_dibf,
                            DGramIndex& dibf) :
                    NFA_(std::move(nfa)),
                    nfa_map_(std::move(nfa_map)),
                    node_count_{NFA_->nodeNum()},
                    ibf_{&ibf},
                    arc_map_{std::move(arc_map)},
                    comp_table_{},
                    rank_map_{},
                    submask_{},
                    kmer_cache_{},
                    gap_map_{std::move(gap_map)},
                    has_dibf_{has_dibf},
                    dibf_{&dibf},
                    cmplx_mtrx_{}
        {
            create_selection_bitmask();
            comp_table_.resize(node_count_);
            determine_top_sort();
            ibf_->spawn_agent(); // Not done by the IBFIndex constructor during deserialization
            if(has_dibf_) dibf_->spawn_agent();
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

    void update_downstream_counts(const uint8_t k, std::vector<int> &local_arr, std::vector<int> &down_arr)
    {
        for(size_t i = 0; i < k; ++i) down_arr[i] += local_arr[i];
    }

    void update_local_counts(std::vector<int> &arr, const uint8_t k, size_t &ttl_cmplx)
    {
        for(size_t i = k-1; i > 0; --i) arr[i] = arr[i-1];
        arr[0] = 1;
        ttl_cmplx += arr[k-1];
    }

    void build_rank_to_id_map()
    {
        for(auto i = 0; i < rank_map_.size(); ++i) rank_to_id_[rank_map_[i]] = i; // I don't need hashing here
    }

    CollectionUtils::cmplx_t compute_complexity(const uint8_t ksize)
    {
        std::vector<std::vector<int>> counts_matrix(node_count_, std::vector<int>(ksize, 0));
        // robin_hood::unordered_map<int, int> rank_to_id;
        rank_to_id_.resize(rank_map_.size());
        build_rank_to_id_map();
        size_t total_complexity = 0;
        int down_idx;
        for(size_t i = 0; i < node_count_; ++i)
        {
            int symbol = (*nfa_map_)[NFA_->nodeFromId(rank_to_id_[i])];
            switch(symbol)
            {
                case Match:
                    break;
                case Ghost:
                    down_idx = rank_map_[NFA_->id(arc_map_.at(rank_to_id_[i]).first)];
                    update_downstream_counts(ksize, counts_matrix[i], counts_matrix[down_idx]);
                    // seqan3::debug_stream << "•:" << counts_matrix[i] << std::endl;
                    break;
                case Split:
                    down_idx = rank_map_[NFA_->id(arc_map_.at(rank_to_id_[i]).first)];
                    update_downstream_counts(ksize, counts_matrix[i], counts_matrix[down_idx]);
                    down_idx = rank_map_[NFA_->id(arc_map_.at(rank_to_id_[i]).second)];
                    update_downstream_counts(ksize, counts_matrix[i], counts_matrix[down_idx]);
                    // seqan3::debug_stream << "Ø:" << counts_matrix[i] << std::endl;
                    break;
                default:
                    update_local_counts(counts_matrix[i], ksize, total_complexity);
                    down_idx = rank_map_[NFA_->id(arc_map_.at(rank_to_id_[i]).first)];
                    update_downstream_counts(ksize, counts_matrix[i], counts_matrix[down_idx]);
                    // seqan3::debug_stream << static_cast<char>(symbol) << ":" << counts_matrix[i] << std::endl;
                    break;
            }
        }
        // DBG(total_complexity);
        // DBG(counts_matrix);
        return counts_matrix;
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

    void update_path(auto &current_state, const int &symbol)
    {
        bitvector hits;
        CollectionUtils::kmer_t canonical_kmer = 0;
        if(current_state.gapped_)
        {
            if(has_dibf_)
            {
                uint64_t dgram = current_state.kmer_;
                dgram += static_cast<uint64_t>(DGramTools::aa_to_num(symbol));
                hits = dibf_->query(dgram);
                current_state.path_ &= hits;
            }
            current_state.kmer_ = 0u;
            current_state.gapped_ = false;
        }
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

    void gap_procedure(const int& id, const CollectionUtils::CollectorsItem& top)
    {
        uint8_t code_a = top.kmer_ & 0b11111;
        size_t gap = (*gap_map_)[NFA_->nodeFromId(id)];
        uint64_t dgram = 0;
        dgram = static_cast<uint64_t>(gap) * 400ULL + static_cast<uint64_t>(code_a) * 20ULL;
        ibf_->set_stores(0u, 0u); // Not sure I need to do this
        node_t next = arc_map_.at(id).first;
        CollectionUtils::CollectorsItem item = {next, NFA_->id(next), 0, dgram, top.path_, true};
        push(item);
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

    size_t sumComplexity(const CollectionUtils::cmplx_t& mtrx)
    {
        size_t complexity = 0;
        for(const auto& window: mtrx) complexity += window.back();
        return complexity;
    }

    void determine_top_sort()
    {
        wmap_t list(*NFA_);
        lemon::topologicalSort(*NFA_, list);
        rank_map_.resize(NFA_->nodeNum());
        size_t rank = 1; // Start node is always at 0
        for(auto &&it: list.order)
        {
            rank_map_[NFA_->id(it)] = rank;
            rank++;
        }
    }

    uint64_t collect_dgram()
    {
        return 0ULL;
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
                    case 36:
                        next = arc_map_.at(id).first;
                        item = {next, NFA_->id(next), top.shift_count_, top.kmer_, top.path_};
                        push(item);
                        break;
                    case Ghost:
                        next = arc_map_.at(id).first;
                        item = {next, NFA_->id(next), top.shift_count_, top.kmer_, top.path_};
                        push(item);
                        break;
                    case Split:
                        split_procedure(id, top);
                        break;
                    case Gap:
                        gap_procedure(id, top);
                        break;
                    default:
                        update_path(top, symbol);
                        if(top.path_.none()) break; // Immediately get rid of deadend paths
                        next = arc_map_.at(id).first;
                        item = {next, NFA_->id(next), top.shift_count_, top.kmer_, top.path_};
                        push(item);
                        // DBG(top.kmer_);
                        break;
                }
            }
        }
        return path_matrix;
    }

    void add_gap(const Catsite& cat, const size_t gapsize)
    {
        node_t gap_node = NFA_->addNode();
        ++node_count_; // Can't forget this!
        (*nfa_map_)[gap_node] = Gap;
        (*gap_map_)[gap_node] = gapsize;
        update_arc_map(*NFA_, *nfa_map_, arc_map_, cat.cleavage_site_, gap_node);
        update_arc_map(*NFA_, *nfa_map_, arc_map_, gap_node, cat.downstream_);
    }

    std::pair<node_t, node_t> add_guard(const Catsite& cat)
    {
        node_t split_node = NFA_->addNode();
        ++node_count_; // Can't forget this!
        node_t ghost_node = NFA_->addNode();
        ++node_count_; // Can't forget this!
        (*nfa_map_)[split_node] = Split;
        (*nfa_map_)[ghost_node] = Ghost;
        update_arc_map(*NFA_, *nfa_map_, arc_map_, cat.cleavage_site_, split_node);
        update_arc_map(*NFA_, *nfa_map_, arc_map_, ghost_node, cat.downstream_);
        return {split_node, ghost_node};
    }

    void analyze_complexity()
    {
        comp_table_.resize(node_count_);
        determine_top_sort();
        cmplx_mtrx_ = compute_complexity(ibf_->k_);
        return;
    }

    robin_hood::unordered_set<size_t> sumGaps(const Catsite& c1, const Catsite& c2)
    {
        robin_hood::unordered_set<size_t> gaps;
        for(auto &g: c1.gaps_)
        {
            for(auto &gg: c2.gaps_)
            {
                gaps.insert(g+gg);
            }
        }
        return gaps;
    }

    void merge_catsites(catsites_t& cats)
    {
        std::sort(cats.begin(), cats.end(), [&](const auto& a, const auto& b) {
            return rank_map_[a.cleavage_start_id_] < rank_map_[b.cleavage_start_id_];
        });

        
        bool done = false;
        catsites_t merged;

        for (int i = 0; i < cats.size(); i++)
        {
            size_t currentStart = rank_map_[cats[i].cleavage_start_id_];
            size_t currentEnd = rank_map_[cats[i].cleavage_end_id_];
            if (merged.empty() || (currentStart-1 != rank_map_[merged.back().cleavage_end_id_])) // Weird off by one thing here
            {
                merged.push_back(cats[i]);
                continue;
            }
            merged.back().cleavage_end_ = cats[i].cleavage_end_;
            merged.back().cleavage_end_id_ = cats[i].cleavage_end_id_;
            merged.back().gaps_ = sumGaps(merged.back(), cats[i]);
            done = true;
        }
        if(done) cats = merged;
    }

    void augment(catsites_t& catsites)
    {
        // This updates only the external arc map that is used for collection
        // It doesn't deal with the internal, lemon arc map
        // This is acceptable for now because we only replace gaps between concatentation operators
        // Lemon algos like DFS are only being called within subgraphs (I think)
        // for(auto &cat: catsites) cat.dumpInfo(rank_map_);
        merge_catsites(catsites);
        // DBG("\n");
        // seqan3::debug_stream << "[AUGMENTING " << catsites.size() << " EDGE(S)]" << std::endl;
        for(auto &&cat: catsites)
        {
            robin_hood::unordered_set<size_t> gaps = cat.gaps_;
            // DBG(gaps);
            cat.complete(*NFA_, arc_map_);
            // cat.dumpInfo(rank_map_);
            if(gaps.size() == 1) add_gap(cat, *gaps.begin());
            else
            {
                std::pair<node_t, node_t> guard = add_guard(cat);
                cat.cleavage_site_ = guard.first; // NOTE: This doesn't update any of the IDs
                cat.downstream_ = guard.second;   // Should be ok...
                for(const auto& gap: gaps) add_gap(cat, gap);
            }
        }
        comp_table_.resize(node_count_);
        determine_top_sort();
    }

    void draw_graph(const catsites_t& cats, const bool& augment)
    {
        print_graph(*NFA_, *nfa_map_, cats, augment);
    }

    size_t get_node_count()
    {
        return node_count_;
    }
};

// "LMAEGLYNHSVRVRSDIEEDEED"
// "LMAEGLYN.......DIEEDEED"
// "LMAEGLYN.{7}DIEEDEED"

