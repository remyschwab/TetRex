#pragma once

#include "index_includes.h"


class HIBFIndex
{

private:
    size_t bin_count_;
    size_t max_bin_size_;
    uint8_t hash_count_;
    std::vector<std::string> user_bins_;
    index_structure::HIBF hibf_;


public:
    index_structure::HIBF_Agent agent_{};

    HIBFIndex() = default;

    explicit HIBFIndex(size_t bs, uint8_t hc, std::vector<std::string> user_bins) :
            bin_count_{user_bins.size()},
            max_bin_size_{bs},
            hash_count_{hc},
            user_bins_{user_bins},
            hibf_
            {
                [&]()
                {
                    auto get_user_bin_data = [&](size_t const user_bin_id, seqan::hibf::insert_iterator it)
                    {
                        for (auto value : user_bins[user_bin_id])
                            it = value;
                    };

                    seqan::hibf::config config{.input_fn = get_user_bin_data, .number_of_user_bins = bin_count_, .maximum_fpr = 0.05};
                    return index_structure::HIBF{config};
                }()
            },
            agent_{hibf_.membership_agent()}
    {
        
    }

    auto getBinCount() const
    {
        return bin_count_;
    }

    size_t getBinSize() const
    {
        assert(hibf_.bin_size() == max_bin_size_);
        return max_bin_size_;
    }

    uint8_t getHashCount() const
    {
        assert(hibf_.hash_function_count() == hash_count_);
        return hash_count_;
    }

    // void spawn_agent()
    // {
    //     agent_ = hibf_.membership_agent();
    // }

    index_structure::HIBF& getIBF()
    {
        return hibf_;
    }

    // void emplace(uint64_t &val, size_t &idx)
    // {
    //     hibf_.emplace(val, seqan::hibf::bin_index{idx});
    // }

    void populate_index(uint8_t &k)
    {

    }

    bitvector &populate_bitvector(auto membership_results)
    {
        bitvector bits{bin_count_};
        for(auto && id: membership_results)
        {
            bits[id] = 0b1;
        }
        return bits;
    }

    const bitvector & query(uint64_t &kmer)
    {
        std::vector<uint64_t> stupid_vector{kmer};
        auto &results = agent_.membership_for(stupid_vector, 1u);
        return populate_bitvector(results);
    }

    template<class Archive>
    void serialize(Archive &archive)
    {
        archive(bin_count_, bin_size_, hash_count_, hibf_);
    }
};
