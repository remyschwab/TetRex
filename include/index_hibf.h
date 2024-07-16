#pragma once

#include <functional>
#include <hibf/config.hpp>
#include "hibf/hierarchical_interleaved_bloom_filter.hpp"


class HIBFIndex
{

private:
    size_t bin_count_{};
    float fpr_{};
    uint8_t hash_count_{};
    std::vector<std::string> user_bins_{};
    std::vector<std::vector<uint64_t>> user_bin_data_{};
    seqan::hibf::hierarchical_interleaved_bloom_filter hibf_{};


public:
    bitvector hits_{};
    std::optional<seqan::hibf::hierarchical_interleaved_bloom_filter::membership_agent_type> agent_{};

    HIBFIndex() = default;
    HIBFIndex& operator=(HIBFIndex&& _rhs)
    {
        bin_count_ = _rhs.bin_count_;
        fpr_ = _rhs.fpr_;
        hash_count_ = _rhs.hash_count_;
        user_bins_ = std::move(_rhs.user_bins_);
        user_bin_data_ = std::move(_rhs.user_bin_data_);
        hibf_ = std::move(_rhs.hibf_);
        hits_ = std::move(_rhs.hits_);
        agent_.emplace(hibf_.membership_agent());
        return *this;
    } 

    explicit HIBFIndex(float fpr, uint8_t hc, std::vector<std::string> user_bins) :
            bin_count_{user_bins.size()},
            fpr_{fpr},
            hash_count_{hc},
            user_bins_{std::move(user_bins)},
            hibf_{},
            agent_{hibf_.membership_agent()}
    {
        user_bin_data_.resize(bin_count_);
    }

    std::pair<size_t, size_t> getShape() const
    {
        std::pair<size_t, size_t> shape = std::make_pair(bin_count_, 1ULL);
        return shape;
    }

    // float getFPR() const
    // {
    //     return hibf_.
    // }

    auto getBinCount() const
    {
        return bin_count_;
    }

    float getFPR() const
    {
        return fpr_;
    }

    uint8_t getHashCount() const
    {
        // assert(hibf_.hash_function_count() == hash_count_);
        return hash_count_;
    }

    auto & getIBF()
    {
        return hibf_;
    }

    void emplace(uint64_t const val, size_t const idx)
    {
        user_bin_data_[idx].push_back(val);
    }

    void populate_index(uint8_t const ksize, auto &decomposer, auto &base_ref)
    {
        gzFile handle{};
        kseq_t *record{};
        int status{};
        size_t seq_count{};
        for(size_t i = 0; i < bin_count_; ++i) // Iterate over bins
        {
            // fprintf(stderr, "\r%llu\n", i+1);
            handle = gzopen(user_bins_[i].c_str(), "r");
            record = kseq_init(handle);
            while ((status = kseq_read(record)) >= 0) // Iterate over bin records
            {
                std::string_view record_view = record->seq.s;
                if(record_view.length() < ksize)
                {
                    seqan3::debug_stream << "RECORD TOO SHORT " << record->comment.s << std::endl;
                    continue;
                }
                seq_count++;
                decomposer.decompose_record(record_view, i, base_ref);
            }
            kseq_destroy(record);
            gzclose(handle);
        }
        auto get_user_bin_data = [&](size_t const user_bin_id, seqan::hibf::insert_iterator it)
            {
                for (auto value : user_bin_data_[user_bin_id])
                    it = value;
            };
        seqan3::debug_stream << "Indexed " << seq_count << " sequences across " << bin_count_ << " bins." << std::endl;
        seqan::hibf::config config{.input_fn = get_user_bin_data,
                                   .number_of_user_bins = bin_count_,
                                   .number_of_hash_functions = hash_count_,
                                   .maximum_fpr = fpr_,};
        hibf_ = seqan::hibf::hierarchical_interleaved_bloom_filter{config};
    }

    bitvector populate_bitvector(auto membership_results)
    {
        bitvector hits(bin_count_);
        for(auto && id: membership_results)
        {
            hits[id] = 0b1;
        }
        return hits;
    }

    bitvector query(uint64_t const kmer)
    {
        std::vector<uint64_t> stupid_vector{kmer};
        auto &results = agent_->membership_for(stupid_vector, 1u);
        return populate_bitvector(results);
    }

    void spawn_agent()
    {
        agent_.emplace(hibf_.membership_agent());
    }

    template<class Archive>
    void serialize(Archive &archive)
    {
        archive(bin_count_, fpr_, hash_count_, user_bins_, hibf_);
    }
};
