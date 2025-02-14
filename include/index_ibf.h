//
// Created by Remy Schwab on 13.09.22.
//

#pragma once

#include "hibf/interleaved_bloom_filter.hpp"

KSEQ_INIT(gzFile, gzread)


class IBFIndex
{

private:
    size_t bin_count_{};
    size_t bin_size_{};
    uint8_t hash_count_{};
    std::vector<std::string> tech_bins_{};
    seqan::hibf::interleaved_bloom_filter ibf_{};


public:
    bitvector hits_{};
    seqan::hibf::interleaved_bloom_filter::membership_agent_type agent_{};

    IBFIndex() = default;

    explicit IBFIndex(size_t bin_size, uint8_t hc, std::vector<std::string> tech_bins, size_t bc) :
            bin_count_{bc},
            bin_size_{bin_size},
            hash_count_{hc},
            tech_bins_{std::move(tech_bins)},
            ibf_(seqan::hibf::bin_count{bin_count_},
             seqan::hibf::bin_size{bin_size_},
              seqan::hibf::hash_function_count{hash_count_}),
            hits_(bin_count_, true),
            agent_{ibf_.membership_agent()}
    { }

    std::pair<size_t, size_t> getShape() const
    {
        return std::make_pair(getBinCount(), getBinSize());
    }

    size_t getBinCount() const
    {
        assert(ibf_.bin_count() == bin_count_);
        return bin_count_;
    }

    size_t getBinSize() const
    {
        assert(ibf_.bin_size() == bin_size_);
        return bin_size_;
    }

    uint8_t getHashCount() const
    {
        assert(ibf_.hash_function_count() == hash_count_);
        return hash_count_;
    }

    auto & getIBF()
    {
        return ibf_;
    }

    void emplace(uint64_t const val, seqan::hibf::bin_index const idx)
    {
        ibf_.emplace(val, idx);
    }

    float getFPR() const
    {
        return 0.05; // Lol gotta fix this I guess...
    }

    void populate_index(uint8_t const ksize, auto &decomposer, auto &base_ref)
    {
        gzFile handle{};
        kseq_t *record{};
        int status{};
        size_t seq_count{};
        for(size_t i = 0; i < tech_bins_.size(); ++i) // Iterate over bins
        {
            seqan::hibf::bin_index idx{i};
            handle = gzopen(tech_bins_[i].c_str(), "r");
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
                decomposer.decompose_record(record_view, idx, base_ref);
            }
            kseq_destroy(record);
            gzclose(handle);
        }
        seqan3::debug_stream << "Indexed " << seq_count << " sequences across " << bin_count_ << " bins." << std::endl;
        if(bin_count_ == 1)
        {
            seqan3::debug_stream << "[WARNING] The indexed reference library was not split into bins. The TetRex runtime will be significantly slower." << std::endl;
        }
    }

    void spawn_agent()
    {
        agent_ = ibf_.membership_agent();
    }

    bitvector & query(uint64_t const kmer)
    {
        hits_ = agent_.bulk_contains(kmer);
        return hits_;
    }

    template<class Archive>
    void serialize(Archive &archive)
    {
        archive(bin_count_, bin_size_, hash_count_, tech_bins_, ibf_);
    }
};
