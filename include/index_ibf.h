//
// Created by Remy Schwab on 13.09.22.
//

#pragma once

#include "hibf/interleaved_bloom_filter.hpp"


class IBFIndex
{

private:
    size_t bin_count_{1u};
    float fpr_{.25};
    uint8_t hash_count_{1u};
    std::vector<std::string> tech_bins_{};
    seqan::hibf::interleaved_bloom_filter ibf_{};
    std::vector<std::vector<uint64_t>> user_bin_data_{};
    size_t bin_size_{};


public:
    bitvector hits_{};
    seqan::hibf::interleaved_bloom_filter::containment_agent_type agent_{};

    IBFIndex() = default;
    ~IBFIndex() = default;

    explicit IBFIndex(size_t bc, float fpr, uint8_t hc, std::vector<std::string> tech_bins) :
            bin_count_{bc},
            fpr_{fpr},
            hash_count_{hc},
            tech_bins_{std::move(tech_bins)},
            ibf_(seqan::hibf::bin_count{bin_count_},
                seqan::hibf::bin_size{bin_size_},
                 seqan::hibf::hash_function_count{hash_count_}),
            hits_(bin_count_, true)
    {
        user_bin_data_.resize(bin_count_);
    }

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

    void emplace(uint64_t const val, size_t const idx)
    {
        user_bin_data_[idx].push_back(val);
    }

    float getFPR() const
    {
        return fpr_;
    }

    size_t find_largest_bin() const
    {
        size_t max_count = 0;
        for(auto &&bin: user_bin_data_) max_count = bin.size() > max_count ? bin.size() : max_count;
        return max_count;
    }

    void init_ibf()
    {
        size_t max_bin = find_largest_bin();
        size_t bin_size_ = compute_bitcount(max_bin);
        ibf_ = seqan::hibf::interleaved_bloom_filter{seqan::hibf::bin_count{bin_count_}, seqan::hibf::bin_size{bin_size_}, seqan::hibf::hash_function_count{hash_count_}};
        // ibf_ = std::move(tmp_ibf);
        for(size_t i = 0; i < user_bin_data_.size(); ++i)
        {
            seqan::hibf::bin_index idx{i};
            for(auto val: user_bin_data_[i]) ibf_.emplace(val, idx);
        }
    }

    void populate_index(uint8_t const ksize, auto &decomposer, auto &base_ref, const uint8_t& threadcount = 1)
    {
        gzFile handle{};
        kseq_t *record{};
        int status{};
        size_t seq_count{};
        for(size_t i = 0; i < tech_bins_.size(); ++i) // Iterate over bins
        {
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
                decomposer.decompose_record(record_view, i, base_ref);
            }
            kseq_destroy(record);
            gzclose(handle);
        }
        init_ibf();
        seqan3::debug_stream << "Indexed " << seq_count << " sequences across " << bin_count_ << " bins." << std::endl;
        if(bin_count_ == 1)
        {
            seqan3::debug_stream << "[WARNING] The indexed reference library was not split into bins. The TetRex runtime will be significantly slower." << std::endl;
        }
    }

    size_t compute_bitcount(const size_t bfn) const
    {
        double numerator = -static_cast<double>(bfn) * std::log(fpr_);
        double denominator = std::pow(std::log(2), 2);
        size_t bitcount = static_cast<size_t>(std::ceil(numerator / denominator));
        return bitcount;
    }

    void spawn_agent()
    {
        agent_ = ibf_.containment_agent();
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
