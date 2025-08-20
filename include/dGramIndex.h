#pragma once
#include <unordered_map>
#include <stdint.h>
#include <stdexcept>

#include "utils.h"
#include "arg_parse.h"
#include "hibf/interleaved_bloom_filter.hpp"
#include <cereal/types/vector.hpp>
#include "cereal/types/string.hpp"


namespace DGramTools {
    struct Dgram
    {
        unsigned char a, b;
        size_t gap;
    };
}; // End DGramTools

class DGramIndex 
{
    private:
        size_t min_gap_{};
        size_t max_gap_{};
        size_t pad_{};
        uint8_t hc_{};
        float fpr_{};
        std::vector<std::string> bins_{};
        seqan::hibf::interleaved_bloom_filter dibf_{};
        size_t bc_;
        std::unordered_map<char, uint8_t> alpha_map_{};
        // std::vector<std::pair<std::string, uint64_t>> dgram_buffer_;
        std::vector<std::vector<uint64_t>> dgram_buffer_;
        // seqan::hibf::interleaved_bloom_filter::membership_agent_type agent_{};

    public:

        // Rule of 5
        DGramIndex() = default;
        DGramIndex(const DGramIndex &) = default;
        DGramIndex & operator=(const DGramIndex &) = default;
        DGramIndex(DGramIndex&&) noexcept = default;
        DGramIndex & operator=(DGramIndex &&) noexcept = default;
        ~DGramIndex() = default;

        explicit DGramIndex(size_t min_gap, size_t max_gap, size_t pad, size_t hc, float fpr, std::vector<std::string> bins)
            : min_gap_(min_gap), max_gap_(max_gap), pad_(pad), hc_(hc), fpr_(fpr), bins_(bins)
        {
            init_alphabet();
            bc_ = bins_.size();
        }

        void populate()
        {
            dgram_buffer_.resize(bc_);
            FLOOP(bc_)
            {
                add_fasta(i);
            }
            init_ibf();
        }

        // Optional accessor for testing
        const std::vector<std::vector<uint64_t>>& dgrams() const
        {
            return dgram_buffer_;
        }

        void init_alphabet()
        {
            std::string alphabet = "ACDEFGHIKLMNPQRSTVWY";
            for (size_t i = 0; i < alphabet.size(); ++i)
            {
                alpha_map_[alphabet[i]] = static_cast<uint8_t>(i);
            }
        }

        size_t find_largest_bin() const
        {
            size_t max_count = 0;
            for(auto &&bin: dgram_buffer_) max_count = bin.size() > max_count ? bin.size() : max_count;
            return max_count;
        }

        size_t compute_bitcount(const size_t bfn) const
        {
            size_t bitcount;
            double numerator = -static_cast<double>(bfn) * std::log(fpr_);
            double denominator = std::pow(std::log(2), 2);
            bitcount = static_cast<size_t>(std::ceil(numerator / denominator));
            return bitcount;
        }

        void init_ibf()
        {
            size_t max_bin = find_largest_bin();
            size_t bin_size_ = compute_bitcount(max_bin);
            dibf_ = seqan::hibf::interleaved_bloom_filter{seqan::hibf::bin_count{bc_}, seqan::hibf::bin_size{bin_size_}, seqan::hibf::hash_function_count{hc_}};
            for(size_t i = 0; i < dgram_buffer_.size(); ++i)
            {
                seqan::hibf::bin_index idx{i};
                for(auto val: dgram_buffer_[i]) dibf_.emplace(val, idx);
            }
        }

        void process_sequence(const std::string& seq, const size_t bin_id)
        {
            for (size_t i = 0; i + min_gap_ + 1 < seq.size(); ++i)
            {
                char a = seq[i];
                if (alpha_map_.find(a) == alpha_map_.end()) continue;

                for (size_t gap = min_gap_; gap <= max_gap_; ++gap)
                {
                    size_t j = i + gap + 1;
                    if (j >= seq.size()) break;

                    char b = seq[j];
                    if (alpha_map_.find(b) == alpha_map_.end()) continue;

                    uint8_t code_a = alpha_map_[a];
                    uint8_t code_b = alpha_map_[b];
                    uint64_t code = static_cast<uint64_t>(gap) * 400ULL
                                + static_cast<uint64_t>(code_a) * 20ULL
                                + static_cast<uint64_t>(code_b);

                    dgram_buffer_[bin_id].emplace_back(code);
                }
            }
        }

        void add_fasta(size_t bin_id)
        {
            std::filesystem::path filename = bins_[bin_id];
            gzFile fp = gzopen(filename.c_str(), "r");
            if (!fp)
            {
                throw std::runtime_error("Failed to open file: ");
            }

            kseq_t* seq = kseq_init(fp);
            while (kseq_read(seq) >= 0)
            {
                process_sequence(seq->seq.s, bin_id);
            }
            kseq_destroy(seq);
            gzclose(fp);
        }

        template<class Archive>
        void serialize(Archive &archive)
        {
            archive(min_gap_, max_gap_, pad_, hc_, fpr_, bins_, dibf_, bc_);
        }
}; // DGramIndex

template <class DGramIndex>
void save_dindex(const DGramIndex& dindex, const std::filesystem::path& opath)
{
    std::ofstream os{opath, std::ios::binary};
    cereal::BinaryOutputArchive oarchive{os};
    oarchive(dindex);
}

template <class DGramIndex>
void load_dindex(DGramIndex& dindex, const std::filesystem::path& ipath)
{
    std::ifstream is{ipath, std::ios::binary};
    cereal::BinaryInputArchive iarchive{is};
    iarchive(dindex);
}

void drive_dindex(const dindex_arguments &cmd_args);
