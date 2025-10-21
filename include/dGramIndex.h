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

    std::vector<std::string> read_input_file_list(const std::filesystem::path &input_file);

    constexpr std::array<unsigned char, UCHAR_MAX+1> make_amino_acid_map()
    {
        std::array<unsigned char, UCHAR_MAX+1> aamap_{};
        aamap_['A'] = 0u;
        aamap_['C'] = 1u;
        aamap_['D'] = 2u;
        aamap_['E'] = 3u;
        aamap_['F'] = 4u;
        aamap_['G'] = 5u;
        aamap_['H'] = 6u;
        aamap_['I'] = 7u;
        aamap_['K'] = 8u;
        aamap_['L'] = 9u;
        aamap_['M'] = 10u;
        aamap_['N'] = 11u;
        aamap_['P'] = 12u;
        aamap_['Q'] = 13u;
        aamap_['R'] = 14u;
        aamap_['S'] = 15u;
        aamap_['T'] = 16u;
        aamap_['V'] = 17u;
        aamap_['W'] = 18u;
        aamap_['Y'] = 19u;
        aamap_['X'] = 20u;
        aamap_['B'] = aamap_['D'];
        aamap_['J'] = aamap_['L'];
        aamap_['O'] = aamap_['X'];
        aamap_['U'] = aamap_['X'];
        aamap_['Z'] = aamap_['E'];
        return aamap_;
    }

    constexpr auto amino_acid_map = make_amino_acid_map();

    constexpr int aa_to_num(char aa)
    {
        return (aa >= 'A' && aa <= 'Z') ? amino_acid_map[aa] : 0;
    }


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
        size_t bc_{1u};
        std::unordered_map<char, uint8_t> alpha_map_{};
        std::vector<std::vector<uint64_t>> dgram_buffer_;

    public:
        bitvector hits_{};
        seqan::hibf::interleaved_bloom_filter::containment_agent_type agent_{};

        // Rule of 5
        DGramIndex() = default;
        DGramIndex(const DGramIndex &) = default;
        DGramIndex & operator=(const DGramIndex &) = default;
        DGramIndex(DGramIndex&&) noexcept = default;
        DGramIndex & operator=(DGramIndex &&) noexcept = default;
        ~DGramIndex() = default;

        explicit DGramIndex(size_t min_gap, size_t max_gap, size_t pad, size_t hc, float fpr, std::vector<std::string> bins) :
                min_gap_(min_gap),
                max_gap_(max_gap),
                pad_(pad),
                hc_(hc),
                fpr_(fpr),
                bins_(bins),
                bc_(bins_.size()),
                hits_(bc_, true)
        {
            init_alphabet();
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

        // void process_sequence(const std::string& seq, const size_t bin_id)
        // {
        //     // Need at least 2 residues on each side
        //     if (seq.size() < min_gap_ + 5) return;

        //     for (size_t i = 1; i + min_gap_ + 2 < seq.size(); ++i)
        //     {
        //         char a1 = seq[i - 1];
        //         char a2 = seq[i];
        //         if (alpha_map_.find(a1) == alpha_map_.end()) continue;
        //         if (alpha_map_.find(a2) == alpha_map_.end()) continue;

        //         for (size_t gap = min_gap_; gap <= max_gap_; ++gap)
        //         {
        //             size_t j = i + gap + 1;
        //             if (j + 1 >= seq.size()) break;

        //             char b1 = seq[j];
        //             char b2 = seq[j + 1];
        //             if (alpha_map_.find(b1) == alpha_map_.end()) continue;
        //             if (alpha_map_.find(b2) == alpha_map_.end()) continue;

        //             uint8_t code_a1 = alpha_map_[a1];
        //             uint8_t code_a2 = alpha_map_[a2];
        //             uint8_t code_b1 = alpha_map_[b1];
        //             uint8_t code_b2 = alpha_map_[b2];

        //             uint64_t code = static_cast<uint64_t>(gap);
        //             code = code * 160000ULL + static_cast<uint64_t>(code_a1) * 8000ULL
        //                             + static_cast<uint64_t>(code_a2) * 400ULL
        //                             + static_cast<uint64_t>(code_b1) * 20ULL
        //                             + static_cast<uint64_t>(code_b2);

        //             dgram_buffer_[bin_id].emplace_back(code);
        //         }
        //     }
        // }

        void process_sequence(const std::string& seq, const size_t bin_id)
        {
            // Need at least 3 residues on each side of the gap
            if (seq.size() < min_gap_ + 7) return;

            for (size_t i = 2; i + min_gap_ + 3 < seq.size(); ++i)
            {
                char a1 = seq[i - 2];
                char a2 = seq[i - 1];
                char a3 = seq[i];

                // Skip if any left-side residues are invalid
                if (alpha_map_.find(a1) == alpha_map_.end()) continue;
                if (alpha_map_.find(a2) == alpha_map_.end()) continue;
                if (alpha_map_.find(a3) == alpha_map_.end()) continue;

                for(size_t gap = min_gap_; gap <= max_gap_; ++gap)
                {
                    size_t j = i + gap + 1;
                    if (j + 2 >= seq.size()) break; // Need 3 residues on right

                    char b1 = seq[j];
                    char b2 = seq[j + 1];
                    char b3 = seq[j + 2];

                    // Skip if any right-side residues are invalid
                    if (alpha_map_.find(b1) == alpha_map_.end()) continue;
                    if (alpha_map_.find(b2) == alpha_map_.end()) continue;
                    if (alpha_map_.find(b3) == alpha_map_.end()) continue;

                    uint8_t code_a1 = alpha_map_[a1];
                    uint8_t code_a2 = alpha_map_[a2];
                    uint8_t code_a3 = alpha_map_[a3];
                    uint8_t code_b1 = alpha_map_[b1];
                    uint8_t code_b2 = alpha_map_[b2];
                    uint8_t code_b3 = alpha_map_[b3];

                    uint64_t code = static_cast<uint64_t>(gap);
                    code = code * 64000000ULL  // 20^6
                        + static_cast<uint64_t>(code_a1) * 3200000ULL // 20^5
                        + static_cast<uint64_t>(code_a2) * 160000ULL  // 20^4
                        + static_cast<uint64_t>(code_a3) * 8000ULL    // 20^3
                        + static_cast<uint64_t>(code_b1) * 400ULL     // 20^2
                        + static_cast<uint64_t>(code_b2) * 20ULL      // 20^1
                        + static_cast<uint64_t>(code_b3);            // 20^0

                    dgram_buffer_[bin_id].emplace_back(code);
                }
            }
        }


        size_t getMinGap() const
        {
            return min_gap_;
        }

        size_t getMaxGap() const
        {
            return max_gap_;
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

        void spawn_agent()
        {
            agent_ = dibf_.containment_agent();
        }

        bitvector & query(uint64_t const dgram)
        {
            hits_ = agent_.bulk_contains(dgram);
            return hits_;
        }

        void dump_info() const
        {
            DBG(min_gap_);
            DBG(max_gap_);
            DBG(pad_);
            DBG(hc_);
            DBG(fpr_);
            DBG(bc_);
        }

        size_t getMinGap()
        {
            return min_gap_;
        }

        size_t getMaxGap()
        {
            return max_gap_;
        }

        template<class Archive>
        void serialize(Archive &archive)
        {
            archive(min_gap_, max_gap_, pad_, hc_, fpr_, bins_, dibf_, bc_, hits_);
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
