#include <array>
#include <string>
#include <iostream>
#include <ranges>

#include "seqan3/core/debug_stream/all.hpp"
#include "arg_parse.h"
#include "robin_hood.h"
#include "cnpy.h"
#include <hibf/config.hpp>
#include "hibf/hierarchical_interleaved_bloom_filter.hpp"

namespace aa {

constexpr unsigned char INVALID_AA = 255;

constexpr std::array<unsigned char, UCHAR_MAX + 1> make_aa_lookup_table() {
    std::array<unsigned char, UCHAR_MAX + 1> table{};

    // Fill with INVALID_AA
    for (size_t i = 1; i <= UCHAR_MAX; ++i) {
        table[i] = INVALID_AA;
    }

    table[0] = 27; // Maybe change this? Probably doesn't matter tho

    // Amino acids in standard order     01234567890123456789
    constexpr const char* amino_acids = "ACDEFGHIKLMNPQRSTVWY"; // No reduction
    // constexpr const char* amino_acids = "ACDEFGHIKLMNPQRSTVWY"; // Li reduction
    // constexpr const char* amino_acids = "ACDEFGHIKLMNPQRSTVWY"; // Murphy reduction

    for (unsigned char i = 0; i < 20; ++i) {
        char aa = amino_acids[i];
        table[static_cast<unsigned char>(aa)] = i;
        table[static_cast<unsigned char>(aa + 32)] = i; // lowercase (ASCII only)
    }

    return table;
}

constexpr auto AA_LOOKUP = make_aa_lookup_table();

} // namespace aa



namespace Antibody_Utilities
{

    struct NumberedResidue
    {
        std::string position;          // IMGT position (e.g., 30, 52)
        std::string insertion_code{""};
        unsigned char residue;         // Amino acid
        std::string region;            // CDR1, FR1, etc. (optional for now)
        std::string chain;             // Heavy or Light
    };

    struct QueryMeta
    {
        char chain;
        std::string scheme;
    };

    struct Index
    {
        std::vector<std::string> bins_;
        seqan::hibf::hierarchical_interleaved_bloom_filter hibf_;
        template<class Archive> void serialize(Archive &archive)
        {
            archive(bins_, hibf_);
        }
    };

    template <typename Index>
    void store_idx(Index const &ibf, std::filesystem::path const &opath)
    {
        std::ofstream os{opath, std::ios::binary};
        cereal::BinaryOutputArchive oarchive{os};
        oarchive(ibf);
    }

    template <typename Index>
    void load_idx(Index &ibf, const std::filesystem::path &ipath)
    {
        std::ifstream is{ipath, std::ios::binary};
        cereal::BinaryInputArchive iarchive{is};
        iarchive(ibf);
    }
    
    #if defined(__AVX2__) && defined(__x86_64__)
        #include <immintrin.h>
        #include <vector>
        #include <cstdint>
        #include <stdexcept>
        
        float compute_sequence_identity(const std::vector<uint8_t>& seq1, const std::vector<uint8_t>& seq2)
        {
            if (seq1.size() != seq2.size()) {
                throw std::invalid_argument("Sequences must be of equal length.");
            }
        
            size_t matched = 0;
            size_t valid_positions = 0;
            size_t i = 0;
            size_t len = seq1.size();
        
            const __m256i zero = _mm256_set1_epi8(0);
            const __m256i skip = _mm256_set1_epi8(124);
        
            for (; i + 31 < len; i += 32) {
                __m256i v1 = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(&seq1[i]));
                __m256i v2 = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(&seq2[i]));
        
                // Compare for equality
                __m256i eq = _mm256_cmpeq_epi8(v1, v2);
        
                // Valid mask: (v1 != 0) & (v1 != 124) & (v2 != 0) & (v2 != 124)
                __m256i v1_neq_zero = _mm256_cmpgt_epi8(_mm256_xor_si256(v1, zero), zero); // v1 != 0
                __m256i v1_neq_124  = _mm256_cmpgt_epi8(_mm256_xor_si256(v1, skip), zero); // v1 != 124
                __m256i v2_neq_zero = _mm256_cmpgt_epi8(_mm256_xor_si256(v2, zero), zero);
                __m256i v2_neq_124  = _mm256_cmpgt_epi8(_mm256_xor_si256(v2, skip), zero);
        
                __m256i valid = _mm256_and_si256(_mm256_and_si256(v1_neq_zero, v1_neq_124),
                                                _mm256_and_si256(v2_neq_zero, v2_neq_124));
        
                // Match only where valid
                __m256i match = _mm256_and_si256(eq, valid);
        
                // Count bits
                uint32_t match_mask = _mm256_movemask_epi8(match);
                uint32_t valid_mask = _mm256_movemask_epi8(valid);
        
                matched += __builtin_popcount(match_mask);
                valid_positions += __builtin_popcount(valid_mask);
            }
        
            // Handle remaining elements
            for (; i < len; ++i) {
                uint8_t a = seq1[i], b = seq2[i];
                if ((a != 0 && a != 124) && (b != 0 && b != 124)) {
                    ++valid_positions;
                    if (a == b) ++matched;
                }
            }
        
            return valid_positions > 0 ? static_cast<float>(matched) / valid_positions : 0.0f;
        }
    
    
    #elif defined(__ARM_NEON) || defined(__aarch64__)
        #include <arm_neon.h>
        #include <vector>
        #include <cstdint>
        #include <stdexcept>
        
        float compute_sequence_identity(const std::vector<uint8_t>& seq1, const std::vector<uint8_t>& seq2)
        {
            if (seq1.size() != seq2.size()) {
                throw std::invalid_argument("Sequences must be of equal length.");
            }
        
            size_t matched = 0;
            size_t valid_positions = 0;
            size_t i = 0;
            size_t len = seq1.size();
        
            const uint8x16_t zero = vdupq_n_u8(0);
            const uint8x16_t skip = vdupq_n_u8(124);
        
            for (; i + 15 < len; i += 16) {
                uint8x16_t v1 = vld1q_u8(&seq1[i]);
                uint8x16_t v2 = vld1q_u8(&seq2[i]);
        
                uint8x16_t eq = vceqq_u8(v1, v2);
        
                uint8x16_t v1_valid = vandq_u8(vcgtq_u8(veorq_u8(v1, zero), zero),
                                            vcgtq_u8(veorq_u8(v1, skip), zero));
                uint8x16_t v2_valid = vandq_u8(vcgtq_u8(veorq_u8(v2, zero), zero),
                                            vcgtq_u8(veorq_u8(v2, skip), zero));
                uint8x16_t valid = vandq_u8(v1_valid, v2_valid);
        
                uint8x16_t match = vandq_u8(eq, valid);
        
                // Count set bits
                for (int j = 0; j < 16; ++j) {
                    if (valid[j]) {
                        ++valid_positions;
                        if (match[j]) ++matched;
                    }
                }
            }
        
            for (; i < len; ++i) {
                uint8_t a = seq1[i], b = seq2[i];
                if ((a != 0 && a != 124) && (b != 0 && b != 124)) {
                    ++valid_positions;
                    if (a == b) ++matched;
                }
            }
        
            return valid_positions > 0 ? static_cast<float>(matched) / valid_positions : 0.0f;
        }
    
    
    #else
        float compute_sequence_identity(const std::vector<uint8_t>& seq1,
                                        const std::vector<uint8_t>& seq2) {
            if (seq1.size() != seq2.size()) {
                throw std::invalid_argument("Sequences must be of equal length.");
            }
        
            size_t matched = 0;
            size_t valid_positions = 0;
        
            for (size_t i = 0; i < seq1.size(); ++i) {
                uint8_t a = seq1[i];
                uint8_t b = seq2[i];
        
                // Ignore positions with 0 or 124 in either sequence
                if ((a != 0 && a != 124) && (b != 0 && b != 124)) {
                    ++valid_positions;
                    if (a == b) {
                        ++matched;
                    }
                }
            }
        
            if (valid_positions == 0) {
                return 0.0f; // Or NaN if you want to signal invalid comparison
            }
        
            return static_cast<float>(matched) / static_cast<float>(valid_positions);
        }
    
    #endif
    
    static constexpr std::array<std::pair<std::string_view, uint8_t>, 200> canonical_pairs
    {{
        {"1", 0}, {"2", 1}, {"3", 2}, {"3A", 3}, {"4", 4}, {"5", 5}, {"6", 6}, {"7", 7}, {"8", 8}, {"9", 9},
        {"10", 10}, {"11", 11}, {"12", 12}, {"13", 13}, {"14", 14}, {"15", 15}, {"16", 16}, {"17", 17}, {"18", 18},
        {"19", 19}, {"20", 20}, {"21", 21}, {"22", 22}, {"23", 23}, {"24", 24}, {"25", 25}, {"26", 26}, {"27", 27},
        {"28", 28}, {"29", 29}, {"30", 30}, {"31", 31}, {"32", 32}, {"32A", 33}, {"32B", 34}, {"33C", 35},
        {"33B", 36}, {"33A", 37}, {"33", 38}, {"34", 39}, {"35", 40}, {"36", 41}, {"37", 42}, {"38", 43},
        {"39", 44}, {"40", 45}, {"40A", 46}, {"41", 47}, {"42", 48}, {"43", 49}, {"44", 50}, {"44A", 51},
        {"45", 52}, {"45A", 53}, {"46", 54}, {"46A", 55}, {"47", 56}, {"47A", 57}, {"48", 58}, {"48A", 59},
        {"48B", 60}, {"49", 61}, {"49A", 62}, {"50", 63}, {"51", 64}, {"51A", 65}, {"52", 66}, {"53", 67},
        {"54", 68}, {"55", 69}, {"56", 70}, {"57", 71}, {"58", 72}, {"59", 73}, {"60", 74}, {"60A", 75},
        {"60B", 76}, {"60C", 77}, {"60D", 78}, {"61E", 79}, {"61D", 80}, {"61C", 81}, {"61B", 82}, {"61A", 83},
        {"61", 84}, {"62", 85}, {"63", 86}, {"64", 87}, {"65", 88}, {"66", 89}, {"67", 90}, {"67A", 91},
        {"67B", 92}, {"68", 93}, {"68A", 94}, {"68B", 95}, {"69", 96}, {"69A", 97}, {"69B", 98}, {"70", 99},
        {"71", 100}, {"71A", 101}, {"71B", 102}, {"72", 103}, {"73", 104}, {"73A", 105}, {"73B", 106},
        {"74", 107}, {"75", 108}, {"76", 109}, {"77", 110}, {"78", 111}, {"79", 112}, {"80", 113}, {"80A", 114},
        {"81", 115}, {"81A", 116}, {"81B", 117}, {"81C", 118}, {"82", 119}, {"82A", 120}, {"83", 121},
        {"83A", 122}, {"83B", 123}, {"84", 124}, {"85", 125}, {"85A", 126}, {"85B", 127}, {"85C", 128},
        {"85D", 129}, {"86", 130}, {"86A", 131}, {"87", 132}, {"88", 133}, {"89", 134}, {"90", 135},
        {"91", 136}, {"92", 137}, {"93", 138}, {"94", 139}, {"95", 140}, {"96", 141}, {"96A", 142},
        {"97", 143}, {"98", 144}, {"99", 145}, {"100", 146}, {"101", 147}, {"102", 148}, {"103", 149},
        {"104", 150}, {"105", 151}, {"106", 152}, {"107", 153}, {"108", 154}, {"109", 155}, {"110", 156},
        {"111", 157}, {"111A", 158}, {"111B", 159}, {"111C", 160}, {"111D", 161}, {"111E", 162}, {"111F", 163},
        {"111G", 164}, {"111H", 165}, {"111I", 166}, {"111J", 167}, {"111K", 168}, {"111L", 169}, {"112L", 170},
        {"112K", 171}, {"112J", 172}, {"112I", 173}, {"112H", 174}, {"112G", 175}, {"112F", 176}, {"112E", 177},
        {"112D", 178}, {"112C", 179}, {"112B", 180}, {"112A", 181}, {"112 ", 182}, {"113", 183}, {"114", 184},
        {"115", 185}, {"116", 186}, {"117", 187}, {"118", 188}, {"119", 189}, {"119A", 190}, {"120", 191},
        {"121", 192}, {"122", 193}, {"123", 194}, {"124", 195}, {"125", 196}, {"126", 197}, {"127", 198}, {"128", 199}
    }};

    const robin_hood::unordered_flat_map<std::string_view, uint8_t> canonical_map = [] {
        robin_hood::unordered_flat_map<std::string_view, uint8_t> map;
        for (const auto& [key, val] : canonical_pairs) {
            map[key] = val;
        }
        return map;
    }();

    std::string compute_sequence_from_alignment(const std::vector<uint8_t> &cnn_alignment, const std::vector<uint8_t> &query_alignment)
    {
        std::string abseq;
        for(auto [ref_res, query_res]: std::views::zip(cnn_alignment, query_alignment))
        {
            if(ref_res == 0)
            {
                if(query_res == 0) continue;
                abseq += static_cast<unsigned char>(ref_res);
            }
            abseq += ref_res;
        }
        return abseq;
    }
    
    struct CanonicalAlignment
    {
        std::array<uint8_t, 200> alignment_data;
    };
    
    std::array<uint8_t, 200> region_codes = {
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,       // FR1
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0,
        1, 1, 1, 1, 1, 1, 1, 1, 1,          // CDR1
        1, 1, 1, 1, 1, 1, 1, 1,
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2,       // FR2
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
        2, 2, 2, 2, 2, 2,
        3, 3, 3, 3, 3, 3, 3, 3, 3, 3,       // CDR2
        3, 3, 3, 3, 3, 3, 3, 3, 3,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4,       // FR3
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4,
        5, 5, 5, 5, 5, 5, 5, 5, 5, 5,       // CDR3
        5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
        5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
        5, 5, 5, 5, 5, 5, 5,
        6, 6, 6, 6, 6, 6, 6, 6, 6, 6,       // FR4
        6, 6
    };

    void print_bits(const uint64_t &kmer)
    {
        seqan3::debug_stream << std::bitset<64>(kmer) << std::endl;
    }

    uint64_t bitmask = 8388607; // 3 bits at most significant end for region encoding and 5x5 for encoding the kmer
    uint64_t preshift = 33554400;

    using numberings_t = std::vector<std::vector<uint8_t>>;
    using cnn_alignment_t = std::vector<uint8_t>;
    using alignment_it_t = cnn_alignment_t::const_iterator;
    using kmer_vec_t = std::vector<uint64_t>;

    cnn_alignment_t test_query = {81, 86, 81, 0, 76, 81, 81, 87, 71, 65,  0, 71, 76, 76, 75, 80, 83,
        69, 84, 76, 83, 76, 84, 67, 65, 86, 89, 71, 71, 83, 70,  0,  0,  0,
         0,  0,  0,  0,  0,  0, 83, 71, 89, 89, 87, 83,  0, 87, 73, 82, 81,
         0, 80,  0, 80,  0, 71,  0, 75,  0,  0, 71,  0, 76, 69,  0, 87, 73,
        71, 69, 73, 78, 72, 83,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
         0, 71, 83, 84, 78, 89,  0,  0, 78,  0,  0, 80,  0,  0, 83, 76,  0,
         0, 75,  0,  0,  0, 83, 82, 86, 84, 73, 83, 86,  0, 68,  0,  0,  0,
        84,  0, 83,  0,  0, 75, 78,  0,  0,  0,  0, 81,  0, 70, 83, 76, 75,
        76, 83, 83, 86, 84, 65,  0, 65, 68, 84, 65, 86, 89, 89, 67, 65, 82,
        71, 71, 71, 80, 86, 80, 65,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
         0,  0,  0,  0,  0,  0,  0,  0,  0, 65, 71, 83, 82, 71, 80, 73, 68,
        89, 87, 71,  0, 81, 71, 84, 76, 86, 84, 86, 83, 83};

    void decompose_alignment(alignment_it_t begin, alignment_it_t end, kmer_vec_t &kmer_vec, const uint8_t region_code)
    {
        uint64_t kmer = 0;
        uint8_t k = 5;
        kmer = (kmer<<k) | aa::AA_LOOKUP[*(begin+0)];
        kmer = (kmer<<k) | aa::AA_LOOKUP[*(begin+1)];
        kmer = (kmer<<k) | aa::AA_LOOKUP[*(begin+2)];
        kmer = (kmer<<k) | aa::AA_LOOKUP[*(begin+3)];
        kmer = (kmer<<k) | aa::AA_LOOKUP[*(begin+4)];
        kmer_vec.push_back(kmer);
        auto it = begin+k;
        while(it < end)
        {
            kmer = (kmer<<k) & preshift; // Create room for new residue on least significant end
            kmer = kmer | aa::AA_LOOKUP[*(it)]; // Add the new kmer
            kmer |= region_code;
            kmer_vec.push_back(kmer);
            ++it;
        }
    }

    void decompose_alignments(kmer_vec_t &bin_kmers, const numberings_t &numberings)
    {
        for(const cnn_alignment_t &alignment: numberings)
        {
            decompose_alignment(alignment.begin() + 0u,  alignment.begin() + 27u,  bin_kmers, 0u);  // Decompose 0s
            decompose_alignment(alignment.begin() + 27u, alignment.begin() + 44u,  bin_kmers, 1u);  // Decompose 1s
            decompose_alignment(alignment.begin() + 44u, alignment.begin() + 70u,  bin_kmers, 2u);  // Decompose 2s
            decompose_alignment(alignment.begin() + 70u, alignment.begin() + 89u,  bin_kmers, 3u);  // Decompose 3s
            decompose_alignment(alignment.begin() + 89u, alignment.begin() + 151u, bin_kmers, 4u);  // Decompose 4s
            decompose_alignment(alignment.begin() + 151u, alignment.begin() + 188u, bin_kmers, 5u); // Decompose 5s
            decompose_alignment(alignment.begin() + 188u, alignment.begin() + 200u, bin_kmers, 6u); // Decompose 6s
        }
    }

    void decompose_query(kmer_vec_t &bin_kmers, const cnn_alignment_t &alignment)
    {
        decompose_alignment(alignment.begin() + 0u,  alignment.begin() + 27u,  bin_kmers, 0u);  // Decompose 0s
        decompose_alignment(alignment.begin() + 27u, alignment.begin() + 44u,  bin_kmers, 1u);  // Decompose 1s
        decompose_alignment(alignment.begin() + 44u, alignment.begin() + 70u,  bin_kmers, 2u);  // Decompose 2s
        decompose_alignment(alignment.begin() + 70u, alignment.begin() + 89u,  bin_kmers, 3u);  // Decompose 3s
        decompose_alignment(alignment.begin() + 89u, alignment.begin() + 151u, bin_kmers, 4u);  // Decompose 4s
        decompose_alignment(alignment.begin() + 151u, alignment.begin() +188u, bin_kmers, 5u);  // Decompose 5s
        decompose_alignment(alignment.begin() + 188u, alignment.begin() +200u, bin_kmers, 6u);  // Decompose 6s
    }

    void drive_antibody_index(const antibody_index_arguments &cmd_args)
    {
        std::vector<std::string> input_npz_paths;
        for(auto path: cmd_args.canonical_alignments)
        {
            input_npz_paths.push_back(std::filesystem::absolute(path));
        }
        seqan3::debug_stream << "[FOUND FILES]" << std::endl;
        std::vector<kmer_vec_t> all_kmers;
        for(auto &path: input_npz_paths)
        {
            kmer_vec_t bin_kmers;

            cnpy::npz_t npz_data = cnpy::npz_load(path); // Keys are idxs and numberings
            cnpy::NpyArray alignments = npz_data["numberings"];
            // Sanity check
            assert(alignments.shape.size() == 2);
            size_t num_rows = alignments.shape[0];
            size_t num_cols = alignments.shape[1];
    
            // Get raw pointer to data
            uint8_t* raw_data = alignments.data<uint8_t>();
    
            // Reconstruct the nested vector
            numberings_t numberings(num_rows, std::vector<uint8_t>(num_cols));
            for(size_t i = 0; i < num_rows; ++i) std::copy(raw_data + i * num_cols, raw_data + (i + 1) * num_cols, numberings[i].begin());

            // Generate kmers
            decompose_alignments(bin_kmers, numberings);
            all_kmers.push_back(bin_kmers);
        }

        seqan3::debug_stream << "[DONE READING]" << std::endl;

        auto get_user_bin_data = [&](size_t const user_bin_id, seqan::hibf::insert_iterator it)
        {
            for (auto value : all_kmers[user_bin_id])
                it = value;
        };

        seqan::hibf::config config{.input_fn = get_user_bin_data, // required
            .number_of_user_bins = input_npz_paths.size(),     // required
            .number_of_hash_functions = 3u,
            .maximum_fpr = 0.05,
            .threads = 1u};
        
        seqan::hibf::hierarchical_interleaved_bloom_filter hibf{config};

        Index idx_struct{input_npz_paths, hibf};
        seqan3::debug_stream << "[WRITING INDEX TO DISK]" << std::endl;
        store_idx(idx_struct, cmd_args.output_index_path);
    }


    void parse_anarci_output(cnn_alignment_t &query_aln)
    {
        std::string line;
        while(std::getline(std::cin, line))
        {
            if(line.starts_with("#")) continue;
            if(line.starts_with("//")) continue;
            std::istringstream iss(line);
            NumberedResidue numbered_res;
            if(iss >> numbered_res.chain >> numbered_res.position >> numbered_res.residue)
            {
                if(numbered_res.residue == 45) continue; // For '-'s
                query_aln[canonical_map.at(numbered_res.position)] = numbered_res.residue;
            }
            else if(iss >> numbered_res.chain >> numbered_res.position >> numbered_res.insertion_code >> numbered_res.residue)
            {
                std::string insertion_id = numbered_res.position+numbered_res.insertion_code;
                query_aln[canonical_map.at(insertion_id)] = numbered_res.residue;
            }
            else
            {
                std::cerr << "Failed to parse numbering: " << line << std::endl;
            }
        }
    }

    void drive_antibody_query(const antibody_query_arguments &cmd_args)
    {
        std::pair<std::string, float> top_hit = std::make_pair("",0.0);
        cnn_alignment_t query_alignment(200);
        float kmer_fraction = 1;
        if(cmd_args.anarci_output == "-")
        {
            parse_anarci_output(query_alignment);
        }
        // seqan3::debug_stream << query_alignment << std::endl;

        Index idx_struct;
        load_idx(idx_struct, cmd_args.idx);
        auto agent = idx_struct.hibf_.membership_agent();
        
        kmer_vec_t query_kmers;
        decompose_query(query_kmers, query_alignment);
        // seqan3::debug_stream << query_kmers.size() << std::endl; // 172 kmers
        
        size_t abs_kmer_threshold = std::floor(kmer_fraction*query_kmers.size());
        auto & result1 = agent.membership_for(query_kmers, abs_kmer_threshold);
        // seqan3::debug_stream << result1 << std::endl;
        std::cout << "Sequence" << "\t" << "Identity" << std::endl;
        for(auto && bin: result1)
        {
            std::filesystem::path path = idx_struct.bins_[bin];
            cnpy::npz_t npz_data = cnpy::npz_load(path); // Keys are idxs and numberings
            cnpy::NpyArray alignments = npz_data["numberings"];
            // Sanity check
            assert(alignments.shape.size() == 2);
            size_t num_rows = alignments.shape[0];
            size_t num_cols = alignments.shape[1];
    
            // Get raw pointer to data
            uint8_t* raw_data = alignments.data<uint8_t>();
    
            // Reconstruct the nested vector
            numberings_t numberings(num_rows, cnn_alignment_t(num_cols));
            for(size_t i = 0; i < num_rows; ++i) std::copy(raw_data + i * num_cols, raw_data + (i + 1) * num_cols, numberings[i].begin());
            for(const cnn_alignment_t &cn_aln: numberings)
            {
                float identity = compute_sequence_identity(query_alignment, cn_aln);
                if(identity > top_hit.second) top_hit = std::make_pair(compute_sequence_from_alignment(cn_aln, query_alignment), identity);
            }
        }
        std::cout << top_hit.first << "\t" << top_hit.second << std::endl;
    }
}
