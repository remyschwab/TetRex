#include <array>
#include <string>
#include <iostream>

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
    constexpr const char* amino_acids = "ACDEFGHIKLMNPQRSTVWY";

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
        int simd_identity(const std::vector<uint8_t>& a, const std::vector<uint8_t>& b) {
            size_t length = a.size();
            if (b.size() != length) {
                throw std::runtime_error("Vectors must be the same length.");
            }
    
            const uint8_t* pa = a.data();
            const uint8_t* pb = b.data();
    
            int matches = 0;
            size_t i = 0;
    
            for (; i + 31 < length; i += 32) {
                __m256i va = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(pa + i));
                __m256i vb = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(pb + i));
                __m256i cmp = _mm256_cmpeq_epi8(va, vb);
                int mask = _mm256_movemask_epi8(cmp);
                matches += __builtin_popcount(mask);
            }
    
            for (; i < length; ++i)
                if (pa[i] == pb[i]) ++matches;
    
            return matches;
        }
    
    #elif defined(__ARM_NEON) || defined(__aarch64__)
        #include <arm_neon.h>
        int simd_identity(const std::vector<uint8_t>& a, const std::vector<uint8_t>& b) {
            size_t length = a.size();
            if (b.size() != length) {
                throw std::runtime_error("Vectors must be the same length.");
            }
    
            const uint8_t* pa = a.data();
            const uint8_t* pb = b.data();
    
            int matches = 0;
            size_t i = 0;
    
            for (; i + 15 < length; i += 16) {
                uint8x16_t va = vld1q_u8(pa + i);
                uint8x16_t vb = vld1q_u8(pb + i);
                uint8x16_t cmp = vceqq_u8(va, vb);
    
                // Accumulate matches (each match is 0xFF)
                uint64x2_t wide = vpaddlq_u32(vpaddlq_u16(vpaddlq_u8(cmp)));
                matches += vgetq_lane_u64(wide, 0) / 255;
                matches += vgetq_lane_u64(wide, 1) / 255;
            }
    
            for (; i < length; ++i)
                if (pa[i] == pb[i]) ++matches;
    
            return matches;
        }
    
    #else
        // Scalar fallback
        int simd_identity(const std::vector<uint8_t>& a, const std::vector<uint8_t>& b) {
            size_t length = a.size();
            if (b.size() != length) {
                throw std::runtime_error("Vectors must be the same length.");
            }
    
            int matches = 0;
            for (size_t i = 0; i < length; ++i)
                if (a[i] == b[i]) ++matches;
    
            return matches;
        }
    #endif
    
    std::array<std::string, 200> canonical_numbering = {"1", "2",
        "3", "3A",
        "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18",
        "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31",
        "32", "32A", "32B", "33C", "33B", "33A", "33",
        "34", "35", "36", "37", "38", "39",
        "40", "40A",
        "41", "42", "43",
        "44", "44A",
        "45", "45A",
        "46", "46A",
        "47", "47A",
        "48", "48A", "48B", 
        "49", "49A",
        "50",
        "51", "51A",
        "52", "53", "54", "55", "56", "57", "58", "59",
        "60", "60A", "60B", "60C", "60D", "61E", "61D", "61C", "61B", "61A", "61",
        "62", "63", "64", "65", "66",
        "67", "67A", "67B", 
        "68", "68A", "68B", 
        "69", "69A", "69B",
        "70", 
        "71", "71A", "71B", 
        "72", 
        "73", "73A", "73B",
        "74", "75", "76", "77", "78", "79", 
        "80", "80A", 
        "81", "81A", "81B", "81C",
        "82", "82A", 
        "83", "83A", "83B",
        "84", 
        "85", "85A", "85B", "85C", "85D", 
        "86", "86A", 
        "87", "88", "89", "90", "91", "92", "93", "94", "95", 
        "96", "96A",
        "97", "98", "99", "100", "101", "102", "103", "104", "105", "106", "107", "108", 
        "109", "110",
        "111", "111A", "111B", "111C", "111D", "111E", "111F", "111G", "111H", "111I", "111J", 
        "111K", "111L", "112L",
        "112K", "112J", "112I", "112H", "112G", "112F", "112E", "112D", "112C", "112B", "112A", "112 ",
        "113","114","115","116","117","118",
        "119", "119A",
        "120","121","122","123","124","125", "126","127","128"
    };

    robin_hood::unordered_map<std::string, uint8_t> make_numbering_map()
    {
        robin_hood::unordered_map<std::string, uint8_t> numbering_map;
        for(size_t i = 0; i < 200; ++i) numbering_map[canonical_numbering[i]] = i;
        return numbering_map;
    }

    std::string compute_sequence_from_alignment(const std::vector<uint8_t> &cnn_alignment)
    {
        std::string abseq;
        for(uint8_t residue_byte: cnn_alignment)
        {
            if(residue_byte == 0) residue_byte = 45;
            abseq += residue_byte;
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
    using numbering_map_t = robin_hood::unordered_map<std::string, uint8_t>;

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


    void parse_anarci_output(cnn_alignment_t &query_aln, const numbering_map_t &numbering_map)
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
                query_aln[(numbering_map.find(numbered_res.position))->second] = numbered_res.residue;
            }
            else if(iss >> numbered_res.chain >> numbered_res.position >> numbered_res.insertion_code >> numbered_res.residue)
            {
                std::string insertion_id = numbered_res.position+numbered_res.insertion_code;
                query_aln[(numbering_map.find(insertion_id))->second] = numbered_res.residue;
            }
            else
            {
                std::cerr << "Failed to parse numbering: " << line << std::endl;
            }
        }
    }

    void drive_antibody_query(const antibody_query_arguments &cmd_args)
    {
        auto numbering_map = make_numbering_map();
        cnn_alignment_t query_alignment(200);
        if(cmd_args.anarci_output == "-")
        {
            parse_anarci_output(query_alignment, numbering_map);
        }
        // seqan3::debug_stream << query_alignment << std::endl;

        Index idx_struct;
        load_idx(idx_struct, cmd_args.idx);
        auto agent = idx_struct.hibf_.membership_agent();
        
        kmer_vec_t query_kmers;
        decompose_query(query_kmers, query_alignment);
        // seqan3::debug_stream << query_kmers.size() << std::endl; // 172 kmers
        
        auto & result1 = agent.membership_for(query_kmers, query_kmers.size());
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
                int identity = simd_identity(query_alignment, cn_aln)/2;
                std::string abv = compute_sequence_from_alignment(cn_aln);
                std::cout << abv << "\t" << identity << std::endl;
            }
        }
    }
}
