#include <array>
#include <string>

#include "seqan3/core/debug_stream/all.hpp"
#include "arg_parse.h"
#include "robin_hood.h"
#include "cnpy.h"
#include <hibf/config.hpp>
#include "hibf/hierarchical_interleaved_bloom_filter.hpp"


namespace Antibody_Utilities
{

    struct NumberedResidue
    {
        int position;        // IMGT position (e.g., 30, 52)
        char insertion_code; // 'A', 'B', etc. for 30A, 30B; or ' ' if none
        char aa;             // Amino acid
        std::string region;  // CDR1, FR1, etc. (optional for now)
    };
    
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

    void set_residue_vals(std::array<unsigned char, UCHAR_MAX+1> &aamap)
    {

    }

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

    cnn_alignment_t test_query = {81, 86, 81,  0, 76, 81, 81, 87, 71, 65,  0, 71, 76, 76, 75, 80, 83,
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
        // Decompose FR1
        uint64_t kmer = 0;
        uint8_t k = 5;
        kmer = (kmer<<5) | cnn_alignment[0]; // Init first kmer
        kmer = (kmer<<5) | cnn_alignment[1]; // Init first kmer
        kmer = (kmer<<5) | cnn_alignment[2]; // Init first kmer
        kmer = (kmer<<5) | cnn_alignment[3]; // Init first kmer
        kmer = (kmer<<5) | cnn_alignment[4]; // Init first kmer
        kmer_vec.push_back(kmer);
        auto it = begin + k;
        while(it < end)
        {
            kmer = (kmer<<5) & preshift; // Create room for new residue on least significant end
            kmer = kmer | cnn_alignment[i]; // Add the new kmer
            kmer |= region_code;
            kmer_vec.push_back(kmer);
        }
    }

    void decompose_alignments(kmer_vec_t &bin_kmers, const numberings_t &numberings)
    {
        for(const cnn_alignment_t &alignment: numberings)
            decompose_alignment(alignment, bin_kmers);
    }

    void drive_antibody_index(const antibody_index_arguments &cmd_args)
    {
        std::vector<std::string> input_npz_paths;
        std::array<unsigned char, UCHAR_MAX+1> aamap = {};
        set_residue_vals(aamap);
        for(auto path: cmd_args.canonical_alignments)
        {
            input_npz_paths.push_back(std::filesystem::absolute(path));
        }
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

        seqan3::debug_stream << "Finished Reading Files" << std::endl;

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

        auto agent = hibf.membership_agent();
        kmer_vec_t query_kmers;
        decompose_alignment(test_query, query_kmers);
        auto & result1 = agent.membership_for(query_kmers, query_kmers.size());
        seqan3::debug_stream << result1.size() << std::endl;
    }
}
