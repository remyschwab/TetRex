//
// Created by Remy Schwab on 20.09.22.
//

#include "query.h"
#include <omp.h>


bitvector query_ibf(uint32_t &bin_count, robin_hood::unordered_map<uint64_t, bitvector> &hash_to_bits, std::vector<std::pair<std::string, uint64_t>> &path)
{   
    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed>::membership_agent::binning_bitvector hit_vector{bin_count};
    std::fill(hit_vector.begin(), hit_vector.end(), true);
    for (auto && kmer : path)
    {
        auto & result = hash_to_bits[kmer.second];
        hit_vector.raw_data() &= result.raw_data();
    }
    return hit_vector;
}

// Helper Function to compute the probability of any kmer for a given k
double compute_k_probability(const uint8_t &k)
{
    return pow(0.25, k);
}


double compute_knut_model(const size_t &query_length, const uint8_t &k, const int &m, const size_t &multiplyer)
{
    double km_probability = compute_k_probability(k)*m;
    double running_probability = km_probability;
    double prefix_probability;
    double multi = static_cast<double>(multiplyer);
    for(size_t prefix_length = k+1; prefix_length <= query_length; prefix_length++)
    {
        prefix_probability = compute_k_probability(prefix_length)*m;
        running_probability = running_probability*(multi*km_probability) + (multi*prefix_probability);
    }
    return running_probability;
}


bitvector drive_query(const query_arguments &cmd_args)
{
    
    // Load index from disk
    seqan3::debug_stream << "Reading Index from Disk... ";
    IndexStructure ibf;
    load_ibf(ibf, cmd_args.idx);
    seqan3::debug_stream << "DONE" << std::endl;

    auto bin_count = ibf.getBinCount();

    // Evaluate and search for Regular Expression
    seqan3::debug_stream << "Querying:" << std::endl;
    uint8_t qlength = ibf.k_;
    std::string query = cmd_args.query;
    std::vector<char> a = getAlphabet(query);

    // Postfix to Thompson NFA
    seqan3::debug_stream << "   - Constructing Thompson NFA from RegEx... ";
    State* nfa = post2nfaE(query);
    seqan3::debug_stream << "DONE" << std::endl;

    // Thompson NFA to Korotkov NFA
    seqan3::debug_stream << "   - Construction kNFA from Thompson NFA... ";
    std::vector<kState *> knfa = nfa2knfa(nfa, qlength);
    seqan3::debug_stream << "DONE" << std::endl;
    deleteGraph(nfa);

    // Create kmer path matrix from kNFA
    seqan3::debug_stream << "   - Computing kmer path matrix from kNFA... ";
    
    // Create auxiliary data structures to avoid redundant kmer lookup
    auto hash_adaptor = seqan3::views::kmer_hash(seqan3::ungapped{qlength});
    robin_hood::unordered_map<uint64_t, bitvector> hash_to_bitvector{};
    // std::vector<bitvector> kmer_bitvex;
    // uint32_t vector_idx = 0;

    // Spawn IBF membership agent in this scope because it is expensive
    auto && ibf_ref = ibf.getIBF();
    auto agent = ibf_ref.membership_agent();
    
    std::vector<std::vector<std::string>> matrix{};
    for(auto i : knfa)
    {
        if(ibf.molecule_ == "na")
        {        
            dfs_na(i, matrix, hash_to_bitvector, agent);
        } else {
            dfs_aa(i, matrix, hash_to_bitvector, agent);
        }
    }
    uMatrix(matrix);
    seqan3::debug_stream << "DONE" << std::endl;

    // Search kmer paths in index
    seqan3::debug_stream << "   - Search kmers in index... ";
    path_vector paths_vector;
    if(ibf.molecule_ == "na")
    {
        extract_matrix_paths<seqan3::dna5>(matrix, paths_vector, hash_adaptor);
    } else {
        extract_matrix_paths<seqan3::aa27>(matrix, paths_vector, hash_adaptor);
    }

    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed>::membership_agent::binning_bitvector hit_vector{ibf.getBinCount()};
    std::fill(hit_vector.begin(), hit_vector.end(), false);

    // MODELING STEP
    int text_length = 32;
    size_t query_length = paths_vector[0].size()+qlength-1;
    size_t multiplyer = paths_vector.size();
    seqan3::debug_stream << "FINAL PROBABILITY: " << compute_knut_model(query_length, qlength, text_length, multiplyer) << std::endl;

    for (auto path : paths_vector)
    {
        auto hits = query_ibf(bin_count, hash_to_bitvector, path);
        hit_vector.raw_data() |= hits.raw_data();
    }
    seqan3::debug_stream << "DONE" << std::endl;
    for(auto && bit: hit_vector)
    {
      std::cout << bit << ',';
    }
    seqan3::debug_stream << "\nWrite dot file... ";
    std::string dotfile = cmd_args.graph;
    dotfile += ".dot";
    printGraph(knfa, dotfile);
    seqan3::debug_stream << "DONE" << std::endl;
    return hit_vector;
}
