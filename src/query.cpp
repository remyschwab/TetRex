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
    // TODO: What if k and query_length are the same?
    double km_probability = compute_k_probability(k)*m;
    double running_probability = km_probability;
    double prefix_probability;
    double multi = static_cast<double>(multiplyer);
    seqan3::debug_stream << "\nQLENGTH: " << query_length << std::endl;
    seqan3::debug_stream << "KSIZE: " << k << std::endl;
    seqan3::debug_stream << "TEXTLENGTH: " << m << std::endl;
    seqan3::debug_stream << "MULTIPLYER: " << multiplyer << std::endl;
    seqan3::debug_stream << "KM_Pr: " << km_probability << std::endl;
    for(size_t prefix_length = k+1; prefix_length <= query_length; prefix_length++)
    {
        prefix_probability = compute_k_probability(prefix_length)*m;
        running_probability = running_probability*(multi*km_probability) + (multi*prefix_probability);
    }
    return running_probability;
}


void standard_lib_search(std::filesystem::path matchpath, std::string query, const int &text_length)
{
    // Standard library RegEx Search
    int hitsNr = 0;
    std::fstream file_std;
    file_std.open("standard_search.txt", std::ios::out);
    file_std.clear();
    std::regex regular_expression(query);
    hitsNr = matches(stream_as_string(matchpath), regular_expression, file_std);
    seqan3::debug_stream << "ACTUAL INSTANCE RATE: " << static_cast<double>(hitsNr)/text_length << std::endl;
    file_std.close();
}


void drive_query(const query_arguments &cmd_args)
{
    // Evaluate and search for Regular Expression
    seqan3::debug_stream << "Querying:" << std::endl;
    uint8_t qlength = cmd_args.k;
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
    
    std::vector<std::vector<std::string>> matrix{};
    for(auto i : knfa)
    {
        dfs_na(i, matrix);
    }

    seqan3::debug_stream << "DONE" << std::endl;

    /////////// MODELING STEP ///////////////
    size_t query_length = matrix[0].size()+qlength-1;
    int text_length = cmd_args.text_length;
    size_t multiplyer = matrix.size();
    double result = compute_knut_model(query_length, qlength, text_length, multiplyer);
    seqan3::debug_stream << "FINAL PROBABILITY: " << result << std::endl;
    /////////// MODELING STEP ///////////////

    standard_lib_search(cmd_args.acid_lib, cmd_args.regex, text_length);
}
