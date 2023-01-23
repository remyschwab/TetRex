//
// Created by Remy Schwab on 20.09.22.
//

#include "query.h"
#include <omp.h>


bitvector query_ibf(uint32_t &bin_count, robin_hood::unordered_map<uint64_t, uint32_t> &hash_to_idx,
  std::vector<bitvector> &kmer_bitvex, std::vector<std::pair<std::string, uint64_t>> &path)
{   
    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed>::membership_agent::binning_bitvector hit_vector{bin_count};
    std::fill(hit_vector.begin(), hit_vector.end(), true);
    for (auto && kmer : path)
    {
        auto & result = kmer_bitvex[hash_to_idx[kmer.second]];
        hit_vector.raw_data() &= result.raw_data();
    }
    return hit_vector;
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
    robin_hood::unordered_map<uint64_t, uint32_t> hash_to_idx{};
    std::vector<bitvector> kmer_bitvex;
    uint32_t vector_idx = 0;
    // Spawn IBF membership agent in this scope because it is expensive
    auto && ibf_ref = ibf.getIBF();
    auto agent = ibf_ref.membership_agent();
    
    std::vector<std::vector<std::string>> matrix{};
    for(auto i : knfa)
    {
        dfs(i, matrix, vector_idx, hash_to_idx, kmer_bitvex, agent);
    }
    uMatrix(matrix);
    seqan3::debug_stream << "DONE" << std::endl;

    // Search kmer paths in index
    seqan3::debug_stream << "   - Search kmers in index... ";
    std::vector<std::vector<std::pair<std::string, uint64_t>>> paths_vector;

    for(auto i : matrix)
    {
        std::vector<std::pair<std::string, uint64_t>> hash_vector;
        for(auto j : i)
        {
            std::vector<seqan3::dna5> acid_vec = convertStringToDNA(j);
            auto digest = acid_vec | hash_adaptor;
            // Create a vector of kmer hashes that correspond
            hash_vector.push_back(std::make_pair(j, digest[0]));
        }
        paths_vector.push_back(hash_vector);
    }

    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed>::membership_agent::binning_bitvector hit_vector{ibf.getBinCount()};
    std::fill(hit_vector.begin(), hit_vector.end(), false);

    int total,eins;
    eins = 0;
    for (auto path : paths_vector)
    {
        auto hits = query_ibf(bin_count, hash_to_idx, kmer_bitvex, path);
        hit_vector.raw_data() |= hits.raw_data();
    }
    seqan3::debug_stream << "DONE" << std::endl;
    for(auto && bit: hit_vector)
    {
      std::cout << bit << ',';
      total++;
      eins+=bit;  
    }
    seqan3::debug_stream << "\nWrite .dot file... ";
    std::string dotfile = cmd_args.graph;
    dotfile += ".dot";
    printGraph(knfa, dotfile);
    seqan3::debug_stream << "DONE" << std::endl;
    return hit_vector;
}
