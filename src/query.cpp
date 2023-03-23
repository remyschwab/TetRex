//
// Created by Remy Schwab on 20.09.22.
//

#include "query.h"


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
    // seqan3::debug_stream << "\nQLENGTH: " << query_length << std::endl;
    // seqan3::debug_stream << "KSIZE: " << k << std::endl;
    // seqan3::debug_stream << "TEXTLENGTH: " << m << std::endl;
    // seqan3::debug_stream << "MULTIPLYER: " << multiplyer << std::endl;
    // seqan3::debug_stream << "KM_Pr: " << km_probability << std::endl;
    for(size_t prefix_length = k+1; prefix_length <= query_length; prefix_length++)
    {
        prefix_probability = compute_k_probability(prefix_length)*m;
        running_probability = running_probability*(multi*km_probability) + (multi*prefix_probability);
    }
    return running_probability;
}


bitvector query_ibf(size_t &bin_count, robin_hood::unordered_map<uint64_t, bitvector> &hash_to_bits, std::vector<std::pair<std::string, uint64_t>> &path)
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


void verify_fasta_hit(const std::filesystem::path &bin_path, re2::RE2 &crx)
{
    seqan3::sequence_file_input<seqan3::sequence_file_input_default_traits_aa> fasta_handle{bin_path};
    std::string match;
    for(auto &[SEQ, ID, QUAL]: fasta_handle)
    {
        re2::StringPiece bin_content(to_string(SEQ));
        while (RE2::FindAndConsume(&bin_content, crx, &match))
        {
            std::cout << ">" << ID << "\t" << match << std::endl;
        }
    }
}


void single_disk_search(const bitvector &hits, std::string &query, IndexStructure &ibf)
{
    size_t bins = hits.size();
    const std::filesystem::path lib_path = ibf.acid_libs_[0];
    std::string bin_text = stream_as_string(lib_path);
    for(size_t i = 0; i < bins; i++)
    {
        if(hits[i])
        {
            // RE2::PartialMatch(bin_text, query);
        }
    }
}


void iter_disk_search(const bitvector &hits, std::string &query, IndexStructure &ibf)
{
    size_t bins = hits.size();
    std::filesystem::path lib_path;

    re2::RE2 compiled_regex(query);
    assert(compiled_regex.ok());
    for(size_t i = 0; i < bins; i++)
    {
        if(hits[i])
        {
            lib_path = ibf.acid_libs_[i];
            verify_fasta_hit(lib_path, compiled_regex);
        }
    }
}


bitvector drive_query(query_arguments &cmd_args, const bool &model)
{
    double t1, t2;
    // Load index from disk
    seqan3::debug_stream << "Reading Index from Disk... ";
    IndexStructure ibf;
    t1 = omp_get_wtime();
    load_ibf(ibf, cmd_args.idx);
    t2 = omp_get_wtime();
    seqan3::debug_stream << "DONE in " << t2-t1 << "s" << std::endl;
    double load_time = t2-t1;
    // std::cout << load_time << ",";

    cmd_args.query = translate(cmd_args.regex);
    auto bin_count = ibf.getBinCount();

    // Evaluate and search for Regular Expression
    seqan3::debug_stream << "Querying:" << std::endl;
    t1 = omp_get_wtime();
    uint8_t qlength = ibf.k_;
    std::string rx = cmd_args.regex;
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

    // Spawn IBF membership agent in this scope because it is expensive
    auto && ibf_ref = ibf.getIBF();
    auto agent = ibf_ref.membership_agent();

    std::vector<std::vector<std::string>> matrix{};
    for(auto i : knfa)
    {
        if(ibf.molecule_ == "na") {
            dfs<seqan3::dna5>(i, matrix, hash_to_bitvector, agent);
        } else {
            dfs<seqan3::aa27>(i, matrix, hash_to_bitvector, agent);
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

    for (auto path : paths_vector)
    {
        auto hits = query_ibf(bin_count, hash_to_bitvector, path);
        hit_vector.raw_data() |= hits.raw_data();
    }
    seqan3::debug_stream << "DONE" << std::endl;

    if(model)
    {
        /////////// MODELING STEP ///////////////
        double hit_count = 0;
        for(auto &&bit: hit_vector)
            hit_count+= bit;
        size_t query_length = matrix[0].size()+qlength-1;
        int text_length = cmd_args.text_length;
        size_t multiplyer = matrix.size();
        double result = compute_knut_model(query_length, qlength, text_length, multiplyer);
        seqan3::debug_stream << "FINAL PROBABILITY: " << result << std::endl;
        seqan3::debug_stream << "ACTUAL RATE: " << hit_count/hit_vector.size() << std::endl;
        /////////// MODELING STEP ///////////////
        return hit_vector; // Modeling doesn't require verification step
    }

    seqan3::debug_stream << "Verifying hits... " << std::endl;

    if(ibf.acid_libs_.size() > 1)
    {
        iter_disk_search(hit_vector, rx, ibf);
    } else
    {
        single_disk_search(hit_vector, rx, ibf);
    }
    t2 = omp_get_wtime();
    double search_time = t2-t1;
    double run_time = search_time + load_time;
    // std::cout << run_time << std::endl;
    seqan3::debug_stream << "DONE in " << run_time << "s" << std::endl;
    return hit_vector;
}
