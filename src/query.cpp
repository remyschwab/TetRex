//
// Created by Remy Schwab on 20.09.22.
//

#include "query.h"
KSEQ_INIT(gzFile, gzread)


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


bitvector query_ibf(size_t &bin_count, robin_hood::unordered_map<uint64_t, bitvector> &hash_to_bits, std::vector<uint64_t> &path)
{
    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed>::membership_agent_type::binning_bitvector hit_vector{bin_count};
    std::fill(hit_vector.begin(), hit_vector.end(), true);
    for (auto && kmer : path)
    {
        auto & result = hash_to_bits[kmer];
        hit_vector.raw_data() &= result.raw_data();
    }
    return hit_vector;
}


void preprocess_query(std::string &rx_query, std::string &postfix_query)
{
    // We don't want to generate kmers from something with anchors

    // We want the entire query to be in one capture group
    // But we need to account for the case where the user tried that themselves
    // size_t query_length = rx_query.length();
    // if(rx_query[0] != "(" && rx_query[query_length-1] != ")") // Default case where there is just a query
    // {
    //     postfix_query = translate(rx_query);
    //     rx_query = "(" + rx_query + ")";
    // }
    // seqan3::debug_stream << rx_query << std::endl;
    postfix_query = translate(rx_query);
    rx_query = "(" + rx_query + ")"; // Capture entire RegEx
}


void verify_fasta_hit(const gzFile &fasta_handle, kseq_t *record, re2::RE2 &crx)
{
    int status;
    int start = 0;
    std::string match;
    record = kseq_init(fasta_handle);
    while((status = kseq_read(record)) >= 0)
    {
        re2::StringPiece bin_content(record->seq.s);
        while (RE2::FindAndConsume(&bin_content, crx, &match))
        {
            std::cout << ">" << record->name.s << "\t" << match << "\t" << start << "-" << start+match.length() << std::endl;
            start++;
        }
    }
}


void iter_disk_search(const bitvector &hits, const std::string &query, IndexStructure &ibf)
{
    size_t bins = hits.size();
    gzFile lib_path;
    kseq_t *record;

    re2::RE2 compiled_regex(query);
    assert(compiled_regex.ok());
    #pragma omp parallel for
    for(size_t i = 0; i < bins; i++)
    {
        if(hits[i])
        {
            lib_path = gzopen(ibf.acid_libs_[i].c_str(), "r");
            verify_fasta_hit(lib_path, record, compiled_regex);
        }
    }
    kseq_destroy(record);
    gzclose(lib_path);
}


bitvector drive_query(query_arguments &cmd_args, const bool &model)
{
    double t1, t2, t3;
    omp_set_num_threads(cmd_args.t);
    // Load index from disk
    seqan3::debug_stream << "Reading Index from Disk... ";
    IndexStructure ibf;
    t1 = omp_get_wtime();
    load_ibf(ibf, cmd_args.idx);

    if(ibf.molecule_ == "na") //TODO: Update the serializing logic
    {
        ibf.create_selection_bitmask();
        ibf.set_left_shift();
    }
    else
    {
        ibf.set_alphabet_maps();
        ibf.compute_powers();
    }

    t2 = omp_get_wtime();
    double load_time = t2-t1;
    seqan3::debug_stream << "DONE in " << load_time << "s" << std::endl;
    // std::cout << load_time << ",";

    // The RegEx is given in InFix notation and is translated to postfix notation for internal use

    // Evaluate and search for Regular Expression
    seqan3::debug_stream << "Querying:" << std::endl;
    auto && bin_count = ibf.getBinCount();
    uint8_t &qlength = ibf.k_;
    std::string &rx = cmd_args.regex;
    std::string &query = cmd_args.query;
    preprocess_query(rx, query);

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
    robin_hood::unordered_map<uint64_t, bitvector> hash_to_bitvector{};

    // Spawn IBF membership agent in this scope because it is expensive
    auto && ibf_ref = ibf.getIBF();
    auto agent = ibf_ref.membership_agent();

    std::vector<std::vector<uint64_t>> matrix{};
    for(auto i : knfa)
    {
        if(ibf.molecule_ == "na")
        {
            dfs(i, matrix, hash_to_bitvector, agent, ibf);
        }
        else
        {
            dfs(i, matrix, hash_to_bitvector, agent, ibf);
        }
    }

    uMatrix(matrix);
    seqan3::debug_stream << "DONE" << std::endl;


    // Search kmer paths in index
    seqan3::debug_stream << "   - Search kmers in index... ";

    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed>::membership_agent_type::binning_bitvector hit_vector{ibf.getBinCount()};
    std::fill(hit_vector.begin(), hit_vector.end(), false);

    // #pragma omp parallel for
    for(size_t i = 0; i < matrix.size(); ++i)
    {
        auto hits = query_ibf(bin_count, hash_to_bitvector, matrix[i]);
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

    iter_disk_search(hit_vector, rx, ibf);
    t3 = omp_get_wtime();
    double run_time = t3-t2;
    // std::cout << run_time << std::endl;
    seqan3::debug_stream << "DONE in " << run_time << "s" << std::endl;
    return hit_vector;
}
