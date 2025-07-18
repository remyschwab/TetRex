//
// Created by Remy Schwab on 20.09.22.
//

#include "query.h"

char comp_tab[] = {
	  0,   1,	2,	 3,	  4,   5,	6,	 7,	  8,   9,  10,	11,	 12,  13,  14,	15,
	 16,  17,  18,	19,	 20,  21,  22,	23,	 24,  25,  26,	27,	 28,  29,  30,	31,
	 32,  33,  34,	35,	 36,  37,  38,	39,	 40,  41,  42,	43,	 44,  45,  46,	47,
	 48,  49,  50,	51,	 52,  53,  54,	55,	 56,  57,  58,	59,	 60,  61,  62,	63,
	 64, 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O',
	'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z',	91,	 92,  93,  94,	95,
	 64, 't', 'v', 'g', 'h', 'e', 'f', 'c', 'd', 'i', 'j', 'm', 'l', 'k', 'n', 'o',
	'p', 'q', 'y', 's', 'a', 'a', 'b', 'w', 'x', 'r', 'z', 123, 124, 125, 126, 127
};


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


std::vector<size_t> compute_set_bins(const bitvector &hits, const std::vector<std::string> &acid_lib)
{
    std::vector<size_t> hit_vector;
    // The IBF will automatically use a multiple of 64 for the number of bins
    // This will cause a segfault if someone just uses one bin
    if(acid_lib.size() == 1u)
    {
        hit_vector.push_back(0u);
        return hit_vector;
    }
    size_t const words = seqan::hibf::divide_and_ceil(hits.size(), 64u);
        uint64_t const * const bit_vector_ptr = hits.data();

        // Jump to the next 1 and return the number of jumped bits in value.
        auto jump_to_next_1bit = [](uint64_t & value)
        {
            auto const zeros = std::countr_zero(value);
            value >>= zeros; // skip number of zeros
            return zeros;
        };

        // Each iteration can handle 64 bits, i.e., one word.
        for (size_t iteration = 0; iteration < words; ++iteration)
        {
            uint64_t current_word = bit_vector_ptr[iteration];

            // For each set bit in the current word, add/subtract 1 to the corresponding bin.
            for (size_t bin = iteration * 64u; current_word != 0u; ++bin, current_word >>= 1)
            {
                // Jump to the next 1
                bin += jump_to_next_1bit(current_word);
                hit_vector.push_back(bin);
            }
        }
    return hit_vector;
}


void reduce_query_alphabet(std::string &regex, const std::array<char, 256> &reduction_map)
{
    for(size_t i = 0; i < regex.length(); ++i)
    {
        char residue = regex.at(i);
        if(std::isalpha(residue))
        {
            regex.at(i) = reduction_map[residue];
        }
    }
}


// Check if the RegEx is smaller than the kmer size
bool validate_regex(const std::string &regex, uint8_t ksize)
{
    size_t dot_count = 1; // One dot means there should be two characters
    for(unsigned char c: regex) if(c == '-') ++dot_count;
    return dot_count >= ksize ? true : false;
}


void reverse_verify_fasta_hit(const gzFile &fasta_handle, const re2::RE2 &crx, std::string const &binid)
{
    int status;
    std::string match;
    kseq_t* record = kseq_init(fasta_handle);
    while((status = kseq_read(record)) >= 0)
    {
        int c0, c1;
        for(size_t i = 0; i < record->seq.l>>1; ++i)
        { // reverse complement sequence
            c0 = comp_tab[(int)record->seq.s[i]];
            c1 = comp_tab[(int)record->seq.s[record->seq.l - 1 - i]];
            record->seq.s[i] = c1;
            record->seq.s[record->seq.l - 1 - i] = c0;
        }
        if (record->seq.l & 1) // complement the remaining base
	        record->seq.s[record->seq.l>>1] = comp_tab[(int)record->seq.s[record->seq.l>>1]];
        re2::StringPiece bin_content(record->seq.s);
        while (RE2::FindAndConsume(&bin_content, crx, &match))
        {
            std::osyncstream(std::cout) << binid << "\t>" << record->name.s << "\t" << match << "\t" << "REVERSE STRAND HIT" << std::endl;
        }
    }
    kseq_destroy(record);
}


void verify_fasta_hit(const gzFile &fasta_handle, const re2::RE2 &crx, std::string const &binid)
{
    int status;
    std::string match;
    kseq_t* record = kseq_init(fasta_handle);
    while((status = kseq_read(record)) >= 0)
    {
        re2::StringPiece bin_content(record->seq.s);
        while (RE2::FindAndConsume(&bin_content, crx, &match))
        {
            std::osyncstream(std::cout) << binid << "\t>" << record->name.s << "\t" << match << std::endl;
        }
    }
    kseq_destroy(record);
}


void verify_reduced_fasta_hit(const gzFile &fasta_handle, const re2::RE2 &crx, std::string const &binid, const uint8_t &reduction, const std::array<char, 256> &residue_map)
{
    int status;
    size_t startpos;
    size_t endpos;
    re2::StringPiece match;
    std::string seq_copy;
    kseq_t* record = kseq_init(fasta_handle);
    while((status = kseq_read(record)) >= 0)
    {
        startpos = 0;
        endpos = record->seq.l;
        seq_copy = record->seq.s;
        for(size_t i = 0; i < record->seq.l; ++i)
        {
            char residue = record->seq.s[i];
            record->seq.s[i] = residue_map[residue];
        }
        re2::StringPiece bin_content(record->seq.s);
        while(startpos < endpos)
        {
            if(crx.Match(bin_content, startpos, bin_content.size(), RE2::UNANCHORED, &match, 1))
            {
                std::osyncstream(std::cout) << binid << "\t>" << record->name.s << "\t" << seq_copy.substr(startpos, match.size()) << std::endl;
                startpos += match.size();
                continue;
            }
            startpos++;
        }
    }
    kseq_destroy(record);
}


void verify_fasta_set(const gzFile &fasta_handle, const RE2::Set &reg_set, std::string const &binid, const std::vector<std::string> &queries)
{
    int status;
    std::string match;
    kseq_t* record = kseq_init(fasta_handle);
    uint8_t count = queries.size();
    while((status = kseq_read(record)) >= 0)
    {
        std::vector<int> matching_rules;
        if (!reg_set.Match(record->seq.s, &matching_rules)) continue;
        if(matching_rules.size() == count)
        {
            std::osyncstream(std::cout) << binid << "\t>" << record->name.s << "\tN --> ";
            for(auto rule_index: matching_rules)
            {
                std::osyncstream(std::cout) << queries[rule_index] << " --> ";
            }
            std::osyncstream(std::cout) << "C" << std::endl;
        }
    }
    kseq_destroy(record);
}


std::vector<std::string> read_regex_from_file(const std::string &file_path)
{
    std::vector<std::string> queries;
    std::ifstream handle{file_path};
    std::string line;
    if (handle.is_open())
    {
        while (getline(handle, line))
        {
            queries.push_back(line);
        }
        handle.close();
    }
    else
    {
        std::cerr << "Filepath not valid" << std::endl;
    }
    return queries;
}

std::vector<std::string> split_string(const std::string& str, char delimiter)
{
    std::vector<std::string> result;
    std::stringstream ss(str);
    std::string token;

    while (std::getline(ss, token, delimiter)) {
        result.push_back(token);
    }
    
    return result;
}


void query_ibf_dna(query_arguments &cmd_args, const bool &model)
{
    omp_set_num_threads(cmd_args.t);
    TetrexIndex<index_structure::IBF, molecules::nucleotide> ibf;
    load_ibf(ibf, cmd_args.idx);
    ibf.spawn_agent();
    if(cmd_args.read_file)
    {
        std::vector<std::string> queries = read_regex_from_file(cmd_args.regex);
        run_multiple_queries(cmd_args, queries, model, ibf);
        return;
    }
    else if(cmd_args.conjunction)
    {
        std::vector<std::string> queries = split_string(cmd_args.regex, ':');
        if(queries.size() == 1)
        {
            seqan3::debug_stream << "Did you use the correct delimiter (:)?" << std::endl;
            return;
        }
        run_multiple_queries(cmd_args, queries, model, ibf);
        return;
    }
    run_collection(cmd_args, model, ibf);
}

void query_ibf_aa(query_arguments &cmd_args, const bool &model)
{
    omp_set_num_threads(cmd_args.t);
    TetrexIndex<index_structure::IBF, molecules::peptide> ibf;
    load_ibf(ibf, cmd_args.idx);
    ibf.spawn_agent();
    if(cmd_args.read_file)
    {
        std::vector<std::string> queries = read_regex_from_file(cmd_args.regex);
        run_multiple_queries(cmd_args, queries, model, ibf);
        return;
    }
    else if(cmd_args.conjunction)
    {
        std::vector<std::string> queries = split_string(cmd_args.regex, ':');
        if(queries.size() == 1)
        {
            seqan3::debug_stream << "Did you use the correct delimiter (:)?" << std::endl;
            return;
        }
        run_multiple_queries(cmd_args, queries, model, ibf);
        return;
    }
    run_collection(cmd_args, model, ibf);
}

void query_hibf_dna(query_arguments &cmd_args, const bool &model)
{
    omp_set_num_threads(cmd_args.t);
    TetrexIndex<index_structure::HIBF, molecules::nucleotide> ibf;
    load_ibf(ibf, cmd_args.idx);
    if(cmd_args.read_file)
    {
        std::vector<std::string> queries = read_regex_from_file(cmd_args.regex);
        run_multiple_queries(cmd_args, queries, model, ibf);
        return;
    }
    else if(cmd_args.conjunction)
    {
        std::vector<std::string> queries = split_string(cmd_args.regex, ':');
        if(queries.size() == 1)
        {
            seqan3::debug_stream << "Did you use the correct delimiter (:)?" << std::endl;
            return;
        }
        run_multiple_queries(cmd_args, queries, model, ibf);
        return;
    }
    run_collection(cmd_args, model, ibf);
}

void query_hibf_aa(query_arguments &cmd_args, const bool &model)
{
    omp_set_num_threads(cmd_args.t);
    TetrexIndex<index_structure::HIBF, molecules::peptide> ibf;
    load_ibf(ibf, cmd_args.idx);
    if(cmd_args.read_file)
    {
        std::vector<std::string> queries = read_regex_from_file(cmd_args.regex);
        run_multiple_queries(cmd_args, queries, model, ibf);
        return;
    }
    else if(cmd_args.conjunction)
    {
        std::vector<std::string> queries = split_string(cmd_args.regex, ':');
        if(queries.size() == 1)
        {
            seqan3::debug_stream << "Did you use the correct delimiter (:)?" << std::endl;
            return;
        }
        run_multiple_queries(cmd_args, queries, model, ibf);
        return;
    }
    run_collection(cmd_args, model, ibf);
}

void drive_query(query_arguments &cmd_args, const bool &model)
{
    if(cmd_args.regex == "-") std::cin >> cmd_args.regex;
    index_params params;
    load_params(params, cmd_args.idx);
    if(!params.is_hibf_ && params.molecule_ == "na")
    {
        query_ibf_dna(cmd_args, model);
    }
    else if(params.is_hibf_ && params.molecule_ == "na")
    {
        query_hibf_dna(cmd_args, model);
    }
    else if(!params.is_hibf_ && params.molecule_ == "aa")
    {
        query_ibf_aa(cmd_args, model);
    }
    else if(params.is_hibf_ && params.molecule_ == "aa")
    {
        query_hibf_aa(cmd_args, model);
    }
}
