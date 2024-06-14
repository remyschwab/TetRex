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


// Function to reverse the regular expression while preserving groups
std::string reverseRegexPreservingGroups(const std::string &regex)
{
    std::stack<std::string> stack;
    std::string reversedRegex;
    int i = 0;
    int n = regex.length();
    
    while (i < n) {
        if (regex[i] == '(')
        {
            std::string group = "(";
            i++;
            while (i < n && regex[i] != ')') {
                group += regex[i];
                i++;
            }
            group += ')';
            stack.push(group);
        }
        else
        {
            stack.push(std::string(1, regex[i]));
        }
        i++;
    }
    
    while(!stack.empty())
    {
        reversedRegex += stack.top();
        stack.pop();
    }
    return reversedRegex;
}


// Function to complement the bases in the regular expression
std::string complementBasesInRegex(const std::string &regex)
{
    std::unordered_map<char, char> complement = {
        {'A', 'T'}, {'T', 'A'}, {'C', 'G'}, {'G', 'C'}
    };
    
    std::string complementedRegex;
    int i = 0;
    int n = regex.length();
    
    while (i < n) {
        if (regex[i] == '(') {
            std::string group = "(";
            i++;
            while (i < n && regex[i] != ')') {
                group += complement.count(regex[i]) ? complement[regex[i]] : regex[i];
                i++;
            }
            group += ')';
            complementedRegex += group;
        } else {
            complementedRegex += complement.count(regex[i]) ? complement[regex[i]] : regex[i];
        }
        i++;
    }
    
    return complementedRegex;
}


std::string compute_reverse_complement(std::string &regex)
{
    std::string forward_reverse_regex = "("+regex+")|(";
    std::string reverse_regex = reverseRegexPreservingGroups(regex);
    forward_reverse_regex += complementBasesInRegex(reverse_regex)+")";
    return forward_reverse_regex;
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
}

// Check if the RegEx is smaller than the kmer size
bool validate_regex(const std::string &regex, uint8_t ksize)
{
    size_t dot_count = 1; // One dot means there should be two characters
    for(unsigned char c: regex) if(c == '-') ++dot_count;
    return dot_count >= ksize ? true : false;
}


void verify_fasta_hit(const gzFile &fasta_handle, kseq_t *record, re2::RE2 &crx, std::string const &binid)
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
            std::cout << binid << "\t>" << record->name.s << "\t" << match << std::endl;
            ++start;
        }
    }
}


void query_ibf_dna(query_arguments &cmd_args, const bool &model)
{
    omp_set_num_threads(cmd_args.t);
    TetrexIndex<index_structure::IBF, molecules::nucleotide> ibf;
    load_ibf(ibf, cmd_args.idx);
    run_collection(cmd_args, model, ibf);
}

void query_ibf_aa(query_arguments &cmd_args, const bool &model)
{
    omp_set_num_threads(cmd_args.t);
    TetrexIndex<index_structure::IBF, molecules::peptide> ibf;
    load_ibf(ibf, cmd_args.idx);
    run_collection(cmd_args, model, ibf);
}

void query_hibf_dna(query_arguments &cmd_args, const bool &model)
{
    omp_set_num_threads(cmd_args.t);
    TetrexIndex<index_structure::HIBF, molecules::nucleotide> ibf;
    load_ibf(ibf, cmd_args.idx);
    run_collection(cmd_args, model, ibf);
}

void query_hibf_aa(query_arguments &cmd_args, const bool &model)
{
    omp_set_num_threads(cmd_args.t);
    TetrexIndex<index_structure::HIBF, molecules::peptide> ibf;
    load_ibf(ibf, cmd_args.idx);
    run_collection(cmd_args, model, ibf);
}

void drive_query(query_arguments &cmd_args, const bool &model)
{
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
