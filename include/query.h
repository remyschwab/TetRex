//
// Created by Remy Schwab on 20.09.22.
//

#pragma once
#include <syncstream>
#include <re2/re2.h>
#include <omp.h>

#include "kseq.h"
#include "utils.h"
#include "index_base.h"
#include "arg_parse.h"
#include "otf_collector.h"
#include "construction_tools.h"
#include "construct_nfa.h"
#include "construct_reduced_nfa.h"

double compute_k_probability(const uint8_t &k);

double compute_knut_model(const size_t &query_length, const uint8_t &k, const int &m, const size_t &multiplyer);

void query_ibf_dna(query_arguments &cmd_args, const bool &model);

void query_ibf_aa(query_arguments &cmd_args, const bool &model);

void query_hibf_dna(query_arguments &cmd_args, const bool &model);

void query_hibf_aa(query_arguments &cmd_args, const bool &model);

void drive_query(query_arguments &cmd_args, const bool &model);

void reduce_query_alphabet(std::string &regex, const std::array<char, 256> &reduction_map);

std::vector<std::string> read_regex_from_file(const std::string &file_path);

template<index_structure::is_valid flavor, molecules::is_peptide mol_type>
void preprocess_query(std::string &rx_query, std::string &postfix_query, const TetrexIndex<flavor, mol_type> &ibf)
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
    if(ibf.reduction_ > 0) reduce_query_alphabet(rx_query, ibf.decomposer_.decomposer_.redmap_);
    postfix_query = translate(rx_query);
}

template<index_structure::is_valid flavor, molecules::is_dna mol_type>
void preprocess_query(std::string &rx_query, std::string &postfix_query, const TetrexIndex<flavor, mol_type> &ibf)
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

std::string compute_reverse_complement(std::string &regex);

std::string complementBasesInRegex(const std::string &regex);

bool validate_regex(const std::string &regex, uint8_t ksize);

void reverse_verify_fasta_hit(const gzFile &fasta_handle, const re2::RE2 &crx, std::string const &binid);

void verify_fasta_hit(const gzFile &fasta_handle, const re2::RE2 &crx, std::string const &binid);

void verify_aa_fasta_hit(const gzFile &fasta_handle, const re2::RE2 &crx, std::string const &binid, const uint8_t &reduction, const std::array<char, 256> &residue_map);

std::vector<size_t> compute_set_bins(const bitvector &hits, const std::vector<std::string> &acid_lib);

template<index_structure::is_valid flavor, molecules::is_dna mol_type>
void iter_disk_search(const bitvector &hits, std::string &query, const TetrexIndex<flavor, mol_type> &ibf)
{
    std::vector<size_t> bins = compute_set_bins(hits, ibf.acid_libs_);
    std::string forward_and_reverse = query;
    forward_and_reverse = "(" + forward_and_reverse + ")"; // Capture entire RegEx

    re2::RE2 compiled_regex(forward_and_reverse);
    assert(compiled_regex.ok());
    #pragma omp parallel for
    for(size_t hit: bins)
    {
        gzFile lib_path = gzopen(ibf.acid_libs_[hit].c_str(), "r");
        if(!lib_path)
        {
            throw std::runtime_error("File not found. Did you move/rename an indexed file?");
        }
        verify_fasta_hit(lib_path, compiled_regex, ibf.acid_libs_[hit]);
	lib_path = gzopen(ibf.acid_libs_[hit].c_str(), "r");
        reverse_verify_fasta_hit(lib_path, compiled_regex, ibf.acid_libs_[hit]);
        gzclose(lib_path);
    }
}

template<index_structure::is_valid flavor, molecules::is_peptide mol_type>
void iter_disk_search(const bitvector &hits, std::string &query, const TetrexIndex<flavor, mol_type> &ibf)
{
    std::vector<size_t> bins = compute_set_bins(hits, ibf.acid_libs_);
    query = "(" + query + ")"; // Capture entire RegEx
    re2::RE2 compiled_regex(query);
    assert(compiled_regex.ok());
    #pragma omp parallel for
    for(size_t hit: bins)
    {
        gzFile lib_path = gzopen(ibf.acid_libs_[hit].c_str(), "r");
        if(!lib_path)
        {
            seqan3::debug_stream << ibf.acid_libs_[hit].c_str() << std::endl;
            throw std::runtime_error("File not found. Did you move/rename an indexed file?");
        }
        if(ibf.reduction_ > Base)
        {
            verify_reduced_fasta_hit(lib_path, compiled_regex, ibf.acid_libs_[hit], ibf.reduction_, ibf.decomposer_.decomposer_.redmap_);
            continue;
        }
        verify_fasta_hit(lib_path, compiled_regex, ibf.acid_libs_[hit]);
        gzclose(lib_path);
    }
}

template<index_structure::is_valid flavor, molecules::is_molecule mol_t>
bitvector process_query(const std::string &regex, TetrexIndex<flavor, mol_t> &ibf)
{
    std::string rx = regex;
    std::string query;
    preprocess_query(rx, query, ibf);
    bool valid = validate_regex(query, ibf.k_);
    bitvector hit_vector(ibf.getBinCount(), true);
    if(valid)
    {
        std::unique_ptr<nfa_t> NFA = std::make_unique<nfa_t>();
        std::unique_ptr<lmap_t> nfa_map =  std::make_unique<lmap_t>(*NFA);
        amap_t arc_map;
        
        if(ibf.reduction_ == Base)
        {
            construct_kgraph(query, *NFA, *nfa_map, arc_map, ibf.k_);
        }
        else
        {
            construct_reduced_kgraph(query, *NFA, *nfa_map, arc_map, ibf.k_);
        }

        // if(cmd_args.draw) print_graph(*NFA, *nfa_map);

        std::vector<int> top_rank_map = run_top_sort(*NFA);
        OTFCollector<flavor, mol_t> collector(std::move(NFA), std::move(nfa_map),
                                            ibf,
                                            std::move(top_rank_map), std::move(arc_map));   
        hit_vector = collector.collect();
    }
    else // if the RegEx is shorter than the index kmer size, then prompt user and trigger linear search
    {
        seqan3::debug_stream << "RegEx is too short to use index. Performing linear scan over whole database" << std::endl;
    }
    return hit_vector;
}


template<index_structure::is_valid flavor, molecules::is_molecule mol_t>
void run_collection(query_arguments &cmd_args, const bool &model, TetrexIndex<flavor, mol_t> &ibf)
{
    double t1, t2;
    std::string &rx = cmd_args.regex;
    t1 = omp_get_wtime();
    bitvector hit_vector(ibf.getBinCount(), true);
    if(ibf.getBinCount() > 1u) // If someone forgot to split up their DB into bins then there's no point in the TetRex algorithm
    {
        hit_vector &= process_query(rx, ibf);
    }
    else
    {
        seqan3::debug_stream << "[WARNING] Index contains only 1 bin. Unable to accelerate search using the TetRex algorithm. Performing Linear Scan" << std::endl;
    }

    if(cmd_args.verbose) seqan3::debug_stream << "Narrowed Search to " << OTFCollector<flavor, mol_t>::sumBitVector(hit_vector) << " possible bins" << std::endl;
    
    if(!hit_vector.none())
    {
        try
        {
            iter_disk_search(hit_vector, rx, ibf);
        }
        catch(const std::exception& e)
        {
            std::cerr << e.what() << '\n';
        }
    }
    t2 = omp_get_wtime();
    seqan3::debug_stream << "Query Time: " << (t2-t1) << std::endl;
}


template<index_structure::is_valid flavor, molecules::is_molecule mol_t>
void run_multiple_queries(query_arguments &cmd_args, const bool &model, TetrexIndex<flavor, mol_t> &ibf)
{
    std::vector<std::string> queries = read_regex_from_file(cmd_args.regex);
    for(auto query: queries)
    {
        cmd_args.regex = query;
        seqan3::debug_stream << "\n" << query << std::endl;
        run_collection(cmd_args, model, ibf);
    }
}
