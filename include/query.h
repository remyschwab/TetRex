//
// Created by Remy Schwab on 20.09.22.
//

#pragma once
#include <syncstream>
#include <re2/re2.h>
#include <re2/set.h>
#include <omp.h>

#include "kseq.h"
#include "utils.h"
#include "index_base.h"
#include "arg_parse.h"
#include "otf_collector.h"
#include "construct_nfa.h"
#include "construct_reduced_nfa.h"
#include "construction_tools.h"

double compute_k_probability(const uint8_t &k);

std::vector<std::string> split_string(const std::string& str, char delimiter);

double compute_knut_model(const size_t &query_length, const uint8_t &k, const int &m, const size_t &multiplyer);

void query_ibf_dna(query_arguments &cmd_args, const bool &model);

void query_ibf_aa(query_arguments &cmd_args, const bool &model);

void query_hibf_dna(query_arguments &cmd_args, const bool &model);

void query_hibf_aa(query_arguments &cmd_args, const bool &model);

void drive_query(query_arguments &cmd_args, const bool &model);

void reduce_query_alphabet(std::string &regex, const std::array<char, 256> &reduction_map);

std::vector<std::string> read_regex_from_file(const std::string &file_path);

std::string compute_reverse_complement(std::string &regex);

std::string complementBasesInRegex(const std::string &regex);

bool validate_regex(const std::string &regex, uint8_t ksize);

void reverse_verify_fasta_hit(const gzFile &fasta_handle, const re2::RE2 &crx, std::string const &binid);

void verify_fasta_hit(const gzFile &fasta_handle, const re2::RE2 &crx, std::string const &binid);

void verify_aa_fasta_hit(const gzFile &fasta_handle, const re2::RE2 &crx, std::string const &binid, const uint8_t &reduction, const std::array<char, 256> &residue_map);

void verify_fasta_set(const gzFile &fasta_handle, const RE2::Set &reg_set, std::string const &binid, const std::vector<std::string> &queries);

std::vector<size_t> compute_set_bins(const bitvector &hits, const std::vector<std::string> &acid_lib);

using alphamap = std::array<char, 256>;

template<index_structure::is_valid flavor, molecules::is_molecule mol_type>
void preprocess_query(std::string &rx_query, std::string &postfix_query, const TetrexIndex<flavor, mol_type> &ibf)
{
    if(ibf.reduction_ > 0) reduce_query_alphabet(rx_query, ibf.decomposer_.decomposer_.redmap_);
    // DBG(rx_query);
    postfix_query = translate(rx_query);
}


template<index_structure::is_valid flavor, molecules::is_dna mol_type>
void preprocess_query(std::string &rx_query, std::string &postfix_query, const TetrexIndex<flavor, mol_type> &ibf)
{
    postfix_query = translate(rx_query);
}


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
    re2::RE2 compiled_regex(query, re2::RE2::POSIX);
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
            gzclose(lib_path);
            continue;
        }
        verify_fasta_hit(lib_path, compiled_regex, ibf.acid_libs_[hit]);
        gzclose(lib_path);
    }
}


template<index_structure::is_valid flavor, molecules::is_molecule mol_type>
void iter_disk_search_set(const bitvector &hits, const std::vector<std::string> &queries, const TetrexIndex<flavor, mol_type> &ibf)
{
    std::vector<size_t> bins = compute_set_bins(hits, ibf.acid_libs_);
    RE2::Set regex_set(RE2::DefaultOptions, RE2::UNANCHORED);
    std::string err;
    for(auto query: queries)
    {
        query = "(" + query + ")"; // Capture entire RegEx
        int index = regex_set.Add(query, &err);
        if (index < 0)
        {
            seqan3::debug_stream << "Failed to add RegEx" << std::endl;
            return;
        }
    }
    if (!regex_set.Compile())
    {
        seqan3::debug_stream << "Failed to Compile RegEx" << std::endl;
        return;
    }
    #pragma omp parallel for
    for(size_t hit: bins)
    {
        gzFile lib_path = gzopen(ibf.acid_libs_[hit].c_str(), "r");
        if(!lib_path)
        {
            seqan3::debug_stream << ibf.acid_libs_[hit].c_str() << std::endl;
            throw std::runtime_error("File not found. Did you move/rename an indexed file?");
        }
        verify_fasta_set(lib_path, regex_set, ibf.acid_libs_[hit], queries);
        gzclose(lib_path);
    }
}

template<index_structure::is_valid flavor, molecules::is_molecule mol_t>
bitvector process_query(std::string &regex, TetrexIndex<flavor, mol_t> &ibf, const bool draw, const bool augment, const bool verbose)
{
    std::string query;
    preprocess_query(regex, query, ibf);
    // DBG(query);
    bitvector hit_vector(ibf.getBinCount(), true);
    std::unique_ptr<nfa_t> NFA = std::make_unique<nfa_t>();
    std::unique_ptr<lmap_t> nfa_map =  std::make_unique<lmap_t>(*NFA);
    amap_t arc_map;
    catsites_t catsites;
    if(ibf.reduction_ == Base) catsites = construct_kgraph(query, *NFA, *nfa_map, arc_map, ibf.k_, verbose);
    else catsites = construct_reduced_kgraph(query, *NFA, *nfa_map, arc_map, ibf.k_);
    std::unique_ptr<gmap_t> gap_map = std::make_unique<gmap_t>(*NFA);
    OTFCollector<flavor, mol_t> collector(std::move(NFA), std::move(nfa_map), ibf, std::move(arc_map), std::move(gap_map));
    collector.analyze_complexity();
    if(augment && catsites.size() > 0) collector.augment(catsites);
    if(draw) collector.draw_graph(catsites, augment);
    hit_vector = collector.collect();
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
        hit_vector &= process_query(rx, ibf, cmd_args.draw, cmd_args.augment, cmd_args.verbose);
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
            // iter_disk_search(hit_vector, rx, ibf);
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
void run_conjunction(query_arguments &cmd_args, const std::vector<std::string> &queries, const bool &model, TetrexIndex<flavor, mol_t> &ibf)
{
    double t1, t2;
    std::string rx;
    t1 = omp_get_wtime();
    bitvector hit_vector(ibf.getBinCount(), true);
    if(ibf.getBinCount() > 1u) // If someone forgot to split up their DB into bins then there's no point in the TetRex algorithm
    {
        for(auto rx: queries) hit_vector &= process_query(rx, ibf,cmd_args.draw, cmd_args.augment, cmd_args.verbose);
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
            iter_disk_search_set(hit_vector, queries, ibf);
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
void run_multiple_queries(query_arguments &cmd_args, const std::vector<std::string> &queries, const bool &model, TetrexIndex<flavor, mol_t> &ibf)
{
    if(cmd_args.conjunction)
    {
        run_conjunction(cmd_args, queries, model, ibf);
        return;
    }
    for(auto query: queries)
    {
        cmd_args.regex = query;
        seqan3::debug_stream << "\n" << query << std::endl;
        seqan3::debug_stream << query;
        run_collection(cmd_args, model, ibf);
    }
}
