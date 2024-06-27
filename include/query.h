//
// Created by Remy Schwab on 20.09.22.
//

#pragma once

#include <re2/re2.h>
#include <omp.h>

#include "kseq.h"
#include "utils.h"
#include "index_base.h"
#include "arg_parse.h"
#include "otf_collector.h"
#include "construct_nfa.h"


// bitvector query_ibf(uint32_t &bin_count, robin_hood::unordered_map<uint64_t, bitvector> &hash_to_bits, std::vector<std::pair<std::string, uint64_t>> &path);

double compute_k_probability(const uint8_t &k);

double compute_knut_model(const size_t &query_length, const uint8_t &k, const int &m, const size_t &multiplyer);

void query_ibf_dna(query_arguments &cmd_args, const bool &model);

void query_ibf_aa(query_arguments &cmd_args, const bool &model);

void query_hibf_dna(query_arguments &cmd_args, const bool &model);

void query_hibf_aa(query_arguments &cmd_args, const bool &model);

void drive_query(query_arguments &cmd_args, const bool &model);

void preprocess_query(std::string &rx_query, std::string &postfix_query);

std::string compute_reverse_complement(std::string &regex);

std::string complementBasesInRegex(const std::string &regex);

bool validate_regex(const std::string &regex, uint8_t ksize);

void verify_fasta_hit(const gzFile &fasta_handle, kseq_t *record, re2::RE2 &crx, std::string const &binid);

template<index_structure::is_valid flavor, molecules::is_dna mol_type>
void iter_disk_search(const bitvector &hits, std::string &query, TetrexIndex<flavor, mol_type> &ibf)
{
    size_t bins = hits.size();
    gzFile lib_path;
    kseq_t *record;

    std::string forward_and_reverse = compute_reverse_complement(query);
    forward_and_reverse = "(" + forward_and_reverse + ")"; // Capture entire RegEx

    re2::RE2 compiled_regex(forward_and_reverse);
    assert(compiled_regex.ok());
    #pragma omp parallel for
    for(size_t i = 0; i < bins; i++)
    {
        if(hits[i])
        {
            lib_path = gzopen(ibf.acid_libs_[i].c_str(), "r");
            verify_fasta_hit(lib_path, record, compiled_regex, ibf.acid_libs_[i]);
        }
    }
    kseq_destroy(record);
    gzclose(lib_path);
}

template<index_structure::is_valid flavor, molecules::is_peptide mol_type>
void iter_disk_search(const bitvector &hits, std::string &query, TetrexIndex<flavor, mol_type> &ibf)
{
    size_t bins = hits.size();
    gzFile lib_path;
    kseq_t *record;
    query = "(" + query + ")"; // Capture entire RegEx
    re2::RE2 compiled_regex(query);
    assert(compiled_regex.ok());
    #pragma omp parallel for
    for(size_t i = 0; i < bins; i++)
    {
        if(hits[i])
        {
            lib_path = gzopen(ibf.acid_libs_[i].c_str(), "r");
            verify_fasta_hit(lib_path, record, compiled_regex, ibf.acid_libs_[i]);
        }
    }
    kseq_destroy(record);
    gzclose(lib_path);
}

template<index_structure::is_valid flavor, molecules::is_molecule mol_t>
void run_collection(query_arguments &cmd_args, const bool &model, TetrexIndex<flavor, mol_t> &ibf)
{
    double t1, t2;
    t1 = omp_get_wtime();
    std::string &rx = cmd_args.regex;
    std::string &query = cmd_args.query;
    preprocess_query(rx, query);
    bool valid = validate_regex(query, ibf.k_);
    // seqan3::debug_stream << query << std::endl;
    bitvector hit_vector(ibf.getBinCount(), true);

    if(valid)
    {
        std::unique_ptr<nfa_t> NFA = std::make_unique<nfa_t>();
        std::unique_ptr<lmap_t> nfa_map =  std::make_unique<lmap_t>(*NFA);
        amap_t arc_map;
        
        construct_kgraph(cmd_args.query, *NFA, *nfa_map, arc_map, ibf.k_);

        // print_node_ids(*NFA, *nfa_map);
        // seqan3::debug_stream << std::endl;
        // print_node_pointers(arc_map, *NFA);
        // seqan3::debug_stream << std::endl;
        // print_kgraph_arcs(*NFA);
        // seqan3::debug_stream << std::endl;

        std::vector<int> top_rank_map = run_top_sort(*NFA);
        OTFCollector<flavor, mol_t> collector(std::move(NFA), std::move(nfa_map),
                                            ibf,
                                            std::move(top_rank_map), std::move(arc_map));
        
        hit_vector = collector.collect();
        if(cmd_args.verbose) seqan3::debug_stream << "Narrowed Search to " << collector.sumBitvector(hit_vector) << " possible bins" << std::endl;
    }
    else // if the RegEx is shorter than the index kmer size, then prompt user and trigger linear search
    {
        seqan3::debug_stream << "RegEx is too short to use index. Performing linear scan over whole database" << std::endl;
    }
    if(!hit_vector.none()) iter_disk_search(hit_vector, rx, ibf);
    t2 = omp_get_wtime();
    seqan3::debug_stream << "Query Time: " << (t2-t1) << std::endl;
}
