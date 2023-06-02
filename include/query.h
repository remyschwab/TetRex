//
// Created by Remy Schwab on 20.09.22.
//

#pragma once

#include <seqan3/alphabet/views/all.hpp>

#include <re2/re2.h>
#include <omp.h>

#include "kseq.h"
#include "utils.h"
#include "index.h"
#include "arg_parse.h"
#include "korotkov_nfa.h"
#include "graphMaker.h"
#include "nfa_pointer.h"


bitvector query_ibf(uint32_t &bin_count, robin_hood::unordered_map<uint64_t, bitvector> &hash_to_bits, std::vector<std::pair<std::string, uint64_t>> &path);

double compute_k_probability(const uint8_t &k);

double compute_knut_model(const size_t &query_length, const uint8_t &k, const int &m, const size_t &multiplyer);

void drive_query(query_arguments &cmd_args, const bool &model);

void preprocess_query(std::string &rx_query, std::string &postfix_query);

void verify_fasta_hit(const std::filesystem::path &bin_path, re2::RE2 &crx);

void iter_disk_search(const bitvector &hits, const std::string &query, IndexStructure &ibf);
