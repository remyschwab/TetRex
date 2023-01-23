//
// Created by Remy Schwab on 20.09.22.
//

#ifndef KBIOREG_QUERY_H
#define KBIOREG_QUERY_H

#include "utils.h"
#include "index.h"
#include "arg_parse.h"
#include "korotkov_nfa.h"
#include "graphMaker.h"
#include "nfa_pointer.h"


bitvector query_ibf(uint32_t &bin_count, robin_hood::unordered_map<uint64_t, uint32_t> &hash_to_idx,
  std::vector<bitvector> &kmer_bitvex, std::vector<std::pair<std::string, uint64_t>> &path);

bitvector drive_query(const query_arguments & cmd_args);

bitvector drive_query_benchmark(const query_arguments & cmd_args, std::fstream &benchmark_table);

#endif //KBIOREG_QUERY_H
