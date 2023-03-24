//
// Created by Remy Schwab on 20.09.22.
//

#pragma once

#include <seqan3/alphabet/views/all.hpp>

#include <re2/re2.h>
#include <omp.h>

#include "utils.h"
#include "index.h"
#include "arg_parse.h"
#include "korotkov_nfa.h"
#include "graphMaker.h"
#include "nfa_pointer.h"


bitvector query_ibf(uint32_t &bin_count, robin_hood::unordered_map<uint64_t, bitvector> &hash_to_bits, std::vector<std::pair<std::string, uint64_t>> &path);

double compute_k_probability(const uint8_t &k);

double compute_knut_model(const size_t &query_length, const uint8_t &k, const int &m, const size_t &multiplyer);

bitvector drive_query(query_arguments &cmd_args, const bool &model);

void verify_fasta_hit(const std::filesystem::path &bin_path, re2::RE2 &crx);

void iter_disk_search(const bitvector &hits, const std::string &query, IndexStructure &ibf);

template <seqan3::alphabet Alphabet>
auto to_string(std::vector<Alphabet> const& sequence) -> std::string {
    auto view = sequence | std::views::transform([](Alphabet a) {
        return a.to_char();
    });
    return std::string{view.begin(), view.end()};
}

template <typename MolType>
void extract_matrix_paths(const std::vector<std::vector<std::string>> &matrix,
 path_vector &paths_vector, auto &hash_adaptor)
{
    for(auto const& i : matrix)
    {
        std::vector<std::pair<std::string, uint64_t>> hash_vector;
        for(auto const& j : i)
        {
            std::vector<MolType> acid_vec = convertStringToAlphabet<MolType>(j);
            auto digest = acid_vec | hash_adaptor;
            // Create a vector of kmer hashes that correspond
            hash_vector.emplace_back(j, digest[0]);
        }
        paths_vector.emplace_back(std::move(hash_vector));
    }
}
