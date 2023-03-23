#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <stack>
#include <vector>
#include <filesystem>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <cstring>
#include <regex>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/views/kmer_hash.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>
#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>


/////////////// Type Declarations ///////////////
template <typename MoleculeType> using record_pair = std::pair<std::string, MoleculeType>;
template <typename MoleculeType> using record_list = std::vector<record_pair<MoleculeType>>;
using bitvector = seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed>::membership_agent::binning_bitvector;
using path_vector = std::vector<std::vector<std::pair<std::string, uint64_t>>>;
/////////////// ****** END ****** ///////////////

char* re2post(char *re);

std::string stream_as_string(const std::string& path);

int matches(const std::string& bin, std::regex reg);

std::string translate(const std::string& str);

std::vector<char> getAlphabet(const std::string& regex);

template <typename Alphabet>
auto convertStringToAlphabet(std::string const &str)
{
    std::vector<Alphabet> alphabet_str;
    alphabet_str.reserve(str.size());
    for (auto s : str) {
        alphabet_str.emplace_back(seqan3::assign_char_to(s, Alphabet{}));
    }
    return alphabet_str;
}

//soll alle nötigen qgramme finden
std::vector<std::string> getQgramAlphabet(const std::vector<std::vector<std::string>>& matrix);

void matrixTotxt(const std::vector<std::vector<std::string>>& matrix, std::string& filename);

void matrixTXT(const std::vector<std::vector<std::string>>& matrix, const std::vector<char>& alphabet);
