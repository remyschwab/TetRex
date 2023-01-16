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


// Types
template <typename MoleculeType> using record_pair = std::pair<std::string, MoleculeType>;
template <typename MoleculeType> using record_list = std::vector<record_pair<MoleculeType>>;

using bitvector = seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed>::membership_agent::binning_bitvector;

char* re2post(char *re);

std::string stream_as_string(const std::string& path);

int matches(const std::string& bin, std::regex reg, std::fstream& writefile);

std::string translate(const std::string& str);

std::vector<char> getAlphabet(const std::string& regex);

std::vector<seqan3::dna5> convertStringToDNA(std::string const &str);

//soll alle nötigen qgramme finden
std::vector<std::string> getQgramAlphabet(const std::vector<std::vector<std::string>>& matrix);

//uint shiftValue(const uint& input);


//uint32_t getHash(const std::vector<char>& alphabet, const std::string& input, const uint& sValue);


void matrixTotxt(const std::vector<std::vector<std::string>>& matrix, std::string& filename);


void matrixTXT(const std::vector<std::vector<std::string>>& matrix, const std::vector<char>& alphabet);
