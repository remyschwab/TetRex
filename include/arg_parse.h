//
// Created by Remy Schwab on 12.09.22.
//

#pragma once

#include <sharg/all.hpp>


struct index_arguments
{
    uint8_t k = 3;
    uint8_t t = 1;
    float fpr = 0.05;
    std::string ibf = "hibf";
    std::string ofile;
    std::vector<std::filesystem::path> acid_libs{};
    size_t bin_size = 10000;
    uint8_t hash_count = 3;
    std::string molecule;
    std::string reduction = "None";
};

inline void initialise_index_parser(sharg::parser &parser, index_arguments &args)
{
    parser.info.author = "Remy Schwab";
    parser.info.version = "1.0.0";
    parser.add_option(args.ibf, sharg::config{.short_id='i', .long_id="idx_struct", .description="IBF or HIBF", .validator=sharg::value_list_validator{"ibf", "hibf"}});
    parser.add_option(args.k, sharg::config{'k', "ksize", "size of kmers"});
    parser.add_option(args.fpr, sharg::config{'p', "fpr", "Bloom Filter False Positive Rate"});
    parser.add_option(args.bin_size, sharg::config{'s', "bin_size", "Size of bins"});
    parser.add_option(args.hash_count, sharg::config{'c', "hash_count", "Number of hash functions. NOTE: MORE THAN 4 IS SLOW"});
    parser.add_option(args.t, sharg::config{'t', "threads", "Number of threads"});
    parser.add_option(args.molecule, sharg::config{.short_id='m', .long_id="molecule", .description="Molecule type of library", .required=true, .validator=sharg::value_list_validator{"na", "aa"}});
    parser.add_option(args.reduction, sharg::config{.short_id='r', .long_id="reduce", .description="Use reduced AA alphabet (Murphy or Li)", .validator=sharg::value_list_validator{"murphy","li"}});
    parser.add_option(args.ofile, sharg::config{'o', "ofile", "Name of index on disk"});
    parser.add_positional_option(args.acid_libs, sharg::config{.description="Nucleic or Amino Acid library to indexed", .validator=sharg::input_file_validator{{"lst","fa", "fa.gz","fasta", "fasta.gz", "fna", "fna.gz"}}});
}

struct query_arguments
{
    uint8_t t = 1;
    bool verbose = 0;
    int text_length;
    std::filesystem::path idx{};
    std::string regex;
    std::string query;
};

inline void initialise_query_parser(sharg::parser &parser, query_arguments &args)
{
    parser.info.author = "Remy Schwab";
    parser.info.version = "1.0.0";
    parser.add_option(args.t, sharg::config{'t', "threads", "Number of threads"});
    parser.add_option(args.verbose, sharg::config{'v', "verbose", "Log verbose output"});
    parser.add_positional_option(args.idx, sharg::config{.description="Path to IBF acid index"});
    parser.add_positional_option(args.regex, sharg::config{.description="Input Regex"});
}

inline void initialise_model_parser(sharg::parser &parser, query_arguments &args)
{
    parser.info.author = "Remy Schwab";
    parser.info.version = "1.0.0";
    parser.add_option(args.t, sharg::config{'t', "threads", "Number of threads"});
    parser.add_option(args.text_length, sharg::config{'l', "length", "Length of text"});
    parser.add_positional_option(args.idx, sharg::config{.description="Path to IBF acid index"});
    parser.add_positional_option(args.regex, sharg::config{.description="Input Regex"});
}

struct inspection_arguments
{
    std::filesystem::path idx{};
};

inline void initialise_inspection_parser(sharg::parser &parser, inspection_arguments &args)
{
    parser.info.author = "Remy Schwab";
    parser.info.version = "1.0.0";
    parser.add_positional_option(args.idx, sharg::config{.description="Path to IBF acid index"});
}
