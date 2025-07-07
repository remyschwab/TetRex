//
// Created by Remy Schwab on 12.09.22.
//

#pragma once

#include <sharg/all.hpp>


struct index_arguments
{
    uint8_t k = 6;
    uint8_t t = 1;
    float fpr = 0.05;
    bool idx{false};
    std::string ibf = idx ? "ibf" : "hibf";
    std::string ofile;
    std::vector<std::filesystem::path> acid_libs{};
    uint8_t hash_count = 3;
    bool dna{false};
    std::string molecule = dna ? "na" : "aa";
    std::string reduction = "None";
};

inline void initialise_index_parser(sharg::parser &parser, index_arguments &args)
{
    parser.info.author = "Remy Schwab";
    parser.info.version = "1.0.0";
    parser.add_option(args.k, sharg::config{'k', "ksize", "size of kmers"});
    parser.add_option(args.fpr, sharg::config{'p', "fpr", "Bloom Filter False Positive Rate"});
    parser.add_option(args.hash_count, sharg::config{'c', "hash_count", "Number of hash functions. NOTE: MORE THAN 4 IS SLOW"});
    parser.add_option(args.t, sharg::config{'t', "threads", "Number of threads"});
    parser.add_flag(args.dna, sharg::config{'n', "nucleic_acid", "Index a library of Nucleic Acids (Default is Amino Acid)"});
    parser.add_flag(args.idx, sharg::config{'i', "ibf", "Index using IBF (Default is HIBF)"});
    parser.add_option(args.reduction, sharg::config{.short_id='r', .long_id="reduce", .description="Use reduced AA alphabet (Murphy or Li)", .validator=sharg::value_list_validator{"murphy","li"}});
    parser.add_positional_option(args.ofile, sharg::config{.description="Name of index on disk"});
    parser.add_positional_option(args.acid_libs, sharg::config{.description="Nucleic or Amino Acid library to indexed", .validator=sharg::input_file_validator{{"lst","fa", "fa.gz","fasta", "fasta.gz", "fna", "fna.gz"}}});
}

struct query_arguments
{
    uint8_t t = 1;
    bool verbose{false};
    bool draw{false};
    bool read_file{false};
    bool conjunction{false};
    bool augment{false};
    int text_length;
    size_t complexity = 10;
    std::filesystem::path idx{};
    std::string regex = "-";
    std::string query;
};

inline void initialise_query_parser(sharg::parser &parser, query_arguments &args)
{
    parser.info.author = "Remy Schwab";
    parser.info.version = "1.0.0";
    parser.add_option(args.t, sharg::config{'t', "threads", "Number of threads"});
    parser.add_option(args.complexity, sharg::config{'n', "threshold", "Node count to trigger augmentation"});
    parser.add_flag(args.draw, sharg::config{'d', "draw", "Write Graph Viz file to disk"});
    parser.add_flag(args.verbose, sharg::config{'v', "verbose", "Log verbose output"});
    parser.add_flag(args.read_file, sharg::config{'f', "file", "Interpret last argument as a file containing RegEx"});
    parser.add_flag(args.conjunction, sharg::config{'c', "conj", "Search multiple queries delimited with ':'"});
    parser.add_flag(args.augment, sharg::config{'a', "augment", "Skip over high complexity regions of RegEx"});
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

struct dindex_arguments
{
    size_t min = 3;
    size_t max = 21;
    size_t pad = 1;
    bool dna{false};
    std::string molecule = dna ? "na" : "aa";
    std::string reduction = "None";
    uint8_t hash_count = 3;
    float fpr = 0.05;
    bool idx{false};
    std::string ibf = idx ? "ibf" : "hibf";
    std::string ofile;
    std::vector<std::filesystem::path> acid_libs{};
};

inline void initialize_dgram_parser(sharg::parser &parser, dindex_arguments &args)
{
    parser.info.author = "Remy Schwab";
    parser.info.version = "1.0.0";
    parser.add_flag(args.dna, sharg::config{'n', "nucleic_acid", "Index a library of Nucleic Acids (Default is Amino Acid)"});
    parser.add_flag(args.idx, sharg::config{'i', "ibf", "Index using IBF (Default is HIBF)"});
    parser.add_option(args.min, sharg::config{'l', "lower", "Lower bound gap size"});
    parser.add_option(args.max, sharg::config{'u', "upper", "Upper bound gap size"});
    // parser.add_option(args.pad, sharg::config{'p', "padding", "Number of characters on either side of gap"});
    parser.add_positional_option(args.ofile, sharg::config{.description="Name of index on disk"});
    parser.add_positional_option(args.acid_libs, sharg::config{.description="Nucleic or Amino Acid library to indexed", .validator=sharg::input_file_validator{{"lst","fa", "fa.gz","fasta", "fasta.gz", "fna", "fna.gz"}}});
}

