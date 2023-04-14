//
// Created by Remy Schwab on 12.09.22.
//

#pragma once

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/argument_parser/validators.hpp>


struct index_arguments
{
    uint8_t k = 3;
    uint8_t t = 1;
    std::string ofile;
    std::vector<std::filesystem::path> acid_libs{};
    size_t bin_size = 10000;
    uint8_t hash_count = 3;
    std::string reduction;
    std::string molecule;
};

inline void initialise_index_parser(seqan3::argument_parser &parser, index_arguments &args)
{
    parser.info.author = "Remy Schwab";
    parser.info.version = "1.0.0";
    parser.add_option(args.k, 'k', "ksize", "size of kmers");
    parser.add_option(args.bin_size, 's', "bin_size", "Size of bins");
    parser.add_option(args.hash_count, 'c', "hash_count", "Number of hash functions. NOTE: MORE THAN 4 IS SLOW");
    parser.add_option(args.t, 't', "threads", "Number of threads");
    parser.add_option(args.molecule, 'm', "molecule", "Molecule type of library", seqan3::option_spec::required, seqan3::value_list_validator{"na", "aa"});
    parser.add_option(args.reduction, 'r', "reduce", "Use reduced AA alphabet (Murphy or Li)");                               
    parser.add_option(args.ofile, 'o', "ofile", "Name of index on disk");
    parser.add_positional_option(args.acid_libs, "Nucleic or Amino Acid library to indexed",
                                seqan3::input_file_validator{{"lst","fa", "fa.gz","fasta", "fasta.gz", "fna", "fna.gz"}});
}

struct query_arguments
{
    uint8_t t = 1;
    int text_length;
    std::filesystem::path graph{};
    std::filesystem::path idx{};
    std::string regex;
    std::string query;
};

inline void initialise_query_parser(seqan3::argument_parser &parser, query_arguments &args)
{
    parser.info.author = "Remy Schwab";
    parser.info.version = "1.0.0";
    parser.add_option(args.t, 't', "threads", "Number of threads");
    parser.add_option(args.graph, 'd', "dot", "Path to dot file");
    parser.add_positional_option(args.idx, "Path to IBF acid index");
    parser.add_positional_option(args.regex, "Input Regex in reverse polish notation");
}

inline void initialise_model_parser(seqan3::argument_parser &parser, query_arguments &args)
{
    parser.info.author = "Remy Schwab";
    parser.info.version = "1.0.0";
    parser.add_option(args.t, 't', "threads", "Number of threads");
    parser.add_option(args.text_length, 'l', "length", "Length of text");
    parser.add_option(args.graph, 'd', "dot", "Path to dot file");
    parser.add_positional_option(args.idx, "Path to IBF acid index");
    parser.add_positional_option(args.regex, "Input Regex in reverse polish notation");
}

struct inspection_arguments
{
    std::filesystem::path idx{};
};

inline void initialise_inspection_parser(seqan3::argument_parser &parser, inspection_arguments &args)
{
    parser.info.author = "Remy Schwab";
    parser.info.version = "1.0.0";
    parser.add_positional_option(args.idx, "Path to IBF acid index");
}
