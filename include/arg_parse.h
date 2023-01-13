//
// Created by Remy Schwab on 12.09.22.
//

#ifndef KBIOREG_ARG_PARSE_H
#define KBIOREG_ARG_PARSE_H

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/argument_parser/validators.hpp>


struct index_arguments
{
    uint8_t k = 3;
    uint8_t t = 1;
    std::string ofile;
    std::filesystem::path acid_lib{};
    uint32_t bin_size = 10000;
    uint8_t hash_count = 2;
    std::string molecule;
};

struct query_arguments
{
    uint8_t t = 1;
    std::filesystem::path idx{};
    std::string regex;
    std::string query;
    std::filesystem::path graph = "graph";
};

struct benchmark_arguments
{
    uint8_t k = 3; //all
    uint8_t t = 1; //all
    uint8_t w = 100; //script
    std::string s = "400000"; //script
    std::string b = "64"; //script
    std::string p = "20"; //script
    std::string regex;
    std::string query;
    std::filesystem::path idx = "benchmark_idx.ibf"; //query 
    std::string ofile = "benchmark_idx"; //index
    std::filesystem::path acid_lib= "64/bins/all_bins.fa"; //index
    uint64_t bin_size = 8192; //index
    uint8_t hash_count = 2; //index
    std::string molecule; //index
    std::filesystem::path graph = "graph";
};
    
inline void initialise_index_parser(seqan3::argument_parser &parser, index_arguments &args)
{
    parser.info.author = "Remy Schwab";
    parser.info.version = "1.0.0";
    parser.add_option(args.k, 'k', "ksize", "size of kmers");
    parser.add_option(args.bin_size, 's', "bin_size", "Size of bins");
    parser.add_option(args.hash_count, 'c', "hash_count", "Number of hash functions. NOTE: MORE THAN 4 IS SLOW");
    parser.add_option(args.t, 't', "threads", "Number of threads");
    parser.add_option(args.molecule, 'm', "molecule", "Molecule type of library", seqan3::option_spec::required,
                                seqan3::value_list_validator{"na", "aa"});
    parser.add_option(args.ofile, 'o', "ofile", "Name of index on disk");
    parser.add_positional_option(args.acid_lib, "Nucleic or Amino Acid library to indexed",
                                seqan3::input_file_validator{{"fq","fastq","fa","fasta", "fna"}});
}

inline void initialise_query_parser(seqan3::argument_parser &parser, query_arguments &args)
{
    parser.info.author = "Remy Schwab";
    parser.info.version = "1.0.0";
    parser.add_option(args.t, 't', "threads", "Number of threads");
    parser.add_positional_option(args.idx, "Path to IBF acid index");
    parser.add_positional_option(args.regex, "Input Regex in reverse polish notation");
    parser.add_positional_option(args.graph, "path for .fot file of the automaton");
}

inline void initialise_benchmark_parser(seqan3::argument_parser &parser, benchmark_arguments &args)
{
    parser.info.author = "Vincent Musch";
    parser.info.version = "1.0.0";
    parser.add_option(args.k, 'k', "ksize", "size of kmers");
    parser.add_option(args.t, 't', "threads", "Number of threads");
    parser.add_option(args.w, 'w', "words", "Number of all random words vor the fake genome");
    parser.add_option(args.s, 'j', "gsize", "size of the fake genome");
    parser.add_option(args.b, 'b', "bin", "number of bins for the ibf");
    parser.add_option(args.p, 'p', "percent", "perzent of hidden words in fake genome");
    parser.add_option(args.bin_size, 's', "bin_size", "Size of bins");
    parser.add_option(args.hash_count, 'c', "hash_count", "Number of hash functions. NOTE: MORE THAN 4 IS SLOW");
    parser.add_option(args.molecule, 'm', "molecule", "Molecule type of library", seqan3::option_spec::required,
                                seqan3::value_list_validator{"na", "aa"});
    parser.add_option(args.ofile, 'o', "ofile", "Name of index on disk");
    parser.add_option(args.graph, 'g', "graph", "Name of dot file");
    parser.add_positional_option(args.regex, "Input Regex in reverse polish notation");
}

#endif //KBIOREG_ARG_PARSE_H
