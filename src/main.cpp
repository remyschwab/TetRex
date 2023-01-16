#include "utils.h"
#include "arg_parse.h"
#include "index.h"
#include "query.h"
#include <fstream>
#include <omp.h>
#include <stdlib.h>
#include <algorithm>


void run_index(seqan3::argument_parser &parser)
{
    index_arguments cmd_args{};
    initialise_index_parser(parser, cmd_args);
    try
    {
        parser.parse();
    }
    catch (seqan3::argument_parser_error const & ext)
    {
        seqan3::debug_stream << "[Indexing Parser Error] " << ext.what() << "\n";
        return;
    }
    drive_index(cmd_args);
}

void run_query(seqan3::argument_parser &parser)
{
    // Parse Arguments
    query_arguments cmd_args{};
    initialise_query_parser(parser, cmd_args);
    try
    {
        parser.parse();
    }
    catch (seqan3::argument_parser_error const & ext)
    {
        seqan3::debug_stream << "[Error kBioReg Query module " << ext.what() << "\n";
        return;
    }
    cmd_args.query = translate(cmd_args.regex);
    // drive_query(cmd_args);
}

void run_benchmark(seqan3::argument_parser &parser)
{
    benchmark_arguments cmd_args{};
    initialise_benchmark_parser(parser, cmd_args);
    try
    {
        parser.parse(); // trigger command line parsing
    }
    catch (seqan3::argument_parser_error const & ext) // catch user errors
    {
        seqan3::debug_stream << "[Error] " << ext.what() << "\n"; // customise your error message
        return;
    }

    cmd_args.query = translate(cmd_args.regex);
    // cmd_args.query = cmd_args.regex;

    seqan3::debug_stream << "=======Start Benchmark========" << std::endl;
    seqan3::debug_stream << "=======Make Genome========" << std::endl;
    uint8_t number = cmd_args.w; //number of words
    //std::string query = cmd_args.query; //regex
    std::string bashC ="../bash_scripts/simulate.sh "+cmd_args.s+" "+cmd_args.b+" "+cmd_args.p;
    
    // Postfix to Thompson NFA
    seqan3::debug_stream <<"   - Constructing Thompson NFA from RegEx... ";
    State* nfa = post2nfaE(cmd_args.query);
    seqan3::debug_stream << "DONE" << std::endl;

    // Make a file
    seqan3::debug_stream << "   - Write words.txt ";
    std::fstream file;
    file.open("words.txt", std::ios::out);
    while(number > 0)
    {
        std::string word = getRandomWord(nfa);
        file << word <<"\n" ;   
        number--;
    }
    file.close();
    seqan3::debug_stream << "DONE" << std::endl;

    // Execute shell wrapper to generate syn genome
    seqan3::debug_stream << "   - start script ";
    int occ = system(bashC.c_str());
    seqan3::debug_stream << "DONE" << std::endl;

    //write parser
    seqan3::debug_stream << "   - write parser";

    index_arguments index{};
    index.k = cmd_args.k;
    index.t = cmd_args.t;
    index.ofile = cmd_args.ofile;
    index.acid_lib = cmd_args.b + "/bins/all_bins.fa";
    index.molecule = cmd_args.molecule;
    index.bin_size = cmd_args.bin_size;
    index.hash_count = cmd_args.hash_count;

    query_arguments query{};
    query.idx = cmd_args.idx;
    query.regex = cmd_args.regex;
    query.query = cmd_args.query;
    query.graph = cmd_args.graph;

    seqan3::debug_stream << "DONE" << std::endl;


    seqan3::debug_stream << "=======Genome Done========" << std::endl;
    
    seqan3::debug_stream << "=======Start Search========" << std::endl;
    //seqan3::argument_parser query_parser = {};
    //indexing ./bin/kbioreg benchmark -k 4 -w 30 -p 10 words.txt "AC.G+.T."
    std::fstream benchmark_table;
    benchmark_table.open("benchmark_table.txt", std::ios_base::app);
    double t1,t2,t3,t4;
    t1 = omp_get_wtime();
    drive_index(index);
    t2 = omp_get_wtime();
    benchmark_table<< t2-t1<<",";   //timeIndex
    seqan3::debug_stream << "[INDEXING TIME]: " << t2-t1 << std::endl;

    // Query
    int hitsNr = 0;
    std::fstream file_kbioreg;
    file_kbioreg.clear();
    file_kbioreg.open("kbioreg_search.txt", std::ios::out);
    std::string str = cmd_args.regex;
    std::regex reg(str);
    size_t bins = std::stoi(cmd_args.b);
    std::string binString = "bin_";
    for(size_t i = 1; i < cmd_args.b.size(); i++)
    {
        binString = binString + "0";
    }

    t1 = omp_get_wtime();
    bitvector hits = drive_query_benchmark(query, benchmark_table);
    t2 = omp_get_wtime();
    benchmark_table<<t2-t1<<",";

    seqan3::debug_stream << "[KBIOREG QUERY TIME]: " << t2-t1 << std::endl;

    size_t j = 10;
    t3 = omp_get_wtime();
    for(size_t i = 0; i < bins; i++)
    {
        if(i >= j)
        {
            binString.pop_back();
            j = j * 10;
        } 
        if(hits[i])
        {
            std::string filepath = cmd_args.b+"/bins/"+binString+std::to_string(i)+".fa";
            hitsNr += matches(stream_as_string(filepath), reg, file_kbioreg);
        }
    }
    
    t4 = omp_get_wtime();
    seqan3::debug_stream << "[KBIOREG STDBIN SEARCH TIME]: " << t4-t3 << std::endl;
    seqan3::debug_stream << "[KBIOREG TOTAL SEARCH TIME]: " << t4-t1 << std::endl;
    file_kbioreg<<"[KBIOREG QUERY TIME]: " << t4-t1<<"\n";
    file_kbioreg<<"[KBIOREG Genome Size]: " << cmd_args.s<<"\n";
    file_kbioreg<<"[KBIOREG Bins]: " << cmd_args.b<<"\n";
    file_kbioreg<<"[KBIOREG Hits]: " << hitsNr<<"\n";
    benchmark_table<< t4-t3<<",";   //time for loop
    benchmark_table<< t4-t1<<",";   //timeGlobal
    benchmark_table<< cmd_args.s<<","; //genom size
    benchmark_table<< hitsNr<<",";    //hitsQ
    benchmark_table<< cmd_args.b<<",";    //bins

    file_kbioreg.close();
    
    // Standard library RegEx Searc
    std::fstream file_std;
    file_std.open("standard_search.txt", std::ios::out);
    file_std.clear();
    t1 = omp_get_wtime();
    hitsNr = matches(stream_as_string(cmd_args.b+"/bins/all_bins.fa"), reg, file_std);
    t2 = omp_get_wtime();
    file_std<<"[STDREGEX QUERY TIME]: " << t2-t1<<"\n";
    file_std<<"[STDREGEX Genom Size]: " << cmd_args.s<<"\n";
    file_std<<"[STDREGEX Bins]: " << cmd_args.b<<"\n";
    file_std<<"[STDREGEX Hits]: " << hitsNr<<"\n";
    file_std.close();

    benchmark_table<< t2-t1<<",";   //timeS
    benchmark_table<< hitsNr<<"\n";    //hitsQ

    benchmark_table.close();
    seqan3::debug_stream<<"[STDREGEX TIME]: "<< t2-t1 << std::endl;
}

int main(int argc, char *argv[])
{
    seqan3::argument_parser top_level_parser{"kbioreg", argc, argv,
                                             seqan3::update_notifications::off,
                                             {"index", "query", "benchmark"}};
    top_level_parser.info.description.push_back("Index a NA|AA FASTA library or search a regular expression.");
    try
    {
        top_level_parser.parse(); // trigger command line parsing
    }
    catch (seqan3::argument_parser_error const & ext) // catch user errors
    {
        seqan3::debug_stream << "[Error] " << ext.what() << "\n"; // customise your error message
        return -1;
    }

    seqan3::argument_parser & sub_parser = top_level_parser.get_sub_parser(); // hold a reference to the sub_parser
    if (sub_parser.info.app_name == std::string_view{"kbioreg-index"})
        run_index(sub_parser);
    else if (sub_parser.info.app_name == std::string_view{"kbioreg-query"})
        run_query(sub_parser);
    else if (sub_parser.info.app_name == std::string_view{"kbioreg-benchmark"})
        run_benchmark(sub_parser);
    else
        std::cout << "Unhandled subparser named " << sub_parser.info.app_name << '\n';
  return 0;
}
