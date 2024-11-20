#include "utils.h"
#include "arg_parse.h"
#include "index_base.h"
#include "query.h"
#include "inspect_idx.h"
#include <fstream>
#include <omp.h>
#include <stdlib.h>
#include <algorithm>


void run_index(sharg::parser &parser)
{
    index_arguments cmd_args{};
    initialise_index_parser(parser, cmd_args);
    try
    {
        parser.parse();
    }
    catch (sharg::parser_error const & ext)
    {
        seqan3::debug_stream << "[Indexing Parser Error] " << ext.what() << "\n";
        return;
    }
    if(cmd_args.molecule == "aa" && cmd_args.k > 12)
    {
        seqan3::debug_stream << "[Indexing Parser Error] Max kmer size for amino acids is 12" << "\n";
        return;
    }
    drive_index(cmd_args);
}

void run_query(sharg::parser &parser)
{
    // Parse Arguments
    query_arguments cmd_args{};
    initialise_query_parser(parser, cmd_args);
    try
    {
        parser.parse();
    }
    catch (sharg::parser_error const & ext)
    {
        seqan3::debug_stream << "[Error TetRex Query module " << ext.what() << "\n";
        return;
    }
    try
    {
        drive_query(cmd_args, false);
    }
    catch(const cereal::Exception& e)
    {
        std::cerr << "Filepath to (H)IBF Index not valid" << std::endl;
        std::cerr << e.what() << '\n';
    }
}

void run_inspection(sharg::parser &parser)
{
    // Parse Arguments
    inspection_arguments cmd_args{};
    initialise_inspection_parser(parser, cmd_args);
    try
    {
        parser.parse();
    }
    catch (sharg::parser_error const & ext)
    {
        seqan3::debug_stream << "[Error TetRex Index Inspection module " << ext.what() << "\n";
        return;
    }
    drive_inspection(cmd_args);
}

void run_model(sharg::parser &parser)
{
    // Parse Arguments
    query_arguments cmd_args{};
    initialise_model_parser(parser, cmd_args);
    try
    {
        parser.parse();
    }
    catch (sharg::parser_error const & ext)
    {
        seqan3::debug_stream << "[Error TetRex Discovery Modeling module " << ext.what() << "\n";
        return;
    }
    drive_query(cmd_args, true);
}

int main(int argc, char *argv[])
{
    sharg::parser top_level_parser{"tetrex", argc, argv,
                                             sharg::update_notifications::off,
                                             {"index", "query", "inspect"}};
    top_level_parser.info.description.push_back("Index a NA|AA FASTA library or search a regular expression.");
    try
    {
        top_level_parser.parse(); // trigger command line parsing
    }
    catch (sharg::parser_error const & ext) // catch user errors
    {
        seqan3::debug_stream << "[Error] " << ext.what() << "\n"; // customise your error message
        return -1;
    }

    sharg::parser &sub_parser = top_level_parser.get_sub_parser(); // hold a reference to the sub_parser
    if (sub_parser.info.app_name == std::string_view{"tetrex-index"})
        run_index(sub_parser);
    else if (sub_parser.info.app_name == std::string_view{"tetrex-query"})
        run_query(sub_parser);
    else if (sub_parser.info.app_name == std::string_view{"tetrex-inspect"})
        run_inspection(sub_parser);
    // else if (sub_parser.info.app_name == std::string_view{"tetrex-model"})
    //     run_model(sub_parser);
    else
        std::cout << "Unhandled subparser named " << sub_parser.info.app_name << '\n';
    return 0;
}
