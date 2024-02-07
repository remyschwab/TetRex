//
// Created by Remy Schwab on 13.09.22.
//

#include "index_base.h"


std::vector<std::string> read_input_file_list(std::filesystem::path input_file)
{
    std::vector<std::string> sequence_files;
    std::ifstream fin{input_file};

    if (!fin.good() || !fin.is_open())
        throw std::runtime_error{"Could not open file " + input_file.string() + " for reading."};

    std::string line;
    while (std::getline(fin, line))
    {
        sequence_files.emplace_back(line);
    }
    return sequence_files;
}


void create_ibf_dna_index(const index_arguments &cmd_args, const std::vector<std::string> &input_bin_files, uint8_t &reduction)
{
    std::string molecule = cmd_args.molecule;
    TetrexIndex<index_structure::IBF, molecules::nucleotide> ibf(cmd_args.k, cmd_args.bin_size, cmd_args.hash_count, molecule, input_bin_files, reduction);
    ibf.populate_index();
    seqan3::debug_stream << "Writing to disk... ";
    std::filesystem::path output_path{cmd_args.ofile+".ibf"};
    store_ibf(ibf, output_path);
    seqan3::debug_stream << "DONE" << std::endl;
}


void create_hibf_dna_index(const index_arguments &cmd_args, const std::vector<std::string> &input_bin_files, uint8_t &reduction)
{
    std::string molecule = cmd_args.molecule;
    TetrexIndex<index_structure::HIBF, molecules::nucleotide> ibf(cmd_args.k, cmd_args.bin_size, cmd_args.hash_count, molecule, input_bin_files, reduction);
    ibf.populate_index();
    seqan3::debug_stream << "Writing to disk... ";
    std::filesystem::path output_path{cmd_args.ofile+".ibf"};
    store_ibf(ibf, output_path);
    seqan3::debug_stream << "DONE" << std::endl;
}


void create_ibf_aa_index(const index_arguments &cmd_args, const std::vector<std::string> &input_bin_files, uint8_t &reduction)
{
    std::string molecule = cmd_args.molecule;
    TetrexIndex<index_structure::IBF, molecules::peptide> ibf(cmd_args.k, cmd_args.bin_size, cmd_args.hash_count, molecule, input_bin_files, reduction);
    ibf.populate_index();
    seqan3::debug_stream << "Writing to disk... ";
    std::filesystem::path output_path{cmd_args.ofile+".ibf"};
    store_ibf(ibf, output_path);
    seqan3::debug_stream << "DONE" << std::endl;
}


void create_hibf_aa_index(const index_arguments &cmd_args, const std::vector<std::string> &input_bin_files, uint8_t &reduction)
{
    std::string molecule = cmd_args.molecule;
    TetrexIndex<index_structure::HIBF, molecules::peptide> ibf(cmd_args.k, cmd_args.bin_size, cmd_args.hash_count, molecule, input_bin_files, reduction);
    ibf.populate_index();
    seqan3::debug_stream << "Writing to disk... ";
    std::filesystem::path output_path{cmd_args.ofile+".ibf"};
    store_ibf(ibf, output_path);
    seqan3::debug_stream << "DONE" << std::endl;
}


void drive_index(const index_arguments &cmd_args)
{
    std::vector<std::string> input_bin_files;
    for (auto file : cmd_args.acid_libs)
    {
        if (file.extension() == ".lst")
        {
            for (auto f : read_input_file_list(file))
                input_bin_files.push_back(f);
        } 
        else
        {
            input_bin_files.push_back(std::filesystem::absolute(file));
        }
    }

    bool flavor_test = (cmd_args.ibf == "ibf");
    bool alpha_test = cmd_args.molecule == "dna";

    robin_hood::unordered_map<std::string, uint8_t> reduction_map = {{"murphy", 1u},{"li", 2u},{"None", 0u}};
    uint8_t reduction = reduction_map[cmd_args.reduction];

    if(flavor_test)
    {
        if(alpha_test)
        {
            create_ibf_dna_index(cmd_args, input_bin_files, reduction);
        }
        else
        {
            create_ibf_aa_index(cmd_args, input_bin_files, reduction);
        }
    }
    else
    {
        if(alpha_test)
        {
            create_hibf_dna_index(cmd_args, input_bin_files, reduction);
        }
        else
        {
            create_hibf_aa_index(cmd_args, input_bin_files, reduction);
        }
    }
}
