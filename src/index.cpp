//
// Created by Remy Schwab on 13.09.22.
//

#include "index.h"



void read_input_file_list(std::vector<std::filesystem::path> & sequence_files, std::filesystem::path input_file)
{
    std::ifstream fin{input_file};

    if (!fin.good() || !fin.is_open())
        throw std::runtime_error{"Could not open file " + input_file.string() + " for reading."};

    std::string line;
    while (std::getline(fin, line))
    {
        sequence_files.push_back(line);
    }
}

void create_index_from_filelist(index_arguments &cmd_args)
{
    std::vector<std::filesystem::path> input_bin_files{};
    std::filesystem::path acid_lib = cmd_args.acid_lib;
    read_input_file_list(input_bin_files, acid_lib);

    std::string molecule = cmd_args.molecule;
    size_t seq_count = 0;
    auto hash_adaptor = seqan3::views::kmer_hash(seqan3::ungapped{cmd_args.k});
    size_t bin_count = input_bin_files.size();
    IndexStructure ibf(cmd_args.k, bin_count, cmd_args.bin_size, cmd_args.hash_count, molecule);
    if(molecule == "na")
    {
        for(size_t bin_idx = 0; bin_idx < bin_count; bin_idx++) // For every file in the list...
        {
            record_list<seqan3::dna5_vector> records;
            seq_count += parse_reference_na(input_bin_files[bin_idx], records); // Parse the records from that file...
            populate_bin(ibf, hash_adaptor, records, bin_idx);
        }
    } else if(molecule == "aa")
    {
        for(size_t bin_idx = 0; bin_idx < bin_count; bin_idx++)
        {
            record_list<seqan3::aa27_vector> records;
            seq_count += parse_reference_aa(input_bin_files[bin_idx], records);
            populate_bin(ibf, hash_adaptor, records, bin_idx);
        }
    }
    seqan3::debug_stream << "Indexed " << seq_count << " sequences across " << bin_count << " bins." << std::endl;
    seqan3::debug_stream << "Writing to disk... ";
    std::filesystem::path output_path{cmd_args.ofile+".ibf"};
    store_ibf(ibf, output_path);
    seqan3::debug_stream << std::endl;
}

void drive_index(index_arguments &cmd_args)
{
    std::filesystem::path acid_lib = cmd_args.acid_lib;
    if(acid_lib.extension() == ".lst")
    {
        create_index_from_filelist(cmd_args);
        return;
    }
    if(cmd_args.molecule == "na")
    {
        record_list<seqan3::dna5_vector> records;
        uint32_t bin_count = parse_reference_na(acid_lib, records);
        
        seqan3::debug_stream << "Indexing " << bin_count << " genomes... ";
        IndexStructure ibf = create_index(records, bin_count, cmd_args);
        seqan3::debug_stream << "DONE" << std::endl;

        seqan3::debug_stream << "Writing to disk... ";
        std::filesystem::path output_path{cmd_args.ofile+".ibf"};
        store_ibf(ibf, output_path);
        seqan3::debug_stream << "DONE" << std::endl;
    } else
    {
        record_list<seqan3::aa27_vector> records;
        uint32_t bin_count = parse_reference_aa(acid_lib, records);
        
        seqan3::debug_stream << "Indexing " << bin_count << " genomes... ";
        IndexStructure ibf = create_index(records, bin_count, cmd_args);
        seqan3::debug_stream << "DONE" << std::endl;

        seqan3::debug_stream << "Writing to disk... ";
        std::filesystem::path output_path{cmd_args.ofile+".ibf"};
        store_ibf(ibf, output_path);
        seqan3::debug_stream << "DONE" << std::endl;
    }
}
