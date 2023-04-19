//
// Created by Remy Schwab on 13.09.22.
//

#include "index.h"
KSEQ_INIT(gzFile, gzread)


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


void decompose_record(std::string_view record_seq, IndexStructure ibf, const size_t &bin_id)
{
    uint64_t initial_encoding = encode_dna(record_seq.substr(0,ibf.k_)); // Encode forward
    uint64_t reverse_complement = revComplement(initial_encoding, ibf.k_); // Compute the reverse compelement
    ibf.set_stores(initial_encoding, reverse_complement); // Remember both strands
    ibf.emplace(( initial_encoding <= reverse_complement ? initial_encoding : reverse_complement ), bin_id); // MinHash
    
    for(size_t i = ibf.k_; i < record_seq.length(); ++i)
    {
        auto symbol = record_seq[i];
        ibf.rollover_hash(symbol, bin_id);
        seqan3::debug_stream << symbol << std::endl;
    }
}


void read_with_kseq(const index_arguments &cmd_args, const std::vector<std::string> &input_bin_files)
{
    gzFile handle;
    kseq_t *record;
    int status;

    uint8_t k = cmd_args.k;
    size_t canon_encoded_kmer;
    std::string molecule = cmd_args.molecule;
    size_t seq_count = 0;
    size_t bin_count = input_bin_files.size();
    IndexStructure ibf(cmd_args.k, bin_count, cmd_args.bin_size, cmd_args.hash_count, molecule, input_bin_files);
    ibf.set_lib_paths(input_bin_files);

    for(size_t i = 0; i < input_bin_files.size(); ++i)
    {
        handle = gzopen(input_bin_files[i].c_str(), "r");
        record = kseq_init(handle);
        while ((status = kseq_read(record)) >= 0) {
            seq_count++;
            std::string_view record_view = record->seq.s;
            decompose_record(record_view, ibf, i);
        }
        kseq_destroy(record);
        gzclose(handle);
    }
    seqan3::debug_stream << "Indexed " << seq_count << " sequences across " << bin_count << " bins." << std::endl;
    seqan3::debug_stream << "Writing to disk... ";
    std::filesystem::path output_path{cmd_args.ofile+".ibf"};
    store_ibf(ibf, output_path);
    seqan3::debug_stream << "DONE" << std::endl;
}


// void create_index_from_filelist(const index_arguments &cmd_args, const std::vector<std::string> &input_bin_files)
// {
//     std::string molecule = cmd_args.molecule;
//     size_t seq_count = 0;
//     auto hash_adaptor = seqan3::views::kmer_hash(seqan3::ungapped{cmd_args.k});
//     size_t bin_count = input_bin_files.size();
//     IndexStructure ibf(cmd_args.k, bin_count, cmd_args.bin_size, cmd_args.hash_count, molecule, input_bin_files);
//     ibf.set_lib_paths(input_bin_files);
//     if(molecule == "na")
//     {
//         for(size_t bin_idx = 0; bin_idx < bin_count; bin_idx++) // For every file in the list...
//         {
//             record_list<seqan3::dna5_vector> records;
//             std::filesystem::path bin_file{input_bin_files[bin_idx]};
//             seq_count += parse_reference_na(bin_file, records); // Parse the records from that file...
//             populate_bin(ibf, hash_adaptor, records, bin_idx);
//         }
//     } else if(molecule == "aa")
//     {
//         for(size_t bin_idx = 0; bin_idx < bin_count; bin_idx++)
//         {
//             record_list<seqan3::aa27_vector> records;
//             std::filesystem::path bin_file{input_bin_files[bin_idx]};
//             seq_count += parse_reference_aa(bin_file, records);
//             populate_bin(ibf, hash_adaptor, records, bin_idx);
//         }
//     }
//     seqan3::debug_stream << "Indexed " << seq_count << " sequences across " << bin_count << " bins." << std::endl;
//     seqan3::debug_stream << "Writing to disk... ";
//     std::filesystem::path output_path{cmd_args.ofile+".ibf"};
//     store_ibf(ibf, output_path);
//     seqan3::debug_stream << "DONE" << std::endl;
// }


void drive_index(const index_arguments &cmd_args)
{
    std::vector<std::string> input_bin_files;
    for (auto file : cmd_args.acid_libs) {
        if (file.extension() == ".lst") {
            for (auto f : read_input_file_list(file)) {
                input_bin_files.push_back(f);
            }
        } else {
            input_bin_files.push_back(file);
        }
    }
    // create_index_from_filelist(cmd_args, input_bin_files);
    read_with_kseq(cmd_args, input_bin_files);
}
