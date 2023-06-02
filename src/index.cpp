//
// Created by Remy Schwab on 13.09.22.
//

#include "index.h"
KSEQ_INIT(gzFile, gzread)


void convertStringToVec(const std::string &kmer, const IndexStructure &ibf, std::vector<unsigned char> &digit_vec)
{
    for(size_t i = 0; i < kmer.length(); ++i)
        digit_vec[i] = ibf.aamap_[kmer[i]];
}


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


void decompose_nucleotide_record(std::string_view record_seq, IndexStructure &ibf, const size_t &bin_id)
{
    uint64_t initial_encoding = encode_dna(record_seq.substr(0,ibf.k_)); // Encode forward
    uint64_t reverse_complement = revComplement(initial_encoding, ibf.k_); // Compute the reverse compelement
    ibf.set_stores(initial_encoding, reverse_complement); // Remember both strands
    ibf.emplace((initial_encoding <= reverse_complement ? initial_encoding : reverse_complement), bin_id); // MinHash
    for(size_t i = ibf.k_; i < record_seq.length(); ++i)
    {
        auto symbol = record_seq[i];
        ibf.rollover_nuc_hash(symbol, bin_id);
    }
}


void decompose_peptide_record(const std::vector<unsigned char> &int_seq, const size_t &begin, IndexStructure &ibf)
{
    size_t res1, res2, res3, res4;
    size_t numbElements = ibf.k_;
    size_t end = begin+ibf.k_;
    switch(numbElements)
    {
        case 6:
            res1 = int_seq[begin+0]*ibf.powers_[0];
            res2 = int_seq[begin+1]*ibf.powers_[1];
            res3 = int_seq[begin+2]*ibf.powers_[2];
            res4 = int_seq[begin+3]*ibf.powers_[3];
            res1 += int_seq[begin+4]*ibf.powers_[4];
            res2 += int_seq[begin+5]*ibf.powers_[5];
            ibf.forward_store_ = res1 + res2 + res3 + res4;
            break;
        case 7:
            res1 = int_seq[begin+0]*ibf.powers_[0];
            res2 = int_seq[begin+1]*ibf.powers_[1];
            res3 = int_seq[begin+2]*ibf.powers_[2];
            res4 = int_seq[begin+3]*ibf.powers_[3];
            res1 += int_seq[begin+4]*ibf.powers_[4];
            res2 += int_seq[begin+5]*ibf.powers_[5];
            res3 += int_seq[begin+6]*ibf.powers_[6];
            ibf.forward_store_ = res1 + res2 + res3 + res4;
            break;
        case 10:
            res1 = int_seq[begin+0]*ibf.powers_[0];
            res2 = int_seq[begin+1]*ibf.powers_[1];
            res3 = int_seq[begin+2]*ibf.powers_[2];
            res4 = int_seq[begin+3]*ibf.powers_[3];
            res1 += int_seq[begin+4]*ibf.powers_[4];
            res2 += int_seq[begin+5]*ibf.powers_[5];
            res3 += int_seq[begin+6]*ibf.powers_[6];
            res4 += int_seq[begin+7]*ibf.powers_[7];
            res1 += int_seq[begin+8]*ibf.powers_[8];
            res2 += int_seq[begin+9]*ibf.powers_[9];
            ibf.forward_store_ = res1 + res2 + res3 + res4;
            break;
        case 14:
            res1 = int_seq[begin+0]*ibf.powers_[0];
            res2 = int_seq[begin+1]*ibf.powers_[1];
            res3 = int_seq[begin+2]*ibf.powers_[2];
            res4 = int_seq[begin+3]*ibf.powers_[3];
            res1 += int_seq[begin+4]*ibf.powers_[4];
            res2 += int_seq[begin+5]*ibf.powers_[5];
            res3 += int_seq[begin+6]*ibf.powers_[6];
            res4 += int_seq[begin+7]*ibf.powers_[7];
            res1 += int_seq[begin+8]*ibf.powers_[8];
            res2 += int_seq[begin+9]*ibf.powers_[9];
            res3 += int_seq[begin+10]*ibf.powers_[10];
            res4 += int_seq[begin+11]*ibf.powers_[11];
            res1 += int_seq[begin+12]*ibf.powers_[12];
            res2 += int_seq[begin+13]*ibf.powers_[13];
            ibf.forward_store_ = res1 + res2 + res3 + res4;
            break;
        default:
            for(size_t i = begin; i < end; i++)
                ibf.forward_store_ += int_seq[i]*ibf.powers_[i-begin];
            break;
    }
}


void create_index(const index_arguments &cmd_args, const std::vector<std::string> &input_bin_files)
{
    gzFile handle;
    kseq_t *record;
    int status;

    std::string molecule = cmd_args.molecule;
    size_t seq_count = 0;
    size_t bin_count = input_bin_files.size();
    IndexStructure ibf(cmd_args.k, bin_count, cmd_args.bin_size, cmd_args.hash_count, molecule, input_bin_files, cmd_args.reduction);
    DGramIndex dibf(bin_count, cmd_args.dbin_size, cmd_args.hash_count);
    ibf.set_lib_paths(input_bin_files);

    for(size_t i = 0; i < input_bin_files.size(); ++i)
    {
        handle = gzopen(input_bin_files[i].c_str(), "r");
        record = kseq_init(handle);
        while ((status = kseq_read(record)) >= 0)
        {
            std::string_view record_view = record->seq.s;
            if(record_view.length() < ibf.k_)
            {
                seqan3::debug_stream << "RECORD TOO SHORT " << record->comment.s << std::endl;
                continue;
            }
            seq_count++;
            dibf.init_tracker();
            if(ibf.molecule_ == "na")
            {
                decompose_nucleotide_record(record_view, ibf, i);
            }
            else
            {
                std::vector<uint8_t> peptide_nums(record_view.length());
                for(size_t o = 0; o < record_view.length(); ++o)
                {
                    dibf.track_record(o, record_view[o], i);
                    peptide_nums[o] = ibf.aamap_[record_view[o]];
                }
                for(size_t o = 0; o < peptide_nums.size()-ibf.k_+1; ++o)
                {
                    decompose_peptide_record(peptide_nums, o, ibf);
                    ibf.emplace(ibf.forward_store_, i);
                }
            }
        }
        kseq_destroy(record);
        gzclose(handle);
        // auto && dibf_ref = dibf.getDIBF();
        // auto agent = dibf_ref.membership_agent();
        // agent.bulk_contains();

    }

    seqan3::debug_stream << "Indexed " << seq_count << " sequences across " << bin_count << " bins." << std::endl;
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
    create_index(cmd_args, input_bin_files);
}
