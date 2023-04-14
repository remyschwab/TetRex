//
// Created by Remy Schwab on 13.09.22.
//

#pragma once


#include "utils.h"
#include "arg_parse.h"
#include "kseq.h"
#include "molecule_encodings.h"

#include <zlib.h>

#include <cereal/types/vector.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/views/kmer_hash.hpp>
#include <seqan3/alphabet/nucleotide/all.hpp>



namespace index_structure
{

    using ibf = seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed>;
    using ibf_compressed = seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed>;

    template <typename return_t, typename input_t>
    concept compressible_from = (std::same_as<return_t, ibf_compressed> && std::same_as<input_t, ibf>);

} // namespace index_structure

class IndexStructure
{

private:
    size_t bin_count_;
    size_t bin_size_{};
    uint8_t hash_count_{};
    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf_;

public:

    uint8_t k_;
    std::string molecule_;
    std::vector<std::string> acid_libs_;

    IndexStructure() = default;

    explicit IndexStructure(uint8_t k,
                            size_t bc,
                            size_t bs,
                            uint8_t hc,
                            std::string molecule,
                            std::vector<std::string> acid_libs) :
            bin_count_{bc},
            bin_size_{bs},
            hash_count_{hc},
            ibf_{seqan3::bin_count{bin_count_},
                 seqan3::bin_size{bin_size_},
                 seqan3::hash_function_count{hash_count_}},
            k_{k},
            molecule_{molecule},
            acid_libs_{acid_libs}
    {
        //static_assert(data_layout_mode == seqan3::data_layout::uncompressed);
    }

    auto getBinCount() const
    {
        assert(ibf_.bin_count() == bin_count_);
        return bin_count_;
    }

    size_t getBinSize() const
    {
        assert(ibf_.bin_size() == bin_size_);
        return bin_size_;
    }

    uint8_t getHashCount() const
    {
        assert(ibf_.hash_function_count() == hash_count_);
        return hash_count_;
    }

    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> getIBF()
    {
        return ibf_;
    }

    void emplace(uint64_t val, size_t idx)
    {
        ibf_.emplace(val, seqan3::bin_index{idx});
    }

    void set_lib_paths(std::vector<std::string> path_collection)
    {
        for(auto && path: path_collection)
        {
            acid_libs_.push_back(path);
        }
    }

    template<class Archive>
    void serialize(Archive & archive)
    {
        archive(bin_count_, bin_size_, hash_count_, ibf_, k_, molecule_, acid_libs_);
    }
};

void read_input_file_list(std::vector<std::filesystem::path> & sequence_files, std::filesystem::path input_file);

template <class IndexStructure>
void store_ibf(IndexStructure const & ibf, std::filesystem::path opath)
{
    std::ofstream os{opath, std::ios::binary};
    cereal::BinaryOutputArchive oarchive{os};
    oarchive(ibf);
}

template <class IndexStructure>
void load_ibf(IndexStructure & ibf, std::filesystem::path ipath)
{
    std::ifstream is{ipath, std::ios::binary};
    cereal::BinaryInputArchive iarchive{is};
    iarchive(ibf);
}

template <typename MolType>
uint32_t parse_reference_na(std::filesystem::path &ref_file, record_list<MolType> &refs)
{
    uint32_t bin_count = 0;
    seqan3::sequence_file_input fin{ref_file};
    record_pair<MolType> record_desc;
    for(auto & record: fin) {
        bin_count++;
        record_desc = std::make_pair(record.id(), record.sequence());
        refs.push_back(record_desc);
    }
    return bin_count;
}

template <typename MolType>
uint32_t parse_reference_aa(std::filesystem::path &ref_file, record_list<MolType> &refs)
{
    uint32_t bin_count = 0;
    using traits_type = seqan3::sequence_file_input_default_traits_aa;
    seqan3::sequence_file_input<traits_type> fin{ref_file};
    record_pair<MolType> record_desc;
    for(auto & record: fin) {
        bin_count++;
        record_desc = std::make_pair(record.id(), record.sequence());
        refs.push_back(record_desc);
    }
    return bin_count;
}

template <typename MolType>
void populate_bin(IndexStructure &ibf, auto &hash_adaptor, record_list<MolType> &records, const size_t &bin_idx)
{
    size_t record_count = records.size();
    for (size_t rec_idx = 0; rec_idx < record_count; rec_idx++) // Generate all the kmers from that file's records
    {
        for (auto && value : records[rec_idx].second | hash_adaptor)
        {
            ibf.emplace(value, bin_idx); // Put all the kmers for one file in one bin
        }
    }
}

void create_index_from_filelist(const index_arguments &cmd_args, const std::vector<std::string> &input_bin_files);

void drive_index(const index_arguments &cmd_args);
