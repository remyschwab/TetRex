//
// Created by Remy Schwab on 13.09.22.
//

#pragma once


#include "utils.h"
#include "arg_parse.h"

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
    uint32_t bin_count_;
    uint64_t bin_size_{};
    uint8_t hash_count_{};
    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf_;

public:

    uint8_t k_;
    std::string molecule_;

    IndexStructure() = default;

    explicit IndexStructure(uint8_t k,
                            uint32_t bc,
                            uint64_t bs,
                            uint8_t hc,
                            std::string molecule) :
            bin_count_{bc},
            bin_size_{bs},
            hash_count_{hc},
            ibf_{seqan3::bin_count{bin_count_},
                 seqan3::bin_size{bin_size_},
                 seqan3::hash_function_count{hash_count_}},
            k_{k},
            molecule_{molecule}
    {
        //static_assert(data_layout_mode == seqan3::data_layout::uncompressed);
    }

    auto getBinCount() const
    {
        assert(ibf_.bin_count() == bin_count_);
        return bin_count_;
    }

    uint64_t getBinSize() const
    {
        assert(ibf_.bin_size() == bin_size_);
        return bin_size_;
    }

    uint8_t getHashCount() const
    {
        assert(ibf_.hash_function_count() == hash_count_);
        return hash_count_;
    }

    // auto makeAgent()
    // {
    //     auto agent = ibf_.membership_agent();
    //     return agent;
    // }

    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> getIBF()
    {
        return ibf_;
    }

    void emplace(uint64_t val, uint32_t idx)
    {
        ibf_.emplace(val, seqan3::bin_index{idx});
    }

    template<class Archive>
    void serialize(Archive & archive)
    {
        archive(bin_count_, bin_size_, hash_count_, ibf_, k_, molecule_);
    }
};

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
IndexStructure create_index(record_list<MolType> &refs, uint32_t &bin_count, index_arguments args)
{
    uint8_t k = args.k;
    IndexStructure ibf(k, bin_count, args.bin_size, args.hash_count, args.molecule);

    auto hash_adaptor = seqan3::views::kmer_hash(seqan3::ungapped{k});
    for (size_t i = 0; i < bin_count; i++)
        for (auto && value : refs[i].second | hash_adaptor)
        {
            ibf.emplace(value, i);
        }
    return ibf;
}

void drive_index(const index_arguments &cmd_args);
