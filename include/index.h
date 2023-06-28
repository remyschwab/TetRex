//
// Created by Remy Schwab on 13.09.22.
//

#pragma once


#include "utils.h"
#include "arg_parse.h"
#include "kseq.h"
#include "molecule_encodings.h"
#include "robin_hood.h"
#include "dGramIndex.h"

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
    std::string reduction_;
    uint8_t alphabet_size_; // Amino Acid Encoding
    size_t * powers_;
    std::vector<uint8_t> aamap_;
    // DGramIndex dibf_;
    // Not serialized
    uint64_t selection_mask_; // DNA Encoding
    uint8_t left_shift_;
    uint64_t forward_store_;
    uint64_t reverse_store_;


    IndexStructure() = default;

    explicit IndexStructure(uint8_t k,
                            size_t bc,
                            size_t bs,
                            uint8_t hc,
                            std::string molecule,
                            std::vector<std::string> acid_libs,
                            std::string reduction) :
            bin_count_{bc},
            bin_size_{bs},
            hash_count_{hc},
            ibf_{seqan3::bin_count{bin_count_},
                 seqan3::bin_size{bin_size_},
                 seqan3::hash_function_count{hash_count_}},
            k_{k},
            molecule_{molecule},
            acid_libs_{acid_libs},
            reduction_{reduction}
    {
        //static_assert(data_layout_mode == seqan3::data_layout::uncompressed);
        if(molecule_ == "na")
        {
            create_selection_bitmask();
            set_left_shift();
        }
        else
        {
            set_alphabet_maps();
            compute_powers();
        }
    }

    void create_selection_bitmask()
    {
        /*
        Example with k=4 creates a bitmask like 
        0b00000000-00000000-00000000-00000000-00000000-00000000-00000000-11111111
        for the rolling hash
        */
        size_t countdown = k_;
        while(countdown > 0)
        {
            selection_mask_ = (selection_mask_<<2) | 0b11;
            countdown--;
        }
    }

    void set_stores(uint64_t forward, uint64_t reverse)
    {
        forward_store_ = forward;
        reverse_store_ = reverse;
    }

    std::tuple<uint64_t, uint64_t> get_stores()
    {
        return {forward_store_, reverse_store_};
    }

    void set_left_shift()
    {
        left_shift_ = ((uint8_t)2 * k_)-2;
    }

    void compute_powers()
    {
        powers_ = new size_t[k_];
        size_t pow = 1;
        for(size_t i = 0; i < k_; ++i)
        {
            powers_[i] = pow;
            pow *= alphabet_size_;
        }
    }

    void set_alphabet_maps()
    {
        uint8_t fall;
        if(reduction_ == "None")
        {
            alphabet_size_ = 20;
            fall = 0;
            create_residue_maps(fall, aamap_);
        }
        else if(reduction_ == "murphy")
        {
            alphabet_size_ = 10;
            fall = 1;
            create_residue_maps(fall, aamap_);
        }
        else if(reduction_ == "li")
        {
            alphabet_size_ = 10;
            fall = 2;
            create_residue_maps(fall, aamap_);
        }
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

    void rollover_nuc_hash(const char base , const size_t &bin_id)
    {
        auto fb = (base>>1)&3; // Encode the new base
        auto cb = (fb^0b10)<<left_shift_; // Get its complement and shift it to the big end of the uint64
        forward_store_ = ((forward_store_<<2)&selection_mask_) | fb; // Update forward store 
        reverse_store_ = ((reverse_store_>>2)&selection_mask_) | cb; // Update reverse store
        emplace(( forward_store_ <= reverse_store_ ? forward_store_ : reverse_store_ ), bin_id);
    }

    void emplace(uint64_t val, size_t idx)
    {
        ibf_.emplace(val, seqan3::bin_index{idx});
    }

    void set_lib_paths(std::vector<std::string> path_collection)
    {
        acid_libs_ = path_collection;
    }

    // void set_dibf(DGramIndex &dibf)
    // {
    //     dibf_ = dibf;
    // }

    template<class Archive>
    void serialize(Archive &archive)
    {
        archive(bin_count_, bin_size_, hash_count_, ibf_, k_, molecule_, acid_libs_, reduction_, alphabet_size_);
    }

    uint8_t map_aa(unsigned char &residue)
    {
        return aamap_[residue];
    }
};

void convertStringToVec(const std::string &kmer, const IndexStructure &ibf, std::vector<unsigned char> &digit_vec);

void read_input_file_list(std::vector<std::filesystem::path> & sequence_files, std::filesystem::path input_file);

void decompose_nucleotide_record(std::string_view record_seq, IndexStructure &ibf, const size_t, const size_t &bin_id);

void decompose_peptide_record(const std::vector<unsigned char> &record_nums, const size_t &begin, IndexStructure &ibf);

void create_index(const index_arguments &cmd_args, const std::vector<std::string> &input_bin_files);

void drive_index(const index_arguments &cmd_args);

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
