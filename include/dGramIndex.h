#pragma once


#include "utils.h"
#include "arg_parse.h"
#include "kseq.h"
#include "molecule_encodings.h"
#include "robin_hood.h"

#include <zlib.h>

#include <cereal/types/vector.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/views/kmer_hash.hpp>
#include <seqan3/alphabet/nucleotide/all.hpp>


class DGramIndex
{
    private:
        size_t bin_count_;
        size_t bin_size_{};
        uint8_t hash_count_{};
        seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> dibf_;

    public:
        std::vector<unsigned char> alphabet_;
        std::vector<size_t> aamap_;


    DGramIndex() = default;


    explicit DGramIndex(size_t bc, size_t bs, size_t hc) :
        bin_count_{bc},
        bin_size_{bs},
        hash_count_{hc},
        dibf_{seqan3::bin_count{bin_count_}, seqan3::bin_size{bin_size_}, seqan3::hash_function_count{hash_count_}}

    {
        make_alphabet();
        aamap_.resize(UCHAR_MAX+1, UCHAR_MAX);
        // init_tracker();
    }

    void make_alphabet()
    {
        std::string alpha_str = "acdefghiklmnpqrstvwy";
        for(auto &residue: alpha_str)
            alphabet_.push_back(toupper(residue));
    }

    void init_tracker()
    {
        aamap_['A'] = -1;
        aamap_['C'] = -1;
        aamap_['D'] = -1;
        aamap_['E'] = -1;
        aamap_['F'] = -1;
        aamap_['G'] = -1;
        aamap_['H'] = -1;
        aamap_['I'] = -1;
        aamap_['K'] = -1;
        aamap_['L'] = -1;
        aamap_['M'] = -1;
        aamap_['N'] = -1;
        aamap_['P'] = -1;
        aamap_['Q'] = -1;
        aamap_['R'] = -1;
        aamap_['S'] = -1;
        aamap_['T'] = -1;
        aamap_['V'] = -1;
        aamap_['W'] = -1;
        aamap_['Y'] = -1;
    }

    void emplace(uint64_t val, size_t idx)
    {
        dibf_.emplace(val, seqan3::bin_index{idx});
    }

    void track_record(const size_t &i, const unsigned char &residue, const size_t &bin_idx)
    {
        uint64_t value_encoding;
        uint32_t distance;
        for(auto &alpha: alphabet_)
        {
            if(aamap_[alpha] > 0)
            {
                value_encoding = (value_encoding & residue)<<8;
                value_encoding = (value_encoding & alpha)<<8;
                distance = i-aamap_[alpha];
                value_encoding = (value_encoding <<32) & distance;
                seqan3::debug_stream << value_encoding << "," << bin_idx << std::endl;
                emplace(value_encoding, bin_idx);
            }
            aamap_[residue] = i;
        }
    }

    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> getDIBF()
    {
        return dibf_;
    }

    uint64_t inline MurmurHash3(size_t key)
    {
        uint64_t k = (uint64_t) key;
        k ^= k >> 33;
        k *= 0xff51afd7ed558ccd;
        k ^= k >> 33;
        k *= 0xc4ceb9fe1a85ec53;
        k ^= k >> 33;
        return k;
    }
};
