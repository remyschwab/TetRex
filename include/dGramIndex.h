#pragma once


class DGramIndex
{
    private:
        size_t bin_count_;
        size_t bin_size_{};
        uint8_t hash_count_{};
        seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> dibf_;

    public:
        std::vector<unsigned char> alphabet_;
        std::vector<int> tracker_;
        std::vector<uint8_t> bin_cache_;


    DGramIndex() = default;


    explicit DGramIndex(size_t bc, size_t bs, uint8_t hc) :
        bin_count_{bc},
        bin_size_{bs},
        hash_count_{hc},
        dibf_{seqan3::bin_count{bin_count_}, seqan3::bin_size{bin_size_}, seqan3::hash_function_count{hash_count_}}
    {
        make_alphabet();
        tracker_.resize(UCHAR_MAX);
        reset_tracker();
        bin_cache_.resize(UINT32_MAX);
    }


    void make_alphabet()
    {
        std::string alpha_str = "acdefghiklmnpqrstvwy";
        for(auto &residue: alpha_str)
            alphabet_.push_back(toupper(residue));
    }


    void reset_tracker()
    {
        tracker_['A'] = -1;
        tracker_['B'] = -1;
        tracker_['C'] = -1;
        tracker_['D'] = -1;
        tracker_['E'] = -1;
        tracker_['F'] = -1;
        tracker_['G'] = -1;
        tracker_['H'] = -1;
        tracker_['I'] = -1;
        tracker_['J'] = -1;
        tracker_['K'] = -1;
        tracker_['L'] = -1;
        tracker_['M'] = -1;
        tracker_['N'] = -1;
        tracker_['O'] = -1;
        tracker_['P'] = -1;
        tracker_['Q'] = -1;
        tracker_['R'] = -1;
        tracker_['S'] = -1;
        tracker_['T'] = -1;
        tracker_['U'] = -1;
        tracker_['V'] = -1;
        tracker_['W'] = -1;
        tracker_['X'] = -1;
        tracker_['Y'] = -1;
        tracker_['Z'] = -1;
    }


    void init_bin_cache()
    {
        std::fill(bin_cache_.begin(), bin_cache_.end(), 0);
    }


    size_t getBinCount()
    {
        return bin_count_;
    }


    void emplace(uint64_t val, size_t idx)
    {
        dibf_.emplace(val, seqan3::bin_index{idx});
    }


    void track_record(const size_t &record_idx, const unsigned char &residue, const size_t &bin_idx)
    {
        uint32_t dgram_encoding = 0;
        uint16_t distance;
        for(unsigned char &alpha: alphabet_) // Only iterate over alphabetical indices
        {
            if(tracker_[alpha] >= 0) // If this residue/base has been seen before
            {
                dgram_encoding = (dgram_encoding | residue)<<8;
                dgram_encoding = (dgram_encoding | alpha)<<8;
                distance = record_idx-tracker_[alpha];
                dgram_encoding = (dgram_encoding<<16) | distance;
                seqan3::debug_stream << static_cast<char>(residue) << " " << static_cast<char>(alpha) << " " << distance << " " << dgram_encoding << std::endl;
                if(bin_cache_[dgram_encoding] == 0) // Avoid computing hashes for dgrams seen in a bin
                {
                    emplace(dgram_encoding, bin_idx);
                    bin_cache_[dgram_encoding] = 1;
                    // seqan3::debug_stream << static_cast<char>(residue) << " " << static_cast<char>(alpha) << " " << distance << std::endl;
                }
            }
        }
        tracker_[residue] = record_idx;
    }


    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> getDIBF()
    {
        return dibf_;
    }


    template<class Archive>
    void serialize(Archive &archive)
    {
        archive(bin_count_, bin_size_, hash_count_, dibf_);
    }
};
