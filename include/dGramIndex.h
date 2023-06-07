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
        std::vector<size_t> tracker_;


    DGramIndex() = default;


    explicit DGramIndex(size_t bc, size_t bs, uint8_t hc) :
        bin_count_{bc},
        bin_size_{bs},
        hash_count_{hc},
        dibf_{seqan3::bin_count{bin_count_}, seqan3::bin_size{bin_size_}, seqan3::hash_function_count{hash_count_}}
    {
        make_alphabet();
        tracker_.resize(UCHAR_MAX+1, UCHAR_MAX);
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
        tracker_['A'] = -1;
        tracker_['C'] = -1;
        tracker_['D'] = -1;
        tracker_['E'] = -1;
        tracker_['F'] = -1;
        tracker_['G'] = -1;
        tracker_['H'] = -1;
        tracker_['I'] = -1;
        tracker_['K'] = -1;
        tracker_['L'] = -1;
        tracker_['M'] = -1;
        tracker_['N'] = -1;
        tracker_['P'] = -1;
        tracker_['Q'] = -1;
        tracker_['R'] = -1;
        tracker_['S'] = -1;
        tracker_['T'] = -1;
        tracker_['V'] = -1;
        tracker_['W'] = -1;
        tracker_['Y'] = -1;
    }


    void emplace(uint64_t val, size_t idx)
    {
        dibf_.emplace(val, seqan3::bin_index{idx});
    }


    void track_record(const size_t &i, const unsigned char &residue, const size_t &bin_idx)
    {
        uint64_t value_encoding = 0;
        uint32_t distance;
        for(unsigned char &alpha: alphabet_)
        {
            if(tracker_[alpha] > 0)
            {
                value_encoding = (value_encoding | residue)<<8;
                value_encoding = (value_encoding | alpha)<<8;
                distance = i-tracker_[alpha];
                value_encoding = (value_encoding<<32) | distance;
                emplace(value_encoding, bin_idx);
            }
            tracker_[residue] = i;
        }
    }


    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> getDIBF()
    {
        return dibf_;
    }
};
