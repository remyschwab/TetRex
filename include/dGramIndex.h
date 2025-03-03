#pragma once


#include "utils.h"
#include "kseq.h"



class DGramIndex
{
    private:
        size_t l_{};
        size_t u_{};
        size_t pad_{};
        uint8_t hash_count_{};
        float fpr_{};
        std::vector<std::filesystem::path> user_bins_{};
        seqan::hibf::hierarchical_interleaved_bloom_filter dibf_{};
        uint32_t dist_{};
        uint32_t span_{};
        

    public:
        std::optional<seqan::hibf::hierarchical_interleaved_bloom_filter::membership_agent_type> agent_{};
        std::array<uint8_t, 256> aamap_{};
        std::array<uint8_t, 128> dmap_{}; 


    DGramIndex() = default;


    explicit DGramIndex(size_t lower, size_t upper, size_t pad, uint8_t hc, float fpr, std::vector<std::filesystem::path> user_bins) :
            l_{lower},
            u_{upper},
            pad_{pad},
            hash_count_{hc},
            fpr_{fpr},
            user_bins_{user_bins},
            dibf_{}
            {
                create_distance_maps();
                create_residue_maps();
            }

    void track_record(std::vector<std::vector<uint64_t>> &dgrams, const std::string &ref_)
    {
        std::vector<uint64_t> inner_list;
        for(size_t i = 0; i < ref_.length()-span_+1; ++i)
        {
            uint64_t dmer = 0;
            auto nterm_residue = aamap_[ref_.at(i)];
            auto cterm_residue = aamap_[ref_.at(i+span_-1)];
            dmer |= nterm_residue;
            dmer = (dmer<<5);
            dmer |= cterm_residue;
            dmer = dmer<<32;
            dmer |= dist_;
            inner_list.push_back(dmer);
        }
        dgrams.push_back(inner_list);
    }


    void track_bins()
    {
        std::vector<std::vector<std::vector<uint64_t>>> multigrams;
        gzFile handle{};
        kseq_t *record{};
        int status{};
        size_t seq_count{};
        for(size_t i = 0; i < user_bins_.size(); ++i) // Iterate over bins
        {
            for(dist_ = l_; dist_ <= u_; ++dist_)
            {
                seqan3::debug_stream << user_bins_[i].c_str() << " : " << dist_ << std::endl;
                span_ = dist_ + pad_ + pad_;
                std::vector<std::vector<uint64_t>> dgrams;
                handle = gzopen(user_bins_[i].c_str(), "r");
                record = kseq_init(handle);
                while ((status = kseq_read(record)) >= 0) // Iterate over bin records
                {
                    if(record->seq.l < span_) continue;
                    track_record(dgrams, record->seq.s);
                }
                kseq_destroy(record);
                gzclose(handle);
            }
        }
    }


    size_t compute_bitcount(const size_t bfn) const
    {
        size_t bitcount;
        double numerator = -static_cast<double>(bfn) * std::log(fpr_);
        double denominator = std::pow(std::log(2), 2);
        bitcount = static_cast<size_t>(std::ceil(numerator / denominator));
        return bitcount;
    }


    void create_residue_maps()
    {
        aamap_['A'] = 0u;
        aamap_['C'] = 1u;
        aamap_['D'] = 2u;
        aamap_['E'] = 3u;
        aamap_['F'] = 4u;
        aamap_['G'] = 5u;
        aamap_['H'] = 6u;
        aamap_['I'] = 7u;
        aamap_['K'] = 8u;
        aamap_['L'] = 9u;
        aamap_['M'] = 10u;
        aamap_['N'] = 11u;
        aamap_['P'] = 12u;
        aamap_['Q'] = 13u;
        aamap_['R'] = 14u;
        aamap_['S'] = 15u;
        aamap_['T'] = 16u;
        aamap_['V'] = 17u;
        aamap_['W'] = 18u;
        aamap_['Y'] = 19u;
        aamap_['X'] = 20u;
        aamap_['B'] = aamap_['D'];
        aamap_['J'] = aamap_['L'];
        aamap_['O'] = aamap_['X'];
        aamap_['U'] = aamap_['X'];
        aamap_['Z'] = aamap_['E'];
    }

    void create_distance_maps()
    {
        // No reduction
        for(size_t i = 0; i < 128; ++i) dmap_[i] = i;
    }

    // seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> getDIBF()
    // {
    //     return dibf_;
    // }


    template<class Archive>
    void serialize(Archive &archive)
    {
        archive(l_, u_, pad_, hash_count_, fpr_);
    }
};
