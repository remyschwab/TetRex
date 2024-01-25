//
// Created by Remy Schwab on 13.09.22.
//

#pragma once

#include "index_includes.h"

KSEQ_INIT(gzFile, gzread)


template<molecules::is_molecule molecule_t>
class IBFIndex
{

private:
    size_t bin_count_;
    size_t bin_size_;
    uint8_t hash_count_;
    std::vector<std::string> tech_bins_;
    index_structure::IBF ibf_;


public:
    index_structure::IBF_Agent agent_{};

    IBFIndex() = default;

    explicit IBFIndex(size_t bin_size, uint8_t hc, std::vector<std::string> tech_bins) :
            bin_count_{tech_bins.size()},
            bin_size_{bin_size},
            hash_count_{hc},
            tech_bins_{tech_bins},
            ibf_(seqan::hibf::bin_count {bin_count_},
             seqan::hibf::bin_size {bin_size_},
              seqan::hibf::hash_function_count{hash_count_})
    {
        
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

    index_structure::IBF& getIBF()
    {
        return ibf_;
    }

    void spawn_agent()
    {
        agent_ = ibf_.membership_agent();
    }

    void emplace(uint64_t &val, size_t &idx)
    {
        ibf_.emplace(val, seqan::hibf::bin_index{idx});
    }

    template<molecules::is_molecule molecule_t>
    void populate_index(uint8_t &k, MoleculeDecomposer<molecule_t> &decomposer)
    {
        gzFile handle;
        kseq_t *record;
        int status;
        size_t seq_count = 0;
        for(size_t i = 0; i < tech_bins_.size(); ++i) // Iterate over bins
        {
            handle = gzopen(tech_bins_[i].c_str(), "r");
            record = kseq_init(handle);
            while ((status = kseq_read(record)) >= 0) // Iterate over bin records
            {
                std::string_view record_view = record->seq.s;
                if(record_view.length() < k)
                {
                    seqan3::debug_stream << "RECORD TOO SHORT " << record->comment.s << std::endl;
                    continue;
                }
                seq_count++;

                decomposer.decompose_record(record_view, ibf_, i);

                // std::vector<uint8_t> peptide_nums(record_view.length());
                // for(size_t o = 0; o < record_view.length(); ++o)
                // {
                //     peptide_nums[o] = ibf.aamap_[record_view[o]];
                // }
                // for(size_t o = 0; o < peptide_nums.size()-ibf.k_+1; ++o)
                // {
                //     decompose_peptide_record(peptide_nums, o, ibf);
                //     ibf.emplace(ibf.forward_store_, i);
                // }
            }
            kseq_destroy(record);
            gzclose(handle);
        }
        seqan3::debug_stream << "Indexed " << seq_count << " sequences across " << bin_count_ << " bins." << std::endl;
    }

    const bitvector & query(uint64_t &kmer)
    {
        return agent_.bulk_contains(kmer);
    }

    template<class Archive>
    void serialize(Archive &archive)
    {
        archive(bin_count_, bin_size_, hash_count_, ibf_);
    }
};
