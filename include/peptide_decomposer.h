#pragma once

#include <cereal/types/array.hpp>


enum {
    Base = 0,
    Murphy = 1,
    Li = 2
};

namespace molecules
{
    class PeptideDecomposer
    {
        private:
            uint8_t ksize_{};
            uint8_t reduction_{};
            uint8_t alphabet_size_;
            uint64_t *access_masks_;
            uint64_t selection_mask_;

        public:
            std::vector<uint8_t> aamap_;

            PeptideDecomposer() = default;
            PeptideDecomposer(uint8_t &k, uint8_t &reduction) : ksize_{k}, reduction_{reduction}
            {
                create_residue_maps(reduction_, aamap_);
                create_selection_bitmask();
                // compute_access_masks();
            }

            void compute_access_masks()
            {
                access_masks_ = new uint64_t[ksize_];
                uint64_t mask = 0b11111;
                access_masks_[0] = mask;
                for(size_t i = 1; i < ksize_; ++i)
                {
                    mask = (mask << 5);
                    access_masks_[i] = mask;
                }
            }

            void create_selection_bitmask()
            {
                /*
                Example with k=4 creates a bitmask like 
                0b00000000-00000000-00000000-00000000-00000000-00000000-00000000-11111111
                for the rolling hash
                */
                size_t countdown = ksize_;
                while(countdown > 0)
                {
                    selection_mask_ = (selection_mask_<<5) | 0b11111;
                    countdown--;
                }
            }

            void create_residue_maps(uint8_t &alphabet, std::vector<uint8_t> &aamap)
            {
                aamap.resize(UCHAR_MAX+1, UCHAR_MAX);
                switch(alphabet)
                {
                    case Murphy:
                        alphabet_size_ = 10;
                        aamap['A'] = 0;
                        aamap['B'] = 1;
                        aamap['C'] = 2;
                        aamap['F'] = 3;
                        aamap['G'] = 4;
                        aamap['H'] = 5;
                        aamap['I'] = 6;
                        aamap['K'] = 7;
                        aamap['P'] = 8;
                        aamap['S'] = 9;
                        aamap['D'] = aamap['B'];
                        aamap['E'] = aamap['B'];
                        aamap['L'] = aamap['I'];
                        aamap['M'] = aamap['I'];
                        aamap['N'] = aamap['B'];
                        aamap['Q'] = aamap['B'];
                        aamap['R'] = aamap['K'];
                        aamap['T'] = aamap['S'];
                        aamap['V'] = aamap['I'];
                        aamap['W'] = aamap['F'];
                        aamap['Y'] = aamap['F'];
                        aamap['J'] = aamap['I'];
                        aamap['O'] = aamap['K'];
                        aamap['U'] = aamap['C'];
                        aamap['X'] = aamap['S'];
                        aamap['Z'] = aamap['B'];
                    case Li:
                        alphabet_size_ = 10;
                        aamap['A'] = 0;
                        aamap['B'] = 1;
                        aamap['C'] = 2;
                        aamap['F'] = 3;
                        aamap['G'] = 4;
                        aamap['H'] = 5;
                        aamap['I'] = 6;
                        aamap['J'] = 7;
                        aamap['K'] = 8;
                        aamap['P'] = 9;
                        aamap['D'] = aamap['B'];
                        aamap['E'] = aamap['B'];
                        aamap['L'] = aamap['J'];
                        aamap['M'] = aamap['J'];
                        aamap['N'] = aamap['H'];
                        aamap['Q'] = aamap['B'];
                        aamap['R'] = aamap['K'];
                        aamap['T'] = aamap['A'];
                        aamap['V'] = aamap['I'];
                        aamap['W'] = aamap['F'];
                        aamap['Y'] = aamap['F'];
                        aamap['S'] = aamap['A'];
                        aamap['O'] = aamap['K'];
                        aamap['U'] = aamap['C'];
                        aamap['X'] = aamap['A'];
                        aamap['Z'] = aamap['B'];
                    default: // Base
                        alphabet_size_ = 20;
                        aamap['A'] = 0;
                        aamap['C'] = 1;
                        aamap['D'] = 2;
                        aamap['E'] = 3;
                        aamap['F'] = 4;
                        aamap['G'] = 5;
                        aamap['H'] = 6;
                        aamap['I'] = 7;
                        aamap['K'] = 8;
                        aamap['L'] = 9;
                        aamap['M'] = 10;
                        aamap['N'] = 11;
                        aamap['P'] = 12;
                        aamap['Q'] = 13;
                        aamap['R'] = 14;
                        aamap['S'] = 15;
                        aamap['T'] = 16;
                        aamap['V'] = 17;
                        aamap['W'] = 18;
                        aamap['Y'] = 19;
                        aamap['X'] = 20;
                        aamap['B'] = aamap['D'];
                        aamap['J'] = aamap['L'];
                        aamap['O'] = aamap['X'];
                        aamap['U'] = aamap['X'];
                        aamap['Z'] = aamap['E'];
                }
            }

            uint64_t compute_hash(const uint64_t &kmer)
            {
                uint64_t hash;
                size_t res1, res2, res3, res4;
                size_t numbElements = ksize_;
                switch(numbElements)
                {
                    case 6:
                        res1 = (kmer&access_masks_[0]);
                        res2 = (kmer&access_masks_[1]);
                        res3 = (kmer&access_masks_[2]);
                        res4 = (kmer&access_masks_[3]);
                        res1 += (kmer&access_masks_[4]);
                        res2 += (kmer&access_masks_[5]);
                        hash = res1 + res2 + res3 + res4;
                        break;
                    case 7:
                        res1 = (kmer&access_masks_[0]);
                        res2 = (kmer&access_masks_[1]);
                        res3 = (kmer&access_masks_[2]);
                        res4 = (kmer&access_masks_[3]);
                        res1 += (kmer&access_masks_[4]);
                        res2 += (kmer&access_masks_[5]);
                        res3 += (kmer&access_masks_[6]);
                        hash = res1 + res2 + res3 + res4;
                        break;
                    case 10:
                        res1 = (kmer&access_masks_[0]);
                        res2 = (kmer&access_masks_[1]);
                        res3 = (kmer&access_masks_[2]);
                        res4 = (kmer&access_masks_[3]);
                        res1 += (kmer&access_masks_[4]);
                        res2 += (kmer&access_masks_[5]);
                        res3 += (kmer&access_masks_[6]);
                        res4 += (kmer&access_masks_[7]);
                        res1 += (kmer&access_masks_[8]);
                        res2 += (kmer&access_masks_[9]);
                        hash = res1 + res2 + res3 + res4;
                        break;
                    default:
                        hash = kmer;
                        break;
                }
                return hash;
            }

            uint64_t encode_peptide(std::string_view kmer)
            {
                uint64_t codemer = 0;
                for(auto && base: kmer)
                {
                    codemer = codemer << 5;
                    codemer += aamap_[base];
                }
                return codemer;
            }

            void rollover_peptide_hash(const char residue, const seqan::hibf::bin_index &bin_id, auto &base_ref)
            {
                auto encoded_residue = aamap_[residue];
                base_ref.forward_store_ = ((base_ref.forward_store_<<5)&selection_mask_) | encoded_residue;
                base_ref.emplace(base_ref.forward_store_, bin_id);
            }

            void decompose_record(std::string_view record_seq, seqan::hibf::bin_index &tech_bin_id, auto &base_ref)
            {
                uint64_t initial_encoding = encode_peptide(record_seq.substr(0, ksize_));
                base_ref.forward_store_ = initial_encoding;
                base_ref.emplace(initial_encoding, tech_bin_id);
                for(size_t i = ksize_; i < record_seq.length(); ++i)
                {
                    auto symbol = record_seq[i];
                    rollover_peptide_hash(symbol, tech_bin_id, base_ref);
                }
            }

            void update_kmer(const int &symbol, uint64_t &kmer)
            {
                uint64_t fb = aamap_[symbol];
                uint64_t forward = ((kmer<<5)&selection_mask_) | fb;
                kmer = forward;
            }

            template<class Archive>
            void serialize(Archive &archive)
            {
                archive(ksize_, reduction_, alphabet_size_, access_masks_, aamap_);
            }
    };

    // void create_molecule_bitmaps()
    // {
    //     aamap['A'] = 0b00000; // 0
    //     aamap['C'] = 0b00001; // 1
    //     aamap['D'] = 0b00010; // 2
    //     aamap['E'] = 0b00011; // 3
    //     aamap['F'] = 0b00100; // 4
    //     aamap['G'] = 0b00101; // 5
    //     aamap['H'] = 0b00110; // 6
    //     aamap['I'] = 0b00111; // 7
    //     aamap['K'] = 0b01000; // 8
    //     aamap['L'] = 0b01001; // 9
    //     aamap['M'] = 0b01010; // 10
    //     aamap['N'] = 0b01011; // 11
    //     aamap['P'] = 0b01100; // 12
    //     aamap['Q'] = 0b01101; // 13
    //     aamap['R'] = 0b01110; // 14
    //     aamap['S'] = 0b01111; // 15
    //     aamap['T'] = 0b10000; // 16
    //     aamap['V'] = 0b10001; // 17
    //     aamap['W'] = 0b10010; // 18
    //     aamap['Y'] = 0b10011; // 19
    //     aamap['B'] = aamap['D'];
    //     aamap['J'] = aamap['L'];
    //     aamap['O'] = aamap['X'];
    //     aamap['U'] = aamap['X'];
    //     aamap['Z'] = aamap['E'];
    // }

} // end namespace molecules
