#pragma once


enum {
    Base = 0u,
    Murphy = 1u,
    Li = 2u
};

namespace molecules
{
    class PeptideDecomposer
    {
        private:
            uint8_t ksize_{};
            uint8_t reduction_{};
            uint8_t alphabet_size_;
            uint64_t selection_mask_{};

        public:
            std::array<uint8_t, 256> aamap_{};

            PeptideDecomposer() = default;
            PeptideDecomposer(uint8_t const k, uint8_t const reduction) : ksize_{k}, reduction_{reduction}
            {
                create_residue_maps();
                create_selection_bitmask();
            }

            // void compute_access_masks()
            // {
            //     access_masks_ = new uint64_t[ksize_];
            //     uint64_t mask = 0b11111;
            //     access_masks_[0] = mask;
            //     for(size_t i = 1; i < ksize_; ++i)
            //     {
            //         mask = (mask << 5);
            //         access_masks_[i] = mask;
            //     }
            // }

            void create_selection_bitmask()
            {
                /*
                Example with k=4 creates a bitmask like 
                0b00000000-00000000-00000000-00000000-00000000-00000000-00000000-11111111
                for the rolling hash
                */
                selection_mask_ = (ksize_ >= 32) ? static_cast<uint64_t>(-1) : (1ULL << (5*(ksize_))) - 1u;
            }

            void print_mask() const
            {
                std::cout << selection_mask_ << std::endl;
            }

            void create_residue_maps()
            {
                switch(reduction_)
                {
                    case Murphy:
                        alphabet_size_ = 10;
                        aamap_['A'] = 0;
                        aamap_['B'] = 1;
                        aamap_['C'] = 2;
                        aamap_['F'] = 3;
                        aamap_['G'] = 4;
                        aamap_['H'] = 5;
                        aamap_['I'] = 6;
                        aamap_['K'] = 7;
                        aamap_['P'] = 8;
                        aamap_['S'] = 9;
                        aamap_['D'] = aamap_['B'];
                        aamap_['E'] = aamap_['B'];
                        aamap_['L'] = aamap_['I'];
                        aamap_['M'] = aamap_['I'];
                        aamap_['N'] = aamap_['B'];
                        aamap_['Q'] = aamap_['B'];
                        aamap_['R'] = aamap_['K'];
                        aamap_['T'] = aamap_['S'];
                        aamap_['V'] = aamap_['I'];
                        aamap_['W'] = aamap_['F'];
                        aamap_['Y'] = aamap_['F'];
                        aamap_['J'] = aamap_['I'];
                        aamap_['O'] = aamap_['K'];
                        aamap_['U'] = aamap_['C'];
                        aamap_['X'] = aamap_['S'];
                        aamap_['Z'] = aamap_['B'];
                    case Li:
                        alphabet_size_ = 10;
                        aamap_['A'] = 0;
                        aamap_['B'] = 1;
                        aamap_['C'] = 2;
                        aamap_['F'] = 3;
                        aamap_['G'] = 4;
                        aamap_['H'] = 5;
                        aamap_['I'] = 6;
                        aamap_['J'] = 7;
                        aamap_['K'] = 8;
                        aamap_['P'] = 9;
                        aamap_['D'] = aamap_['B'];
                        aamap_['E'] = aamap_['B'];
                        aamap_['L'] = aamap_['J'];
                        aamap_['M'] = aamap_['J'];
                        aamap_['N'] = aamap_['H'];
                        aamap_['Q'] = aamap_['B'];
                        aamap_['R'] = aamap_['K'];
                        aamap_['T'] = aamap_['A'];
                        aamap_['V'] = aamap_['I'];
                        aamap_['W'] = aamap_['F'];
                        aamap_['Y'] = aamap_['F'];
                        aamap_['S'] = aamap_['A'];
                        aamap_['O'] = aamap_['K'];
                        aamap_['U'] = aamap_['C'];
                        aamap_['X'] = aamap_['A'];
                        aamap_['Z'] = aamap_['B'];
                    default: // Base
                        alphabet_size_ = 20;
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
            }

            // uint64_t compute_hash(const uint64_t &kmer)
            // {
            //     uint64_t hash;
            //     size_t res1, res2, res3, res4;
            //     size_t numbElements = ksize_;
            //     switch(numbElements)
            //     {
            //         case 6:
            //             res1 = (kmer&access_masks_[0]);
            //             res2 = (kmer&access_masks_[1]);
            //             res3 = (kmer&access_masks_[2]);
            //             res4 = (kmer&access_masks_[3]);
            //             res1 += (kmer&access_masks_[4]);
            //             res2 += (kmer&access_masks_[5]);
            //             hash = res1 + res2 + res3 + res4;
            //             break;
            //         case 7:
            //             res1 = (kmer&access_masks_[0]);
            //             res2 = (kmer&access_masks_[1]);
            //             res3 = (kmer&access_masks_[2]);
            //             res4 = (kmer&access_masks_[3]);
            //             res1 += (kmer&access_masks_[4]);
            //             res2 += (kmer&access_masks_[5]);
            //             res3 += (kmer&access_masks_[6]);
            //             hash = res1 + res2 + res3 + res4;
            //             break;
            //         case 10:
            //             res1 = (kmer&access_masks_[0]);
            //             res2 = (kmer&access_masks_[1]);
            //             res3 = (kmer&access_masks_[2]);
            //             res4 = (kmer&access_masks_[3]);
            //             res1 += (kmer&access_masks_[4]);
            //             res2 += (kmer&access_masks_[5]);
            //             res3 += (kmer&access_masks_[6]);
            //             res4 += (kmer&access_masks_[7]);
            //             res1 += (kmer&access_masks_[8]);
            //             res2 += (kmer&access_masks_[9]);
            //             hash = res1 + res2 + res3 + res4;
            //             break;
            //         default:
            //             hash = kmer;
            //             break;
            //     }
            //     return hash;
            // }

            uint64_t encode_peptide(std::string_view const kmer) const
            {
                uint64_t codemer{};
                for(auto && base: kmer)
                {
                    codemer <<= 5;
                    codemer |= aamap_[base];
                }
                return codemer;
            }

            template<typename T>
            void rollover_peptide_hash(const unsigned char residue, T const bin_id, auto &base_ref) const
            {
                auto encoded_residue = aamap_[residue];
                base_ref.forward_store_ = ((base_ref.forward_store_<<5)&selection_mask_) | encoded_residue;
                base_ref.emplace(base_ref.forward_store_, bin_id);
            }

            template<typename T>
            void decompose_record(std::string_view const record_seq, T const tech_bin_id, auto &base_ref)
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

            uint64_t update_kmer(const int &symbol, uint64_t &kmer)
            {
                uint64_t fb = aamap_[symbol];
                uint64_t forward = ((kmer<<5)&selection_mask_) | fb;
                kmer = forward;
                return forward;
            }

            template<class Archive>
            void serialize(Archive &archive)
            {
                archive(ksize_, reduction_, alphabet_size_, selection_mask_, aamap_);
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
