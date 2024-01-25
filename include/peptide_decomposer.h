#pragma once


#include <vector>


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
            size_t * powers_;

        public:
            std::vector<uint8_t> aamap_;

            PeptideDecomposer() = default;
            PeptideDecomposer(uint8_t &k, uint8_t &reduction) : ksize_{k}, reduction_{reduction}
            {
                create_residue_maps(reduction_, aamap_);
                compute_powers();
            }

            void compute_powers()
            {
                powers_ = new size_t[ksize_];
                size_t pow = 1;
                for(size_t i = 0; i < ksize_; ++i)
                {
                    powers_[i] = pow;
                    pow *= alphabet_size_;
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

            template<index_structure::is_valid ibf_flavor>
            void decompose(const std::vector<unsigned char> &int_seq, const size_t &begin, TetrexIndex<ibf_flavor, molecules::peptide> &ibf)
            {
                size_t res1, res2, res3, res4;
                size_t numbElements = ibf.k_;
                size_t end = begin+ibf.k_;
                switch(numbElements)
                {
                    case 6:
                        res1 = int_seq[begin+0]*ibf.powers_[0];
                        res2 = int_seq[begin+1]*ibf.powers_[1];
                        res3 = int_seq[begin+2]*ibf.powers_[2];
                        res4 = int_seq[begin+3]*ibf.powers_[3];
                        res1 += int_seq[begin+4]*ibf.powers_[4];
                        res2 += int_seq[begin+5]*ibf.powers_[5];
                        ibf.forward_store_ = res1 + res2 + res3 + res4;
                        break;
                    case 7:
                        res1 = int_seq[begin+0]*ibf.powers_[0];
                        res2 = int_seq[begin+1]*ibf.powers_[1];
                        res3 = int_seq[begin+2]*ibf.powers_[2];
                        res4 = int_seq[begin+3]*ibf.powers_[3];
                        res1 += int_seq[begin+4]*ibf.powers_[4];
                        res2 += int_seq[begin+5]*ibf.powers_[5];
                        res3 += int_seq[begin+6]*ibf.powers_[6];
                        ibf.forward_store_ = res1 + res2 + res3 + res4;
                        break;
                    case 10:
                        res1 = int_seq[begin+0]*ibf.powers_[0];
                        res2 = int_seq[begin+1]*ibf.powers_[1];
                        res3 = int_seq[begin+2]*ibf.powers_[2];
                        res4 = int_seq[begin+3]*ibf.powers_[3];
                        res1 += int_seq[begin+4]*ibf.powers_[4];
                        res2 += int_seq[begin+5]*ibf.powers_[5];
                        res3 += int_seq[begin+6]*ibf.powers_[6];
                        res4 += int_seq[begin+7]*ibf.powers_[7];
                        res1 += int_seq[begin+8]*ibf.powers_[8];
                        res2 += int_seq[begin+9]*ibf.powers_[9];
                        ibf.forward_store_ = res1 + res2 + res3 + res4;
                        break;
                    case 14:
                        res1 = int_seq[begin+0]*ibf.powers_[0];
                        res2 = int_seq[begin+1]*ibf.powers_[1];
                        res3 = int_seq[begin+2]*ibf.powers_[2];
                        res4 = int_seq[begin+3]*ibf.powers_[3];
                        res1 += int_seq[begin+4]*ibf.powers_[4];
                        res2 += int_seq[begin+5]*ibf.powers_[5];
                        res3 += int_seq[begin+6]*ibf.powers_[6];
                        res4 += int_seq[begin+7]*ibf.powers_[7];
                        res1 += int_seq[begin+8]*ibf.powers_[8];
                        res2 += int_seq[begin+9]*ibf.powers_[9];
                        res3 += int_seq[begin+10]*ibf.powers_[10];
                        res4 += int_seq[begin+11]*ibf.powers_[11];
                        res1 += int_seq[begin+12]*ibf.powers_[12];
                        res2 += int_seq[begin+13]*ibf.powers_[13];
                        ibf.forward_store_ = res1 + res2 + res3 + res4;
                        break;
                    default:
                        for(size_t i = begin; i < end; i++)
                            ibf.forward_store_ += int_seq[i]*ibf.powers_[i-begin];
                        break;
                }
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
