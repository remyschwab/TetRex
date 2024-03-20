#pragma once

#include <simde/x86/ssse3.h>

namespace molecules
{
    class NucleotideDecomposer
    {
        private:
            uint8_t k_{};
            uint8_t reduction_{};

        public:
            uint8_t left_shift_{};
            uint64_t selection_mask_{};

            NucleotideDecomposer() = default;
            NucleotideDecomposer(uint8_t const k, uint8_t const reduction) : k_{k}, reduction_{reduction}
            {
                create_selection_bitmask();
                set_left_shift();
            }

            void set_left_shift()
            {
                left_shift_ = ((uint8_t)2 * k_)-2;
            }

            void create_selection_bitmask()
            {
                /*
                Example with k=4 creates a bitmask like 
                0b00000000-00000000-00000000-00000000-00000000-00000000-00000000-11111111
                for the rolling hash
                */
                selection_mask_ = (k_ >= 32) ? static_cast<uint64_t>(-1) : (1ULL << (2*(k_))) - 1u;
            }

            void print_mask() const
            {
                std::cout << selection_mask_ << std::endl;
            }

            uint64_t revComplement(const uint64_t kmer, const int k)
            {
                // broadcast 64bit to 128 bit
                simde__m128i x = simde_mm_cvtsi64_si128(kmer);

                // create lookup (set 16 bytes in 128 bit)
                // a lookup entry at the index of two nucleotides (4 bit) describes the reverse
                // complement of these two nucleotides in the higher 4 bits (lookup1) or in the
                // lower 4 bits (lookup2)
                simde__m128i lookup1 = simde_mm_set_epi8(0x50,0x10,0xD0,0x90,0x40,0x00,0xC0,0x80,0x70,
                                                        0x30,0xF0,0xB0,0x60,0x20,0xE0,0xA0);
                simde__m128i lookup2 = simde_mm_set_epi8(0x05,0x01,0x0D,0x09,0x04,0x00,0x0C,0x08,0x07,
                                                        0x03,0x0F,0x0B,0x06,0x02,0x0E,0x0A);
                // set upper 8 bytes to 0 and revert order of lower 8 bytes
                simde__m128i upper = simde_mm_set_epi8(0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0,1,2,3,4,5,6,7);

                simde__m128i kmer1 = simde_mm_and_si128(x, simde_mm_set1_epi8(0x0F)); // get lower 4 bits
                simde__m128i kmer2 = simde_mm_and_si128(x, simde_mm_set1_epi8(0xF0)); // get higher 4 bits

                // shift right by 2 nucleotides
                kmer2 >>= 4;

                // use _mm_shuffle_epi8 to look up reverse complement
                kmer1 = simde_mm_shuffle_epi8(lookup1, kmer1);
                kmer2 = simde_mm_shuffle_epi8(lookup2, kmer2);

                // _mm_or_si128: bitwise OR
                x = simde_mm_or_si128(kmer1, kmer2);

                // set upper 8 bytes to 0 and revert order of lower 8 bytes
                x = simde_mm_shuffle_epi8(x, upper);

                // shift out the unused nucleotide positions (1 <= k <=32 )
                // broadcast 128 bit to 64 bit
                return (((uint64_t)simde_mm_cvtsi128_si64(x)) >> (uint64_t)(64-2*k));
            }

            uint64_t encode_dna(std::string_view const kmer) const
            {
                uint64_t codemer{};
                for(auto && base: kmer)
                {
                    codemer <<= 2;
                    codemer |= (base>>1)&3;
                }
                return codemer;
            }

            void rollover_nuc_hash(const char base , seqan::hibf::bin_index const bin_id, auto &base_ref) const
            {
                auto fb = (base>>1)&3; // Encode the new base
                auto cb = (fb^0b10)<<left_shift_; // Get its complement and shift it to the big end of the uint64
                base_ref.forward_store_ = ((base_ref.forward_store_<<2)&selection_mask_) | fb; // Update forward store 
                base_ref.reverse_store_ = ((base_ref.reverse_store_>>2)&selection_mask_) | cb; // Update reverse store
                uint64_t canon_kmer = ( base_ref.forward_store_ <= base_ref.reverse_store_ ? base_ref.forward_store_ : base_ref.reverse_store_ );
                base_ref.emplace(canon_kmer, bin_id);
            }

            void decompose_record(std::string_view const record_seq, seqan::hibf::bin_index const tech_bin_id, auto &base_ref)
            {
                uint64_t initial_encoding = encode_dna(record_seq.substr(0, k_)); // Encode forward
                uint64_t reverse_complement = revComplement(initial_encoding, k_); // Compute the reverse compelement
                base_ref.set_stores(initial_encoding, reverse_complement); // Remember both strands
                base_ref.emplace((initial_encoding <= reverse_complement ? initial_encoding : reverse_complement), tech_bin_id); // MinHash
                for(auto symbol : record_seq)
                {
                    rollover_nuc_hash(symbol, tech_bin_id, base_ref);
                }
            }

            uint64_t update_kmer(const int &symbol, uint64_t &kmer)
            {
                uint64_t fb = (symbol>>1)&3;
                uint64_t forward = ((kmer<<2)&selection_mask_) | fb;
                kmer = forward;
                uint64_t reverse = revComplement(forward, k_);
                uint64_t canon_kmer = forward <= reverse ? forward : reverse;
                return canon_kmer;
            }
    
            template<class Archive>
            void serialize(Archive &archive)
            {
                archive(k_, reduction_, left_shift_, selection_mask_);
            }
    };
}