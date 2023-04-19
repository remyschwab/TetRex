#include "molecule_encodings.h"


uint64_t encode_dna(std::string_view kmer)
{
    uint64_t codemer = 0;
    for(auto && base: kmer)
    {
        codemer = codemer << 2;
        // codemer += std::popcount(-33 & (int)base - 65);
        codemer += (base>>1)&3;
        
    }
    // uint64_t revcomp = revComplement(codemer, kmer.length());
    return codemer;
    // return codemer <= revcomp ? codemer : revcomp;
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
