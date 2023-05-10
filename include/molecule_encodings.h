#pragma once

#include <string>
#include <cstddef>
#include <cstdint>

#include <robin_hood.h>
#include <simde/x86/ssse3.h>
#include <seqan3/core/debug_stream.hpp>


/*
    ASCII TABLE
    
    NUCLEOTIDES: (base>>1)&3 also try std::popcount(-33 & base-65)
        - A	065	01000-00-1
        - C	067	01000-01-1
        - G	071	01000-11-1
        - T	084	01010-10-0
        - U	085	01010-10-1
    
    AMINO ACIDS:
        - A	065	01000001 Alanine
        - B	066	01000010 Aspartic Acid or Asparagine
        - C	067	01000011 Cysteine
        - D	068	01000100 Aspartic Acid
        - E	069	01000101 Glutamic Acid
        - F	070	01000110 Phenylalanine
        - G	071	01000111 Glycine
        - H	072	01001000 Histidine
        - I	073	01001001 Isoleucine
        - J	074	01001010 Leucine or Isoleucine
        - K	075	01001011 Lysine
        - L	076	01001100 Leucine
        - M	077	01001101 Methionine
        - N	078	01001110 Asparagine
        - O	079	01001111 Pyrrolysine
        - P	080	01010000 Proline
        - Q	081	01010001 Glutamine
        - R	082	01010010 Arginine
        - S	083	01010011 Serine
        - T	084	01010100 Threonine
        - U	085	01010101 Selenocysteine
        - V	086	01010110 Valine
        - W	087	01010111 Tryptophan
        - X	088	01011000 ANY
        - Y	089	01011001 Tyrosine
        - Z	090	01011010 Glutamic Acid or Glutamine
*/


uint64_t encode_dna(std::string_view kmer);

uint64_t revComplement(const uint64_t kmer, const int k);

void create_residue_maps(uint8_t &alphabet_size, std::vector<uint8_t> &aamap);
