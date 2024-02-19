#pragma once


#include <string>
#include <cstddef>
#include <cstdint>
#include <type_traits>
#include "nucleotide_decomposer.h"
#include "peptide_decomposer.h"

namespace molecules
{
    using nucleotide = NucleotideDecomposer;
    using peptide = PeptideDecomposer;

    template<typename molecule_t>
    concept is_dna = std::same_as<molecule_t, molecules::nucleotide>;

    template<typename molecule_t>
    concept is_peptide = std::same_as<molecule_t, molecules::peptide>;

    template<typename molecule_t>
    concept is_molecule = is_dna<molecule_t> || is_peptide<molecule_t>;
}

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


template<molecules::is_molecule DecomposerType>
class MoleculeDecomposer
{
    private:
        uint8_t ksize_;

    public:
        uint8_t lshift_;
        uint8_t rmask_;
        DecomposerType decomposer_;
        static bool constexpr is_nucleotide{molecules::is_dna<DecomposerType>};
        MoleculeDecomposer() = default;
        explicit MoleculeDecomposer(uint8_t &ksize, uint8_t &reduction) : ksize_{ksize}
        {
            if constexpr(is_nucleotide)
            {
                decomposer_ = molecules::NucleotideDecomposer(ksize_, reduction);
                lshift_ = 2;
                rmask_ = 0b11;
            }
            else
            {
                decomposer_ = molecules::PeptideDecomposer(ksize_, reduction);
                lshift_ = 5;
                rmask_ = 0b11111;
            }
        }

    void decompose_record(std::string_view record, seqan::hibf::bin_index &tech_bin_id, auto &base_ref)
    {
        decomposer_.decompose_record(record, tech_bin_id, base_ref);
    }

    uint64_t update_kmer(const int &symbol, uint64_t &kmer)
    {
        return decomposer_.update_kmer(symbol, kmer);
    }

    void create_selection_bitmask()
    {
        decomposer_.create_selection_bitmask();
    }

    void print_mask()
    {
        decomposer_.print_mask();
    }

    template<class Archive>
    void serialize(Archive &archive)
    {
        archive(ksize_, lshift_, rmask_, decomposer_);
    }
};
