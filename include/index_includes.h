#ifndef index_includes_h
#define index_includes_h

#include <zlib.h>
#include <cereal/types/vector.hpp>
#include "cereal/types/string.hpp"

#include "utils.h"
#include "arg_parse.h"
#include "kseq.h"
#include "robin_hood.h"

#include "molecule_decomposer.h"
#include "index_ibf.h"
#include "index_hibf.h"


namespace index_structure
{
    using SeqAnIBF = seqan::hibf::interleaved_bloom_filter;
    using IBF = IBFIndex;
    using IBF_Agent = seqan::hibf::interleaved_bloom_filter::membership_agent_type;
    using SeqAnHIBF = seqan::hibf::hierarchical_interleaved_bloom_filter;
    using HIBF = HIBFIndex;
    using HIBF_Agent = seqan::hibf::hierarchical_interleaved_bloom_filter::membership_agent_type;

    template <typename index_t>
    concept is_ibf = std::same_as<index_t, IBF>;

    template <typename index_t>
    concept is_hibf = std::same_as<index_t, HIBF>;

    template <typename index_t>
    concept is_valid = is_ibf<index_t> || is_hibf<index_t>;
} // namespace index_structure

#endif
