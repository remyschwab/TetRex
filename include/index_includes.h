#pragma once

#include <zlib.h>

#include <cereal/types/vector.hpp>

#include "utils.h"
#include "arg_parse.h"
#include "kseq.h"
#include "robin_hood.h"
#include "molecule_decomposer.h"

namespace index_structure
{
    using IBF = seqan::hibf::interleaved_bloom_filter;
    using IBF_Agent = seqan::hibf::interleaved_bloom_filter::membership_agent_type;
    using HIBF = seqan::hibf::hierarchical_interleaved_bloom_filter;
    using HIBF_Agent = seqan::hibf::hierarchical_interleaved_bloom_filter::membership_agent_type;

    template <typename index_t>
    concept is_ibf = std::same_as<index_t, index_structure::IBF>;

    template <typename index_t>
    concept is_hibf = std::same_as<index_t, index_structure::HIBF>;

    template <typename index_t>
    concept is_valid = is_ibf<index_t> || is_hibf<index_t>;
} // namespace index_structure
