#pragma once

#include <iostream>
#include <type_traits>
#include <string>
#include <vector>
#include <algorithm>
#include <math.h>
#include <fstream>
#include <stack>
#include <filesystem>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <cstring>
#include <zlib.h>

#include "kseq.h"
#include "hibf/interleaved_bloom_filter.hpp"
#include "hibf/hierarchical_interleaved_bloom_filter.hpp"
#include "hibf/misc/bit_vector.hpp"
#include <seqan3/core/debug_stream.hpp>

KSEQ_INIT(gzFile, gzread)

/////////////// Type Declarations ///////////////
using bitvector = seqan::hibf::bit_vector;
/////////////// ****** END ****** ///////////////

char* re2post(char *re);

std::string translate(const std::string& str);
