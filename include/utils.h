#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
#include <stack>
#include <vector>
#include <filesystem>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <cstring>


#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>




/////////////// Type Declarations ///////////////
using bitvector = seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed>::membership_agent_type::binning_bitvector;
/////////////// ****** END ****** ///////////////

char* re2post(char *re);

std::string translate(const std::string& str);
