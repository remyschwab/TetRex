#include <array>

#include "seqan3/core/debug_stream/all.hpp"
#include "arg_parse.h"
#include "robin_hood.h"
#include "cnpy.h"


namespace Antibody_Utilities
{

    struct NumberedResidue
    {
        int position;        // IMGT position (e.g., 30, 52)
        char insertion_code; // 'A', 'B', etc. for 30A, 30B; or ' ' if none
        char aa;             // Amino acid
        std::string region;  // CDR1, FR1, etc. (optional for now)
    };
    
    struct CanonicalAlignment
    {
        std::array<uint8_t, 200> alignment_data;
    };
    
    
    std::array<uint8_t, 200> region_codes = {
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,       // FR1
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0,
        1, 1, 1, 1, 1, 1, 1, 1, 1,          // CDR1
        1, 1, 1, 1, 1, 1, 1, 1,
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2,       // FR2
        2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
        2, 2, 2, 2, 2, 2,
        3, 3, 3, 3, 3, 3, 3, 3, 3, 3,       // CDR2
        3, 3, 3, 3, 3, 3, 3, 3, 3,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4,       // FR3
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4,
        5, 5, 5, 5, 5, 5, 5, 5, 5, 5,       // CDR3
        5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
        5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
        5, 5, 5, 5, 5, 5, 5,
        6, 6, 6, 6, 6, 6, 6, 6, 6, 6,       // FR4
        6, 6
    };

}

void drive_antibody_index(const antibody_index_arguments &cmd_args);
