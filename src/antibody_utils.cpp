#include "antibody_utils.h"


void drive_antibody_index(const antibody_index_arguments &cmd_args)
{
    std::vector<std::string> input_npz_paths;
    for(auto path: cmd_args.canonical_alignments)
    {
        input_npz_paths.push_back(std::filesystem::absolute(path));
    }
}