#include "mash.h"


void collect_kmers(std::vector<uint64_t> &kmers)
{

}


std::vector<uint64_t> construction_wrapper(std::string const &query, uint8_t const &k)
{
    nfa_t NFA;
    lmap_t nfa_map(NFA);
    amap_t arc_map;
    std::vector<uint64_t> query_kmers;

    construct_kgraph(query, NFA, nfa_map, arc_map, k);
    collect_kmers(query_kmers);
    return query_kmers;
}


void drive_mash(mash_arguments const &cmd_args)
{
    std::ifstream pattern_file_stream(cmd_args.pattern_file);
    std::vector<std::string> pattern_vector;
    std::string line_pattern;
    while(pattern_file_stream >> line_pattern)
    {
        pattern_vector.push_back(line_pattern);
    }
    std::vector<std::string> postfix_vector;
    for(auto p: pattern_vector)
    {
        postfix_vector.push_back(translate(p));
    }
    for(auto c: postfix_vector) std::cout << c << std::endl;
}
