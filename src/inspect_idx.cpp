#include "inspect_idx.h"



void drive_inspection(const inspection_arguments &cmd_args)
{
    double t1, t2;
    // Load index from disk
    seqan3::debug_stream << "Reading Index from Disk... ";
    IndexStructure ibf;
    t1 = omp_get_wtime();
    load_ibf(ibf, cmd_args.idx);
    t2 = omp_get_wtime();
    seqan3::debug_stream << "DONE in " << t2-t1 << "s" << std::endl;

    // Print Info
    std::cout << "BIN COUNT (BFs): " << ibf.getBinCount() << std::endl;
    std::cout << "BIN SIZE (bits): " << ibf.getBinSize() << std::endl;
    std::cout << "HASH COUNT (hash functions): " << unsigned(ibf.getHashCount()) << std::endl;
    std::cout << "KMER LENGTH (bases): " << unsigned(ibf.k_) << std::endl;
    std::string acid_type;
    if(ibf.molecule_ == "na")
    {
        acid_type = "Nucleic Acid";
    } else
    {
        acid_type = "Amino Acid";
    }
    std::cout << "MOLECULE TYPE (alphabet): " << acid_type << std::endl;
    std::cout << "ACID LIBRARY (filepaths): " << std::endl;
    for(auto && path: ibf.acid_libs_)
        std::cout << "\t- " << path << std::endl;

    seqan3::debug_stream << "DONE" << std::endl;
}
