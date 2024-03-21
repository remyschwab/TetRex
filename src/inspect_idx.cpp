#include "inspect_idx.h"


void inspect_dna_ibf(inspection_arguments const &cmd_args)
{
    double t1, t2;
    // Load index from disk
    seqan3::debug_stream << "Reading Index from Disk... ";
    TetrexIndex<index_structure::IBF, molecules::nucleotide> ibf;
    t1 = omp_get_wtime();
    load_ibf(ibf, cmd_args.idx);
    t2 = omp_get_wtime();
    seqan3::debug_stream << "DONE in " << t2-t1 << "s" << std::endl;

    std::pair<size_t, size_t> shape = ibf.getShape();

    // Print Info to std out
    std::cout << "BIN COUNT (BFs): " << shape.first << std::endl;
    std::cout << "BIN SIZE (bits): " << shape.second << std::endl;
    std::cout << "HASH COUNT (hash functions): " << unsigned(ibf.getHashCount()) << std::endl;
    std::cout << "KMER LENGTH (bases): " << unsigned(ibf.k_) << std::endl;
    std::cout << "MOLECULE TYPE (alphabet): Nucleic Acid [REDUCTION=" << ibf.reduction_ << "]" << std::endl;
    std::cout << "ACID LIBRARY (filepaths):" << std::endl;
    for(auto && path: ibf.acid_libs_)
        std::cout << "\t- " << path << std::endl;

    seqan3::debug_stream << "DONE" << std::endl;
}


void inspect_dna_hibf(inspection_arguments const &cmd_args)
{
    (void)cmd_args;
    return;
}


void inspect_aa_ibf(inspection_arguments const &cmd_args)
{
    double t1, t2;
    // Load index from disk
    seqan3::debug_stream << "Reading Index from Disk... ";
    TetrexIndex<index_structure::IBF, molecules::nucleotide> ibf;
    t1 = omp_get_wtime();
    load_ibf(ibf, cmd_args.idx);
    t2 = omp_get_wtime();
    seqan3::debug_stream << "DONE in " << t2-t1 << "s" << std::endl;

    std::pair<size_t, size_t> shape = ibf.getShape();

    // Print Info to std out
    std::cout << "BIN COUNT (BFs): " << shape.first << std::endl;
    std::cout << "BIN SIZE (bits): " << shape.second << std::endl;
    std::cout << "HASH COUNT (hash functions): " << unsigned(ibf.getHashCount()) << std::endl;
    std::cout << "KMER LENGTH (bases): " << unsigned(ibf.k_) << std::endl;
    std::cout << "MOLECULE TYPE (alphabet): Amino Acid [REDUCTION=" << ibf.reduction_ << "]" << std::endl;
    std::cout << "ACID LIBRARY (filepaths):" << std::endl;
    for(auto && path: ibf.acid_libs_)
        std::cout << "\t- " << path << std::endl;

    seqan3::debug_stream << "DONE" << std::endl;
}


void inspect_aa_hibf(inspection_arguments const &cmd_args)
{
    (void)cmd_args;
    return;
}


void drive_inspection(inspection_arguments const &cmd_args)
{
    index_params params;
    load_params(params, cmd_args.idx);
    if(!params.is_hibf_ && params.molecule_ == "na")
    {
        inspect_dna_ibf(cmd_args);
    }
    else if(params.is_hibf_ && params.molecule_ == "na")
    {
        inspect_dna_hibf(cmd_args);
    }
    else if(!params.is_hibf_ && params.molecule_ == "aa")
    {
        inspect_aa_ibf(cmd_args);
    }
    else if(params.is_hibf_ && params.molecule_ == "aa")
    {
        inspect_aa_hibf(cmd_args);
    }
}
