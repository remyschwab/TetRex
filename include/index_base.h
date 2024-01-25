#pragma once

#include "index_ibf.h"
#include "index_hibf.h"
#include "index_includes.h"


template <index_structure::is_valid ibf_flavor, molecules::is_molecule molecule_t>
class TetrexIndex
{
    private:
        std::vector<std::vector<std::string>> bin_path_;
        double fpr_;
        constexpr bool is_hibf_{index_structure::is_hibf<ibf_flavor>};
        ibf_flavor ibf_;
        molecule_t decomposer_;

    public:
        uint8_t k_;
        std::string molecule_;
        std::vector<std::string> acid_libs_;
        std::string reduction_;
        uint64_t forward_store_{};
        uint64_t reverse_store_{};


        TetrexIndex() = default;

        explicit TetrexIndex(uint8_t k,
                            size_t bin_size,
                            uint8_t hc,
                            std::string molecule,
                            std::vector<std::string> acid_libs,
                            std::string reduction) :
            k_{k},
            molecule_{molecule},
            acid_libs_{acid_libs},
            reduction_{reduction},
            MoleculeDecomposer<molecule_t> decompser_(k_, reduction_)
        {
            if constexpr(!is_hibf_)
            {
                IBFIndex ibf_(bin_size, hc, acid_libs_);
            }
            else
            {
                HIBFIndex ibf_(bin_size, hc, acid_libs_);
            }
        }

        void populate_index()
        {
            ibf_.populate_index(k_, decomposer_);
        }

        void emplace(uint64_t &val, size_t &idx)
        {
            ibf_.emplace(val, idx);
        }

        const bitvector & query(uint64_t &kmer)
        {
            return ibf_.query(kmer);
        }

        std::vector<std::string> const & bin_path() const
        {
            return bin_path_;
        }

        // double fpr() const
        // {
        //     return fpr_;
        // }

        bool is_hibf() const
        {
            return is_hibf_;
        }

        ibf_flavor & ibf()
        {
            return ibf_;
        }

        void set_lib_paths(std::vector<std::string> path_collection)
        {
            acid_libs_ = path_collection;
        }

        uint8_t map_aa(unsigned char &residue)
        {
            return aamap_[residue];
        }

        void set_stores(uint64_t forward, uint64_t reverse)
        {
            forward_store_ = forward;
            reverse_store_ = reverse;
        }

        std::tuple<uint64_t, uint64_t> get_stores()
        {
            return {forward_store_, reverse_store_};
        }

        // template <seqan3::cereal_archive archive_t>
        // void CEREAL_SERIALIZE_FUNCTION_NAME(archive_t & archive)
        // {
        //     uint32_t parsed_version{TetrexIndex<>::version};
        //     archive(parsed_version);
        //     if (parsed_version == TetrexIndex<>::version)
        //     {
        //         try
        //         {
        //             archive(window_size_);
        //             archive(shape_);
        //             archive(parts_);
        //             archive(bin_path_);
        //             archive(fpr_);
        //             archive(is_hibf_);
        //             archive(ibf_);
        //         }
        //         // GCOVR_EXCL_START
        //         catch (std::exception const & e)
        //         {
        //             throw sharg::parser_error{"Cannot read index: " + std::string{e.what()}};
        //         }
        //         // GCOVR_EXCL_STOP
        //     }
        //     else
        //     {
        //         throw sharg::parser_error{"Unsupported index version. Check raptor upgrade."}; // GCOVR_EXCL_LINE
        //     }
        // }

        // template <seqan3::cereal_input_archive archive_t>
        // void load_parameters(archive_t & archive)
        // {
        //     uint32_t parsed_version{};
        //     archive(parsed_version);
        //     if (parsed_version == version)
        //     {
        //         try
        //         {
        //             archive(window_size_);
        //             archive(shape_);
        //             archive(parts_);
        //             archive(bin_path_);
        //             archive(fpr_);
        //             archive(is_hibf_);
        //         }
        //         // GCOVR_EXCL_START
        //         catch (std::exception const & e)
        //         {
        //             throw sharg::parser_error{"Cannot read index: " + std::string{e.what()}};
        //         }
        //         // GCOVR_EXCL_STOP
        //     }
        //     else
        //     {
        //         throw sharg::parser_error{"Unsupported index version. Check raptor upgrade."}; // GCOVR_EXCL_LINE
        //     }
        // }

        template<class Archive>
        void serialize(Archive &archive)
        {
            archive(ibf_, decomposer_, k_, molecule_, acid_libs_, reduction_);
        }
}; // TetrexIndex

void read_input_file_list(std::vector<std::filesystem::path> & sequence_files, std::filesystem::path input_file);

void create_ibf_dna_index(const index_arguments &cmd_args, const std::vector<std::string> &input_bin_files);

void create_hibf_dna_index(const index_arguments &cmd_args, const std::vector<std::string> &input_bin_files);

void create_ibf_aa_index(const index_arguments &cmd_args, const std::vector<std::string> &input_bin_files);

void create_hibf_aa_index(const index_arguments &cmd_args, const std::vector<std::string> &input_bin_files);

void drive_index(const index_arguments &cmd_args);

template <class IndexStructure>
void store_ibf(IndexStructure const & ibf, std::filesystem::path opath)
{
    std::ofstream os{opath, std::ios::binary};
    cereal::BinaryOutputArchive oarchive{os};
    oarchive(ibf);
}

template <class IndexStructure>
void load_ibf(IndexStructure & ibf, std::filesystem::path ipath)
{
    std::ifstream is{ipath, std::ios::binary};
    cereal::BinaryInputArchive iarchive{is};
    iarchive(ibf);
}
