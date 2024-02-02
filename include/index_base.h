#pragma once

#include "index_includes.h"


template <index_structure::is_valid ibf_flavor, molecules::is_molecule molecule_t>
class TetrexIndex
{
    private:
        double fpr_;
        bool is_hibf_{index_structure::is_hibf<ibf_flavor>};
        ibf_flavor ibf_;
        MoleculeDecomposer<molecule_t> decomposer_;

    public:
        uint8_t k_;
        std::string molecule_;
        std::vector<std::string> acid_libs_;
        uint8_t reduction_;
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
            reduction_{index_structure::reduction_map[reduction]},
            decomposer_(k, reduction_)
        {
            if (!is_hibf_)
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
            ibf_.populate_index(k_, decomposer_, *this);
        }

        void emplace(uint64_t &val, const seqan::hibf::bin_index &idx)
        {
            ibf_.emplace(val, idx);
        }

        const bitvector & query(uint64_t &kmer)
        {
            return ibf_.query(kmer);
        }

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

        void set_stores(uint64_t forward, uint64_t reverse)
        {
            forward_store_ = forward;
            reverse_store_ = reverse;
        }

        std::tuple<uint64_t, uint64_t> get_stores()
        {
            return {forward_store_, reverse_store_};
        }

        void update_kmer(const int &symbol, uint64_t &kmer)
        {
            decomposer_.update_kmer(symbol, kmer);
        }

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

template <class TetrexIndex>
void store_ibf(TetrexIndex const & ibf, std::filesystem::path opath)
{
    std::ofstream os{opath, std::ios::binary};
    cereal::BinaryOutputArchive oarchive{os};
    oarchive(ibf);
}

template <class TetrexIndex>
void load_ibf(TetrexIndex &ibf, std::filesystem::path ipath)
{
    std::ifstream is{ipath, std::ios::binary};
    cereal::BinaryInputArchive iarchive{is};
    iarchive(ibf);
}
