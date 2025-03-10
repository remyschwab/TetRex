#pragma once

#include "index_includes.h"

struct index_params
{
    uint8_t k_{};
    std::string molecule_{};
    bool is_hibf_{};

    template<class Archive>
    void serialize(Archive &archive)
    {
        archive(k_, molecule_, is_hibf_);
    }
};


template <index_structure::is_valid ibf_flavor, molecules::is_molecule molecule_t>
class TetrexIndex
{
    public:
        uint8_t k_{};
        std::string molecule_{};
        bool permanent_hibf_status_{};
        std::vector<std::string> acid_libs_{};
        uint8_t reduction_{};
        ibf_flavor ibf_{};
        MoleculeDecomposer<molecule_t> decomposer_{};
        static bool constexpr is_hibf_{index_structure::is_hibf<ibf_flavor>};
        uint64_t forward_store_{};
        uint64_t reverse_store_{};


        TetrexIndex() = default;

        explicit TetrexIndex(uint8_t k,
                            float fpr,
                            uint8_t hc,
                            std::string molecule,
                            std::vector<std::string> acid_libs,
                            uint8_t reduction) requires is_hibf_:
            k_{k},
            molecule_{std::move(molecule)},
            acid_libs_{std::move(acid_libs)},
            reduction_{reduction},
            decomposer_(k_, reduction_),
            permanent_hibf_status_{is_hibf_}
        {
            permanent_hibf_status_ = true;
            ibf_ = HIBFIndex(fpr, hc, acid_libs_);
        }

        explicit TetrexIndex(uint8_t k,
                            float fpr,
                            uint8_t hc,
                            std::string molecule,
                            std::vector<std::string> acid_libs,
                            uint8_t reduction) requires (!is_hibf_):
            k_{k},
            molecule_{std::move(molecule)},
            acid_libs_{std::move(acid_libs)},
            reduction_{reduction},
            decomposer_(k_, reduction_),
            permanent_hibf_status_{is_hibf_}
        {
            permanent_hibf_status_ = false;
            size_t bin_count = acid_libs_.size();
            ibf_ = IBFIndex(bin_count, fpr, hc, acid_libs_);
        }

        int calculate_m(size_t n, float p)
        {
            double numerator = n * std::log(p);
            double denominator = std::log(1.0 / std::pow(2, std::log(2)));
            return std::ceil(numerator / denominator);
        }

        std::pair<size_t, size_t> getShape() const
        {
            return ibf_.getShape();
        }

        size_t getBinCount() const
        {
            return ibf_.getBinCount();
        }

        size_t getHashCount() const
        {
            return ibf_.getHashCount();
        }

        void populate_index()
        {
            ibf_.populate_index(k_, decomposer_, *this);
        }

        void emplace(uint64_t const val, size_t const idx)
        {
            ibf_.emplace(val, idx);
        }

        bitvector query(uint64_t const kmer)
        {
            return ibf_.query(kmer); // Hmmm ask about this
        }

        bool is_hibf() const
        {
            return is_hibf_;
        }

        ibf_flavor & getIBF()
        {
            return ibf_;
        }

        float getFPR() const
        {
            return ibf_.getFPR();
        }

        void set_lib_paths(std::vector<std::string> path_collection)
        {
            acid_libs_ = std::move(path_collection);
        }

        void set_stores(uint64_t const forward, uint64_t const reverse)
        {
            forward_store_ = forward;
            reverse_store_ = reverse;
        }

        std::tuple<uint64_t, uint64_t> get_stores() const
        {
            return {forward_store_, reverse_store_};
        }

        uint64_t update_kmer(int const symbol, uint64_t &kmer)
        {
            return decomposer_.update_kmer(symbol, kmer);
        }

        void spawn_agent()
        {
            ibf_.spawn_agent();
        }

        void create_selection_bitmask()
        {
            decomposer_.create_selection_bitmask();
        }

        void print_mask() const
        {
            decomposer_.print_mask();
        }

        template<class Archive>
        void serialize(Archive &archive)
        {
            archive(k_, molecule_, permanent_hibf_status_);
            archive(acid_libs_, reduction_, ibf_, decomposer_);
        }
        
}; // TetrexIndex

std::vector<std::string> read_input_file_list(const std::filesystem::path &input_file);

void create_ibf_dna_index(const index_arguments &cmd_args, const std::vector<std::string> &input_bin_files, uint8_t const reduction);

void create_hibf_dna_index(const index_arguments &cmd_args, const std::vector<std::string> &input_bin_files, uint8_t const reduction);

void create_ibf_aa_index(const index_arguments &cmd_args, const std::vector<std::string> &input_bin_files, uint8_t const reduction);

void create_hibf_aa_index(const index_arguments &cmd_args, const std::vector<std::string> &input_bin_files, uint8_t const reduction);

void drive_index(const index_arguments &cmd_args);

template <class TetrexIndex>
void store_ibf(TetrexIndex const &ibf, std::filesystem::path const &opath)
{
    std::ofstream os{opath, std::ios::binary};
    cereal::BinaryOutputArchive oarchive{os};
    oarchive(ibf);
}

template <class TetrexIndex>
void load_ibf(TetrexIndex &ibf, const std::filesystem::path &ipath)
{
    std::ifstream is{ipath, std::ios::binary};
    cereal::BinaryInputArchive iarchive{is};
    iarchive(ibf);
}

inline void load_params(index_params &params, const std::filesystem::path &ipath)
{
    std::ifstream is{ipath, std::ios::binary};
    cereal::BinaryInputArchive iarchive{is};
    iarchive(params);
}
