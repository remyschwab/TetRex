#include "dGramIndex.h"


    std::vector<std::string> DGramTools::read_input_file_list(const std::filesystem::path &input_file)
    {
        std::vector<std::string> sequence_files;
        std::ifstream fin{input_file};

        if (!fin.good() || !fin.is_open())
            throw std::runtime_error{"Could not open file " + input_file.string() + " for reading."};

        std::string line;
        while (std::getline(fin, line))
        {
            sequence_files.emplace_back(line);
        }
        return sequence_files;
    }

void drive_dindex(const dindex_arguments &cmd_args)
{
    std::vector<std::string> input_bin_files;
    for (auto file : cmd_args.acid_libs)
    {
        if (file.extension() == ".lst")
        {
            for (auto f : DGramTools::read_input_file_list(file))
                input_bin_files.push_back(f);
        } 
        else
        {
            input_bin_files.push_back(std::filesystem::absolute(file));
        }
    }
    DGramIndex dindex(cmd_args.min, cmd_args.max, cmd_args.pad, cmd_args.hash_count, cmd_args.fpr, input_bin_files);
    dindex.populate();
    save_dindex(dindex, cmd_args.ofile);
}
