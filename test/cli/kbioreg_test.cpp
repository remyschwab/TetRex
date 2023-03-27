#include <ranges>     // range comparisons
#include <string>     // strings
#include <vector>     // vectors

#include <seqan3/alphabet/detail/debug_stream_alphabet.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/core/debug_stream/tuple.hpp>
#include <seqan3/io/sequence_file/input.hpp>

// Include the EXPECT_RANGE_EQ macro for better information if range elements differ.
#include <seqan3/test/expect_range_eq.hpp>

#include "cli_test.hpp"

// To prevent issues when running multiple CLI tests in parallel, give each CLI test unique names:
struct kbioreg : public cli_test {};

TEST_F(kbioreg, no_options)
{
    cli_test_result result = execute_app("kbioreg");
    std::string expected
    {
        "kbioreg\n"
        "=======\n"
        "    Try -h or --help for more information.\n"
    };
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(kbioreg, build_index)
{
    cli_test_result result = execute_app("kbioreg", "index", "-k", "3", "-c", "3", "-s", "64", "-m", "na",
                                         "-o", "ibf_idx", data("file1.fa"), data("file2.fa"));
    std::string expected
    {
        "Indexed 4 sequences across 2 bins.\n"
        "Writing to disk... DONE\n"
    };
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected);
}

TEST_F(kbioreg, query)
{
    std::system("ln -s ../.. test"); // !HACK !TODO how to link to the original path?
    cli_test_result result = execute_app("kbioreg", "query", data("ibf_idx.ibf"),
                                         "\"(AC+G)\"");
    std::string expected
    {
        ">Snippet1.1\tACCG\n"
        ">Snippet1.2\tACG\n"
    };
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, expected);
//    EXPECT_EQ(result.err, std::string{}); // !TODO how to deal with timing outputs
}
