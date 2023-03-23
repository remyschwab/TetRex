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
struct fastq_to_fasta : public cli_test {};

TEST_F(fastq_to_fasta, no_options)
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
