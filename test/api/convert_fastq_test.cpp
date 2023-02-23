#include <gtest/gtest.h>

#include <seqan3/alphabet/nucleotide/dna5.hpp>




TEST(group1, out_empty)
{
    std::vector<seqan3::dna5> vector;
    EXPECT_EQ(vector.size(), 0);
//    std::string expected{"> seq1\nACGTTTGATTCGCG\n> seq2\nTCGGGGGATTCGCG\n"};
//    testing::internal::CaptureStdout();
//    convert_fastq(DATADIR"in.fastq", "");
//    std::string std_cout = testing::internal::GetCapturedStdout();
//    EXPECT_RANGE_EQ(expected, std_cout);
}
