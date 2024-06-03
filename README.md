# Search for Regular Expressions in large datasets
Despite the efficiency of modern day tools for Regular Expression search, their runtime is often dominated by the size of the text. We present TetRex, a novel algorithm for regular expression matching that leverages the (Hierarchical) Interleaved Bloom Filter as an index. Regular Expressions are given as input in the command line and support the following operations:

1. **|** - Or
2. __*__ - Zero or more repetitions
3. **+** - One or More repetitions
4. **?** - Optional Character


## Installation

1. Clone the repository with
```git clone --recurse-submodules git@github.com:remyschwab/TetRex.git```
2. Descend into the home directory and input:
```mkdir build && cd build```
3. Configure with cmake ```cmake -DCMAKE_CXX_COMPILER=/path/to/g++-11 ..```
4. Build with make ```make```

## Usage
Tetrix offers two main commands [index & query]:
```bash
## Index Nucleic Acid DB, with kmer size = 7, 3 hash functions, & 644830 bits per Bloom Filter, each input file represents a bin
tetrex index -m na -k 3 -o test data/dna_example_split/*.fa
## Query RegEx
tetrex query test.ibf "A(C+|G+)T" 
#TetRex/data/dna_example_split/sequence1.fa >Sequence1      ACT
#TetRex/data/dna_example_split/sequence1.fa >Sequence1      ACT
#TetRex/data/dna_example_split/sequence1.fa >Sequence1      ACT
#TetRex/data/dna_example_split/sequence2.fa >Sequence2      ACT
#TetRex/data/dna_example_split/sequence2.fa >Sequence2      AGT
#TetRex/data/dna_example_split/sequence4.fa >Sequence4      ACCCT
```

## Notes
This app was generated from the [SeqAn App Template](https://github.com/seqan/app-template) and makes heavy use of the [SeqAn library](https://github.com/seqan/seqan3/tree/4668203ee1526b4ac3dbdc47869bee72253f684c).
