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
Tetrix offers two main commands [index & query] and one utility command [inspect]:
```bash
## Construct an HIBF index of nucleic acids, with kmer size = 3, 3 hash functions, & a FPR of 0.05, each input file represents a bin
tetrex index -m na -k 3 -o test data/dna_example_split/*.fa
## Query RegEx
tetrex query test.ibf "A(C+|G+)T" 
#TetRex/data/dna_example_split/sequence1.fa >Sequence1      ACT
#TetRex/data/dna_example_split/sequence1.fa >Sequence1      ACT
#TetRex/data/dna_example_split/sequence1.fa >Sequence1      ACT
#TetRex/data/dna_example_split/sequence2.fa >Sequence2      ACT
#TetRex/data/dna_example_split/sequence2.fa >Sequence2      AGT
#TetRex/data/dna_example_split/sequence4.fa >Sequence4      ACCCT
## Report some info about a constructed index
tetrex inspect test.ibf
# Reading Index from Disk... DONE in 3.1e-05s
# INDEX TYPE: HIBF
# FALSE POSITIVE RATE: 0.05
# HASH COUNT (hash functions): 3
# KMER LENGTH (bases): 3
# MOLECULE TYPE (alphabet): Nucleic Acid [REDUCTION=NONE]
# ACID LIBRARY (filepaths):
#         - /Users/rschwab/repos/TetRex/./data/dna_example_split/sequence1.fa
#         - /Users/rschwab/repos/TetRex/./data/dna_example_split/sequence2.fa
#         - /Users/rschwab/repos/TetRex/./data/dna_example_split/sequence3.fa
#         - /Users/rschwab/repos/TetRex/./data/dna_example_split/sequence4.fa
#         - /Users/rschwab/repos/TetRex/./data/dna_example_split/sequence5.fa
# DONE
```

## Notes
This app was generated from the [SeqAn App Template](https://github.com/seqan/app-template) and makes heavy use of the [SeqAn library](https://github.com/seqan/seqan3/tree/4668203ee1526b4ac3dbdc47869bee72253f684c).
