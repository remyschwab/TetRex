# Search for Regular Expressions in large datasets
Despite the many efficient tools implemented for finding Regular Expressions, the runtime needed to search for them is always dominated by the size of the text. This tool employs a novel algorithm for regular expression matching that leverages the Interleaved Bloom Filter. Regular Expressions are given as input in the command line and support the following operations:

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
Tetrix offers two main commands [index & query] and two utility commands:
```bash
## Index Nucleic Acid DB, with kmer size = 7, 3 hash functions, & 644830 bits per Bloom Filter, each input file represents a bin
tetrix index -k 7 -c 3 -s 644830 -m na -o ibf_idx file1.fna file2.fna
## Query RegEx
tetrix query ibf_idx.ibf "ACGTA(C|G)CC(A|G|T)A"
## Print meta info for an index
tetrix inspect ibf_idx.ibf
## Compute the probability of finding your RegEx
tetrix model -l 1000 ibf_idx.ibf "ACGTA(C|G)CC(A|G|T)A"
```

## Notes
This app was generated from the [SeqAn App Template](https://github.com/seqan/app-template) and makes heavy use of the [SeqAn library](https://github.com/seqan/seqan3/tree/4668203ee1526b4ac3dbdc47869bee72253f684c).
