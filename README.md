# Search for Regular Expressions in Large Datasets
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
3. Configure with cmake ```cmake -DCMAKE_CXX_COMPILER=/path/to/g++-14 ..```
4. Build with make ```make```
5. We recommend you add the ```tetrex``` binary to a location on your PATH ie ```mv ./bin/tetrex /usr/local/bin```

Unforunately, with the latest updates to macOS, configuration with cmake may produce the following error:

```shell
CMake Error at lib/seqan3/build_system/seqan3-config.cmake:106 (message):
  The SDSL library is required, but wasn't found.  Get it from
  https://github.com/xxsds/sdsl-lite
Call Stack (most recent call first):
  lib/seqan3/build_system/seqan3-config.cmake:266 (seqan3_config_error)
  CMakeLists.txt:28 (find_package)
```

As a temporary solution we recommend setting the SDKROOT environment variable via:
```shell
export SDKROOT="$(g++-11 -v 2>&1 | sed -n 's@.*--with-sysroot=\([^ ]*\).*@\1@p')"
```

## Usage
### A Small Example
Tetrix offers two main commands [index & query] and one utility command [inspect]:
```shell
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

### Indexing & Searching the Swissprot Database
A workflow to download, split, index, and query the Swissprot DB from Uniprot (in a Unix-like Environment)
```shell
## Retrieve the Peptide Sequence Library
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gzip -d uniprot_sprot.fasta.gz

## Make a directory to store bins
mkdir sprot_bins && cd sprot_bins

## Split the Swissprot FASTA Library into equal sized bins
gsplit -dn 1024 --additional-suffix .fa ../uniprot_sprot.fasta swissprot_bin_

## Create an HIBF index over the DB
cd ..
tetrex index -k 6 -o sprot_split -i hibf -m aa sprot_bins/*.fa

## Query a domain
tetrex query sprot_split.ibf "LMA(E|Q)GLYN"
/Users/rschwab/Desktop/swissprot_split/swissprot_bin_0346.fa	>sp|Q04896|HME1A_DANRE	LMAQGLYN
/Users/rschwab/Desktop/swissprot_split/swissprot_bin_0346.fa	>sp|P31538|HME1B_XENLA	LMAQGLYN
/Users/rschwab/Desktop/swissprot_split/swissprot_bin_0346.fa	>sp|Q05916|HME1_CHICK	LMAQGLYN
/Users/rschwab/Desktop/swissprot_split/swissprot_bin_0346.fa	>sp|Q05925|HME1_HUMAN	LMAQGLYN
/Users/rschwab/Desktop/swissprot_split/swissprot_bin_0346.fa	>sp|P09065|HME1_MOUSE	LMAQGLYN
/Users/rschwab/Desktop/swissprot_split/swissprot_bin_0346.fa	>sp|P09015|HME2A_DANRE	LMAQGLYN
/Users/rschwab/Desktop/swissprot_split/swissprot_bin_0346.fa	>sp|P52729|HME2A_XENLA	LMAQGLYN
/Users/rschwab/Desktop/swissprot_split/swissprot_bin_0346.fa	>sp|P31533|HME2B_DANRE	LMAQGLYN
/Users/rschwab/Desktop/swissprot_split/swissprot_bin_0346.fa	>sp|P52730|HME2B_XENLA	LMAQGLYN
/Users/rschwab/Desktop/swissprot_split/swissprot_bin_0346.fa	>sp|Q05917|HME2_CHICK	LMAQGLYN
/Users/rschwab/Desktop/swissprot_split/swissprot_bin_0346.fa	>sp|P19622|HME2_HUMAN	LMAQGLYN
/Users/rschwab/Desktop/swissprot_split/swissprot_bin_0346.fa	>sp|P09066|HME2_MOUSE	LMAQGLYN
/Users/rschwab/Desktop/swissprot_split/swissprot_bin_0346.fa	>sp|P09076|HME30_APIME	LMAQGLYN
/Users/rschwab/Desktop/swissprot_split/swissprot_bin_0346.fa	>sp|P09075|HME60_APIME	LMAQGLYN
/Users/rschwab/Desktop/swissprot_split/swissprot_bin_0346.fa	>sp|O02491|HMEN_ANOGA	LMAQGLYN
/Users/rschwab/Desktop/swissprot_split/swissprot_bin_0346.fa	>sp|Q05640|HMEN_ARTSF	LMAQGLYN
/Users/rschwab/Desktop/swissprot_split/swissprot_bin_0346.fa	>sp|P27609|HMEN_BOMMO	LMAQGLYN
/Users/rschwab/Desktop/swissprot_split/swissprot_bin_0346.fa	>sp|P02836|HMEN_DROME	LMAQGLYN
/Users/rschwab/Desktop/swissprot_split/swissprot_bin_0346.fa	>sp|P09145|HMEN_DROVI	LMAQGLYN
/Users/rschwab/Desktop/swissprot_split/swissprot_bin_0346.fa	>sp|P23397|HMEN_HELTR	LMAQGLYN
/Users/rschwab/Desktop/swissprot_split/swissprot_bin_0346.fa	>sp|P09532|HMEN_TRIGR	LMAQGLYN
/Users/rschwab/Desktop/swissprot_split/swissprot_bin_0346.fa	>sp|P27610|HMIN_BOMMO	LMAQGLYN
/Users/rschwab/Desktop/swissprot_split/swissprot_bin_0346.fa	>sp|P05527|HMIN_DROME	LMAQGLYN
/Users/rschwab/Desktop/swissprot_split/swissprot_bin_0811.fa	>sp|Q26601|SMOX2_SCHMA	LMAEGLYN
Query Time: 0.007119
```

Note that, for now, we leave the task of preprocessing the database up to the user. The above example simply splits the database into equal sized bins. You may choose to cluster your sequences or split the bins into variable sizes. You can replace the `gsplit` command to `split` or any program of your choosing. Please note that the database must be split into bins in order to see any significant runtime improvements from the `TetRex` algorithm.

### Using the Prosite Pattern to RegEx Converter
You can use the python executable tetrex_tools.py to convert patterns in the Prosite syntax to POSIX style RegEx's
```shell
./tetrex_tools.py '[LIVMFGAC]-[LIVMTADN]-[LIVFSA]-D-[ST]-G-[STAV]-[STAPDENQ]-{GQ}-[LIVMFSTNC]-{EGK}-[LIVMFGTA]'
(L|I|V|M|F|G|A|C)(L|I|V|M|T|A|D|N)(L|I|V|F|S|A)D(S|T)G(S|T|A|V)(S|T|A|P|D|E|N|Q)(C|I|S|H|T|K|P|A|M|V|Y|E|W|F|R|N|L|D)(L|I|V|M|F|S|T|N|C)(C|I|M|Q|S|Y|V|W|F|H|R|T|N|L|P|A|D)(L|I|V|M|F|G|T|A)
```

## Notes
This app was generated from the [SeqAn App Template](https://github.com/seqan/app-template) and makes heavy use of the [SeqAn library](https://github.com/seqan/seqan3/tree/4668203ee1526b4ac3dbdc47869bee72253f684c).
