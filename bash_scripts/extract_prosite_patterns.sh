#!/usr/bin/env bash

# Process the prosite.dat file

# Extract the patterns from Prosite.dat file (should be 1311)
# Lines containing a prosite pattern start with PA followed by some whitespace
# Some patterns are multiline but they always end with a '.'
grep -E "^PA" $1 | tr -s ' ' | cut -f2 -d' ' | sed -z 's/-\n/-/g' | sed -z 's/\.//g'>patterns.tmp

# Get the Name ID thing
grep "PATTERN" $1 | tr -s ' ' | cut -f2 -d' ' | sed 's/;//g'>names.tmp

# Get the Accession ID
grep "PATTERN" --no-group-separator -A1 $1 | grep -v "PATTERN" | tr -s ' ' | cut -f2 -d' ' | sed -z 's/;//g'>accessions.tmp

# Make a CSV file from the columns
paste names.tmp accessions.tmp patterns.tmp>prosite_patterns.tsv
rm *.tmp # Cleanup
