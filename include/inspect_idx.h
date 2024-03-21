/*
* This module is part of the TetRex search engine for regular expressions
*/

#include <omp.h>
#include "index_base.h"
#include "arg_parse.h"


void inspect_aa_hibf(inspection_arguments const &cmd_args);
void inspect_aa_ibf(inspection_arguments const &cmd_args);
void inspect_dna_hibf(inspection_arguments const &cmd_args);
void inspect_dna_ibf(inspection_arguments const &cmd_args);
void drive_inspection(const inspection_arguments &args);
