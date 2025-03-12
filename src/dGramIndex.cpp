#include "arg_parse.h"
#include "dGramIndex.h"



void drive_dindex(const dindex_arguments &cmd_args)
{
    DGramIndex dindex(cmd_args.min, cmd_args.max, cmd_args.pad, cmd_args.hash_count, cmd_args.fpr, cmd_args.acid_libs);
    dindex.track_bins();
    dindex.spawn_agent();
    // uint64_t const dgram = 687194767365;
    // dindex.query(dgram);
}
