#include "construct_nfa.h"


void default_procedure(nfa_t &nfa, nfa_stack_t &stack, lmap_t &node_map, const int &symbol)
{
    node_t node = nfa.addNode();
    node_map[node] = symbol;
    Subgraph twins = {node, node, 0u, 1u, {1u}, Default}; // Start and End nodes are the same
    stack.push(twins);
}


void detect_bad_graphs(const Subgraph& sg1, const Subgraph& sg2, const Subgraph& sg_new, const nfa_t& nfa, const arc_t& carc, catsites_t& cats)
{
    if(sg2.paths >= 15)
    {
        Catsite catsite{sg1.end, sg2.start, sg2.end, carc};
        catsite.addIDs(nfa); // These need to be added now so they can be merged first
        catsite.gaps_ = sg2.lengths;
        cats.push_back(catsite);
        // subgraph.paths = (subgraph.paths/subgraph2.paths)*subgraph2.lengths.size();
    }
    else if(sg_new.paths >= 690000u && sg2.start != sg2.end)
    {
        Catsite catsite{sg1.end, sg2.start, sg2.end, carc};
        catsite.addIDs(nfa); // These need to be added now so they can be merged first
        catsite.gaps_ = sg2.lengths;
        cats.push_back(catsite);
        // subgraph.paths = (subgraph.paths/subgraph2.paths)*subgraph2.lengths.size();
    }
    // subgraph.dumpInfo();
}


void concat_procedure(nfa_t &nfa, lmap_t &node_map, nfa_stack_t &stack, amap_t &arc_map, catsites_t& cats)
{
    Subgraph subgraph2 = stack.top();
    stack.pop();
    Subgraph subgraph1 = stack.top();
    stack.pop();
    arc_t cat_arc = update_arc_map(nfa, node_map, arc_map, subgraph1.end, subgraph2.start);
    Subgraph subgraph{subgraph1.start, subgraph2.end};
    subgraph.concatInfo(subgraph1, subgraph2);
    detect_bad_graphs(subgraph1, subgraph2, subgraph, nfa, cat_arc, cats);
    stack.push(subgraph);
}


void union_procedure(nfa_t &nfa, nfa_stack_t &stack, lmap_t &node_map, amap_t &arc_map)
{
    Subgraph subgraph2 = stack.top();
    stack.pop();
    Subgraph subgraph1 = stack.top();
    stack.pop();

    node_t split_node = nfa.addNode();
    node_map[split_node] = Split;
    update_arc_map(nfa, node_map, arc_map, split_node, subgraph1.start);
    update_arc_map(nfa, node_map, arc_map, split_node, subgraph2.start);

    node_t ghost_node = nfa.addNode();
    node_map[ghost_node] = Ghost;
    update_arc_map(nfa, node_map, arc_map, subgraph1.end, ghost_node);
    update_arc_map(nfa, node_map, arc_map, subgraph2.end, ghost_node);

    Subgraph new_subgraph{split_node, ghost_node};
    new_subgraph.unionInfo(subgraph1, subgraph2);
    stack.push(new_subgraph);
}


void optional_procedure(nfa_t &nfa, nfa_stack_t &stack, lmap_t &node_map, amap_t &arc_map)
{
    Subgraph subgraph = stack.top();
    stack.pop();
    
    node_t split_node = nfa.addNode();
    node_map[split_node] = Split;
    // update_arc_map(nfa, node_map, arc_map, split_node, node.first);
    update_arc_map(nfa, node_map, arc_map, split_node, subgraph.start);

    node_t ghost_node = nfa.addNode();
    node_map[ghost_node] = Ghost;
    update_arc_map(nfa, node_map, arc_map, split_node, ghost_node);
    // update_arc_map(nfa, node_map, arc_map, node.second, ghost_node);
    update_arc_map(nfa, node_map, arc_map, subgraph.end, ghost_node);
    
    Subgraph new_subgraph{split_node, ghost_node};
    new_subgraph.optionInfo(subgraph);
    // new_subgraph.dumpInfo();
    stack.push(new_subgraph);
}


void kleene_procedure(nfa_t &nfa, nfa_stack_t &stack, lmap_t &node_map, const uint8_t &k, amap_t &arc_map)
{
    Subgraph subgraph = stack.top();
    stack.pop();
    size_t subgraph_node_count = subgraph.split_run_count;

    // Create the Split Node which will be the entrance to the subgraph
    node_t split_node = nfa.addNode();
    node_map[split_node] = Split;
    // update_arc_map(nfa, node_map, arc_map, split_node, subgraph.first);
    update_arc_map(nfa, node_map, arc_map, split_node, subgraph.start);

    // Create the Ghost node which will be the exit of the subgraph
    node_t ghost_node = nfa.addNode();
    node_map[ghost_node] = Ghost;
    update_arc_map(nfa, node_map, arc_map, split_node, ghost_node); // * allows to skip the operand completely

    // node_t back_node = subgraph.second;
    node_t back_node = subgraph.end;
    for(uint8_t i = 1; i < (k-1); ++i) // I iterate starting at 1 to represent how I already linearized one cycle
    {
        node_t inner_split = nfa.addNode();
        ++subgraph_node_count;
        node_map[inner_split] = Split;
        update_arc_map(nfa, node_map, arc_map, inner_split, ghost_node);

        Subgraph new_subgraph;
        copy_subgraph(subgraph, nfa, node_map, new_subgraph, arc_map);
        // Copy the subgraph before connecting it to anything downstream
        // so that the DFS in the copy routine doesn't go beyond the subgraph
        update_arc_map(nfa, node_map, arc_map, back_node, inner_split);
        // update_arc_map(nfa, node_map, arc_map, inner_split, new_subgraph.first);
        update_arc_map(nfa, node_map, arc_map, inner_split, new_subgraph.start);

        if(i == (k-2)) // Is this right?
        {
            // update_arc_map(nfa, node_map, arc_map, new_subgraph.second, ghost_node);
            update_arc_map(nfa, node_map, arc_map, new_subgraph.end, ghost_node);
            break;
        }
        // back_node = new_subgraph.second;
        back_node = new_subgraph.end;
    }
    Subgraph stack_subgraph{split_node, ghost_node};
    stack_subgraph.kleeneInfo(subgraph, k);
    stack.push(stack_subgraph);
}


void plus_procedure(nfa_t &nfa, nfa_stack_t &stack, lmap_t &node_map, const uint8_t &k, amap_t &arc_map)
{
    Subgraph subgraph = stack.top();
    stack.pop();
    size_t subgraph_node_count = subgraph.split_run_count;

    node_t ghost_node = nfa.addNode();
    node_map[ghost_node] = Ghost;

    node_t back_node = subgraph.end;
    for(uint8_t i = 1; i < (k-1); ++i) // I iterate starting at 1 because one subgraph already exists in the graph
    {
        Subgraph new_subgraph;
        node_t inner_split = nfa.addNode();
        node_map[inner_split] = Split;
        copy_subgraph(subgraph, nfa, node_map, new_subgraph, arc_map);
        update_arc_map(nfa, node_map, arc_map, back_node, inner_split);
        update_arc_map(nfa, node_map, arc_map, inner_split, ghost_node);
        update_arc_map(nfa, node_map, arc_map, inner_split, new_subgraph.start);
        if(i == (k-2)) // Is this right?
        {
            update_arc_map(nfa, node_map, arc_map, new_subgraph.end, ghost_node);
            break;
        }
        back_node = new_subgraph.end;
    }
    Subgraph stack_subgraph{subgraph.start, ghost_node};
    stack.push(stack_subgraph);
}


bool quant_procedure(nfa_t &nfa, nfa_stack_t &stack, lmap_t &node_map, const uint8_t &k, amap_t &arc_map, const size_t min, const size_t max, catsites_t& cats)
{
    bool skip = false;
    if(min == 0) // Quants like {0,4}
    {
        kleene_procedure(nfa, stack, node_map, (max+1), arc_map);
        if(stack.size() != 1)
        {
            concat_procedure(nfa, node_map, stack, arc_map, cats);
            skip = true;
        } 
        return skip;
    }
    Subgraph subgraph = stack.top();
    if(stack.size() != 1) // If there's one item on the stack then the subgraph being copied is the first subgraph in the regex
    {
        concat_procedure(nfa, node_map, stack, arc_map, cats);
        skip = true; // If it's not first then skip the next concat operator
    }
    size_t extra = (max == 0) ? 0 : (max-min);
    for(size_t i = 1; i < min; ++i)
    {
        Subgraph new_subgraph;
        copy_subgraph(subgraph, nfa, node_map, new_subgraph, arc_map);
        stack.push(new_subgraph);
        concat_procedure(nfa, node_map, stack, arc_map, cats);
    }
    for(size_t i = 0; i < extra; ++i)
    {
        Subgraph new_subgraph;
        copy_subgraph(subgraph, nfa, node_map, new_subgraph, arc_map);
        stack.push(new_subgraph);
        optional_procedure(nfa, stack, node_map, arc_map);
        concat_procedure(nfa, node_map, stack, arc_map, cats);
    }
    return skip;
}


catsites_t construct_kgraph(const std::string &postfix, nfa_t &nfa, lmap_t &node_map, amap_t &arc_map, const uint8_t &k, const bool& verbose)
{
    nfa_stack_t stack;
    node_t start_node = nfa.addNode(); // I don't know why, but a buffer node is necessary for top sort...
    node_map[start_node] = Ghost;
    std::pair<size_t, size_t> min_max;
    catsites_t catsites;
    bool skip = false;
    for(size_t i = 0; i < postfix.size(); i++)
    {
        int symbol = postfix[i];
        if(std::isdigit(symbol)) continue;
        switch(symbol)
        {
            default: // Character
                default_procedure(nfa, stack, node_map, symbol);
                break;
            case '-': // Concat
                if(skip)
                {
                    skip = false;
                    continue;
                }
                concat_procedure(nfa, node_map, stack, arc_map, catsites);
                break;
            case '|': // Or
                union_procedure(nfa, stack, node_map, arc_map);
                break;
            case '?': // Optional
                optional_procedure(nfa, stack, node_map, arc_map);
                break;
            case '*': // Kleene
                kleene_procedure(nfa, stack, node_map, k, arc_map);
                break;
            case '+': // One or More
                plus_procedure(nfa, stack, node_map, k, arc_map);
                break;
            case '{': // Start of Quantifier
                min_max = parse_quant(postfix, i);
                if(min_max == OPT_QUANT) // A special case
                {
                    optional_procedure(nfa, stack, node_map, arc_map);
                    break;
                }
                skip = quant_procedure(nfa, stack, node_map, k, arc_map, min_max.first, min_max.second, catsites);
                break;
            case '}': // End of Quantifier
            case ',':
                break;
        }
    }
    // Cap the graph at both ends
    // With a start node
    // node_t not_head = stack.top().first;
    node_t not_head = stack.top().start;
    update_arc_map(nfa, node_map, arc_map, start_node, not_head);
    
    // and with an accepting node
    node_t match_node = nfa.addNode();
    node_map[match_node] = Match;
    // node_t tail_node = stack.top().second;
    node_t tail_node = stack.top().end;
    update_arc_map(nfa, node_map, arc_map, tail_node, match_node);
    // if(verbose) stack.top().dumpInfo();
    stack.pop();

    // Last but not least...
    // stack.pop();
    assert(stack.size() == 0);
    return catsites;
}
