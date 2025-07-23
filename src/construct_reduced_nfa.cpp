#include "construct_nfa.h"


bool redundancy_test(buffer_t &buffer)
{
    int sym2 = buffer.top();
    buffer.pop();
    int sym1 = buffer.top();
    buffer.push(sym2);
    bool redundant = (sym1 == sym2);
    return redundant;
}


bool twin_test(Subgraph &node_pair)
{
    return node_pair.start == node_pair.end;
}


Subgraph add_node(nfa_t &nfa, lmap_t &node_map, const int &symbol)
{
    node_t node = nfa.addNode();
    node_map[node] = symbol;
    Subgraph twins{node, node}; // Source and target are the same node
    return twins;
}


void twin_procedure(Subgraph &node_pair, buffer_t &buffer, nfa_t &nfa, nfa_stack_t &stack, lmap_t &node_map)
{
    int sym = buffer.top();
    buffer.pop();
    node_pair = add_node(nfa, node_map, sym);
}


void default_procedure(buffer_t &buffer, const int symbol, nfa_stack_t &stack)
{
    buffer.push(symbol);
    node_t node;
    Subgraph twins = {node, node}; // Source and target are the same node
    stack.push(twins); // Basically just push a pair of dummy twins onto the stack
}


void concat_procedure(nfa_t &nfa, lmap_t &node_map, nfa_stack_t &stack, amap_t &arc_map, buffer_t &buffer, catsites_t& cats)
{
    Subgraph subgraph2 = stack.top();
    stack.pop();
    Subgraph subgraph1 = stack.top();
    stack.pop();
    if(twin_test(subgraph1)) twin_procedure(subgraph1, buffer, nfa, stack, node_map);
    if(twin_test(subgraph2)) twin_procedure(subgraph2, buffer, nfa, stack, node_map);
    arc_t cat_arc = update_arc_map(nfa, node_map, arc_map, subgraph1.end, subgraph2.start);
    if(subgraph2.split_run_count >= 20)
    {
        Catsite catsite{subgraph1.end, subgraph2.start, subgraph2.end, cat_arc};
        catsite.addIDs(nfa); // These need to be added now so they can be merged first
        cats.push_back(catsite);
    }
    Subgraph node_pair{subgraph1.start, subgraph2.end};
    stack.push(node_pair);
}


void union_procedure(nfa_t &nfa, nfa_stack_t &stack, lmap_t &node_map, amap_t &arc_map, buffer_t &buffer)
{
    Subgraph subgraph2 = stack.top();
    stack.pop();
    Subgraph subgraph1 = stack.top();
    stack.pop();
    size_t split_run_count = subgraph1.split_run_count + subgraph2.split_run_count;
    if(split_run_count == 0) ++split_run_count;

    bool redundant = redundancy_test(buffer);
    if(twin_test(subgraph2) && twin_test(subgraph1) && redundant) // Check for redundant twins [in reduced space]
    {
        int sym = buffer.top();
        buffer.pop();
        buffer.pop();
        default_procedure(buffer, sym, stack);
        return;
    }
    if(twin_test(subgraph1)) twin_procedure(subgraph1, buffer, nfa, stack, node_map);
    if(twin_test(subgraph2)) twin_procedure(subgraph2, buffer, nfa, stack, node_map);

    node_t split_node = nfa.addNode();
    node_map[split_node] = Split;
    update_arc_map(nfa, node_map, arc_map, split_node, subgraph1.start);
    update_arc_map(nfa, node_map, arc_map, split_node, subgraph2.start);

    node_t ghost_node = nfa.addNode();
    node_map[ghost_node] = Ghost;
    update_arc_map(nfa, node_map, arc_map, subgraph1.end, ghost_node);
    update_arc_map(nfa, node_map, arc_map, subgraph2.end, ghost_node);

    Subgraph node_pair{split_node, ghost_node, split_run_count, 0u, Union};
    stack.push(node_pair);
}


void optional_procedure(nfa_t &nfa, nfa_stack_t &stack, lmap_t &node_map, amap_t &arc_map, buffer_t &buffer)
{
    Subgraph node = stack.top();
    stack.pop();
    if(twin_test(node)) twin_procedure(node, buffer, nfa, stack, node_map);

    node_t split_node = nfa.addNode();
    node_map[split_node] = Split;
    update_arc_map(nfa, node_map, arc_map, split_node, node.start);

    node_t ghost_node = nfa.addNode();
    node_map[ghost_node] = Ghost;
    update_arc_map(nfa, node_map, arc_map, split_node, ghost_node);
    update_arc_map(nfa, node_map, arc_map, node.end, ghost_node);

    Subgraph node_pair{split_node, ghost_node};
    stack.push(node_pair);
}


void kleene_procedure(nfa_t &nfa, nfa_stack_t &stack, lmap_t &node_map, const uint8_t &k, amap_t &arc_map, buffer_t &buffer)
{
    Subgraph subgraph = stack.top();
    stack.pop();
    if(twin_test(subgraph)) twin_procedure(subgraph, buffer, nfa, stack, node_map);

    // Create the Split Node which will be the entrance to the subgraph
    node_t split_node = nfa.addNode();
    node_map[split_node] = Split;
    update_arc_map(nfa, node_map, arc_map, split_node, subgraph.start);

    // Create the Ghost node which will be the exit of the subgraph
    node_t ghost_node = nfa.addNode();
    node_map[ghost_node] = Ghost;
    update_arc_map(nfa, node_map, arc_map, split_node, ghost_node); // * allows to skip the operand completely


    node_t back_node = subgraph.end;
    for(uint8_t i = 1; i < (k-1); ++i) // I iterate starting at 1 to represent how I already linearized one cycle
    {
        node_t inner_split = nfa.addNode();
        node_map[inner_split] = Split;
        update_arc_map(nfa, node_map, arc_map, inner_split, ghost_node);

        Subgraph new_subgraph;
        copy_subgraph(subgraph, nfa, node_map, new_subgraph, arc_map);
        // Copy the subgraph before connecting it to anything downstream
        // so that the DFS in the copy routine doesn't go beyond the subgraph
        update_arc_map(nfa, node_map, arc_map, back_node, inner_split);
        update_arc_map(nfa, node_map, arc_map, inner_split, new_subgraph.start);

        if(i == (k-2)) // Is this right?
        {
            update_arc_map(nfa, node_map, arc_map, new_subgraph.end, ghost_node);
            break;
        }
        back_node = new_subgraph.end;
    }
    Subgraph node_pair{split_node, ghost_node};
    stack.push(node_pair);
}


void plus_procedure(nfa_t &nfa, nfa_stack_t &stack, lmap_t &node_map, const uint8_t &k, amap_t &arc_map, buffer_t &buffer)
{
    Subgraph subgraph = stack.top();
    stack.pop();
    if(twin_test(subgraph)) twin_procedure(subgraph, buffer, nfa, stack, node_map);

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
    Subgraph node_pair{subgraph.start, ghost_node};
    stack.push(node_pair);
}


catsites_t construct_reduced_kgraph(const std::string &postfix, nfa_t &nfa, lmap_t &node_map, amap_t &arc_map, const uint8_t &k)
{
    nfa_stack_t stack;
    buffer_t buffer;
    node_t start_node = nfa.addNode(); // I don't know why, but a buffer node is necessary for top sort...
    node_map[start_node] = Ghost;
    catsites_t catsites;
    for(size_t i = 0; i < postfix.size(); i++)
    {
        int symbol = postfix[i];
        switch(symbol)
        {
            default: // Character
                // default_procedure(nfa, stack, node_map, symbol);
                default_procedure(buffer, symbol, stack);
                break;
            case '-': // Concat
                concat_procedure(nfa, node_map, stack, arc_map, buffer, catsites);
                break;
            case '|': // Or
                union_procedure(nfa, stack, node_map, arc_map, buffer);
                break;
            case '?': // Optional
                optional_procedure(nfa, stack, node_map, arc_map, buffer);
                break;
            case '*': // Kleene
                kleene_procedure(nfa, stack, node_map, k, arc_map, buffer);
                break;
            case '+': // One or More
                plus_procedure(nfa, stack, node_map, k, arc_map, buffer);
                break;
        }
    }
    // Cap the graph at both ends
    // With a start node
    node_t not_head = stack.top().start;
    update_arc_map(nfa, node_map, arc_map, start_node, not_head);
    
    // and with an accepting node
    node_t match_node = nfa.addNode();
    node_map[match_node] = Match;
    node_t tail_node = stack.top().end;
    update_arc_map(nfa, node_map, arc_map, tail_node, match_node);
    stack.pop();
    
    // Last but not least...
    stack.pop();
    return catsites;
}
