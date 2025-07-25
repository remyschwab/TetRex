#include "construct_nfa.h"


void default_procedure(nfa_t &nfa, nfa_stack_t &stack, lmap_t &node_map, const int &symbol)
{
    node_t node = nfa.addNode();
    node_map[node] = symbol;
    Subgraph twins = {node, node, 0u, 0u, Default}; // Start and End nodes are the same
    stack.push(twins);
}


void concat_procedure(nfa_t &nfa, lmap_t &node_map, nfa_stack_t &stack, amap_t &arc_map, catsites_t& cats)
{
    Subgraph subgraph2 = stack.top();
    stack.pop();
    Subgraph subgraph1 = stack.top();
    stack.pop();
    arc_t cat_arc = update_arc_map(nfa, node_map, arc_map, subgraph1.end, subgraph2.start);
    if(subgraph2.split_run_count >= 18)
    {
        Catsite catsite{subgraph1.end, subgraph2.start, subgraph2.end, cat_arc};
        catsite.addIDs(nfa); // These need to be added now so they can be merged first
        cats.push_back(catsite);
    }
    size_t new_split_runcount = std::max(subgraph1.split_run_count, subgraph2.split_run_count);
    Subgraph subgraph{subgraph1.start, subgraph2.end, new_split_runcount, 0u, Concat};
    stack.push(subgraph);
}


void union_procedure(nfa_t &nfa, nfa_stack_t &stack, lmap_t &node_map, amap_t &arc_map)
{
    Subgraph subgraph2 = stack.top();
    stack.pop();
    Subgraph subgraph1 = stack.top();
    stack.pop();
    size_t split_run_count = subgraph1.split_run_count + subgraph2.split_run_count;
    if(split_run_count == 0) ++split_run_count;

    node_t split_node = nfa.addNode();
    node_map[split_node] = Split;
    update_arc_map(nfa, node_map, arc_map, split_node, subgraph1.start);
    update_arc_map(nfa, node_map, arc_map, split_node, subgraph2.start);

    node_t ghost_node = nfa.addNode();
    node_map[ghost_node] = Ghost;
    update_arc_map(nfa, node_map, arc_map, subgraph1.end, ghost_node);
    update_arc_map(nfa, node_map, arc_map, subgraph2.end, ghost_node);

    if(node_map[subgraph1.start] == Split) ++split_run_count; // For the Split
    // DBG(split_run_count);
    Subgraph new_subgraph{split_node, ghost_node, split_run_count, 0u, Union};
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
    
    size_t new_split_run_count = node_map[subgraph.start] == Split ? (subgraph.split_run_count+1) : subgraph.split_run_count;
    Subgraph new_subgraph{split_node, ghost_node, new_split_run_count, 0u, Optional};
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
        subgraph_node_count += copy_subgraph(subgraph, nfa, node_map, new_subgraph, arc_map);
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
    Subgraph stack_subgraph{split_node, ghost_node, subgraph_node_count, 0u, Kleene};
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
        ++subgraph_node_count;
        node_map[inner_split] = Split;
        subgraph_node_count += copy_subgraph(subgraph, nfa, node_map, new_subgraph, arc_map);
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
    Subgraph stack_subgraph{subgraph.start, ghost_node, subgraph_node_count, 0u, Plus};
    stack.push(stack_subgraph);
}


std::pair<size_t, size_t> parse_quant(const std::string& postfix, size_t quant_start)
{ // Parse quantifiers that look like {3,4} or {4}
    std::pair<size_t, size_t> min_max(0,0);
    size_t comma = postfix.find(',', quant_start);
    size_t end = postfix.find('}', quant_start);
    if(comma == std::string::npos) // No comma means no max
    {
        min_max.first = std::stoi(postfix.substr(quant_start+1, end - quant_start));
        return min_max;
    }
    min_max.first = std::stoi(postfix.substr(quant_start+1, comma - quant_start));
    min_max.second = std::stoi(postfix.substr(comma + 1, end - comma - 1));
    return min_max;
}


void quant_procedure(nfa_t &nfa, nfa_stack_t &stack, lmap_t &node_map, const uint8_t &k, amap_t &arc_map, const size_t min, const size_t max, catsites_t& cats)
{
    Subgraph subgraph = stack.top();
    size_t extra = (max == 0) ? 0 : (max-min);
    if(min == 0)
    {
        optional_procedure(nfa, stack, node_map, arc_map);
        --extra;
    }
    for(size_t i = 1; i < min; ++i)
    {
        Subgraph new_subgraph;
        copy_subgraph(subgraph, nfa, node_map, new_subgraph, arc_map);
        concat_procedure(nfa, node_map, stack, arc_map, cats);
        stack.push(new_subgraph);
    }
    for(size_t i = 0; i < extra; ++i)
    {
        Subgraph new_subgraph;
        copy_subgraph(subgraph, nfa, node_map, new_subgraph, arc_map);
        stack.push(new_subgraph);
        optional_procedure(nfa, stack, node_map, arc_map);
        concat_procedure(nfa, node_map, stack, arc_map, cats);
    }
}


catsites_t construct_kgraph(const std::string &postfix, nfa_t &nfa, lmap_t &node_map, amap_t &arc_map, const uint8_t &k)
{
    nfa_stack_t stack;
    node_t start_node = nfa.addNode(); // I don't know why, but a buffer node is necessary for top sort...
    node_map[start_node] = Ghost;
    std::pair<size_t, size_t> min_max;
    catsites_t catsites;
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
                quant_procedure(nfa, stack, node_map, k, arc_map, min_max.first, min_max.second, catsites);
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
    stack.pop();
    
    // Last but not least...
    stack.pop();
    return catsites;
}
