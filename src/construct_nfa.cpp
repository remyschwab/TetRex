#include "construct_nfa.h"


void default_procedure(nfa_t &nfa, nfa_stack_t &stack, lmap_t &node_map, const int &symbol)
{
    node_t node = nfa.addNode();
    node_map[node] = symbol;
    // node_pair_t twins = std::make_pair(node, node); // Source and target are the same node
    node_pair_t twins = std::make_tuple(node, node, 1u);
    stack.push(twins);
}


void concat_procedure(nfa_t &nfa, lmap_t &node_map, nfa_stack_t &stack, amap_t &arc_map, catsites_t& cats)
{
    node_pair_t node2 = stack.top();
    stack.pop();
    node_pair_t node1 = stack.top();
    stack.pop();
    size_t combined_node_count = std::get<2>(node2)+std::get<2>(node1);
    // arc_t cat_arc = update_arc_map(nfa, node_map, arc_map, node1.second, node2.first);
    arc_t cat_arc = update_arc_map(nfa, node_map, arc_map, std::get<1>(node1), std::get<0>(node2));
    if(std::get<2>(node2) > 10) cats.push_back(std::make_tuple(cat_arc, std::get<1>(node2)));
    // node_pair_t node_pair = std::make_tuple(node1.first, node2.second, combined_node_count);
    node_pair_t node_pair = std::make_tuple(std::get<0>(node1), std::get<1>(node2), combined_node_count);
    stack.push(node_pair);
}


void union_procedure(nfa_t &nfa, nfa_stack_t &stack, lmap_t &node_map, amap_t &arc_map)
{
    node_pair_t node2 = stack.top();
    stack.pop();
    node_pair_t node1 = stack.top();
    stack.pop();
    size_t combined_node_count = std::get<2>(node2)+std::get<2>(node1);

    node_t split_node = nfa.addNode();
    node_map[split_node] = Split;
    // update_arc_map(nfa, node_map, arc_map, split_node, node1.first);
    update_arc_map(nfa, node_map, arc_map, split_node, std::get<0>(node1));
    // update_arc_map(nfa, node_map, arc_map, split_node, node2.first);
    update_arc_map(nfa, node_map, arc_map, split_node, std::get<0>(node2));

    node_t ghost_node = nfa.addNode();
    node_map[ghost_node] = Ghost;
    // update_arc_map(nfa, node_map, arc_map, node1.second, ghost_node);
    update_arc_map(nfa, node_map, arc_map, std::get<1>(node1), ghost_node);
    update_arc_map(nfa, node_map, arc_map, std::get<1>(node2), ghost_node);

    combined_node_count += 2; // For the Split and Ghost node
    node_pair_t node_pair = std::make_tuple(split_node, ghost_node, combined_node_count);
    stack.push(node_pair);
}


void optional_procedure(nfa_t &nfa, nfa_stack_t &stack, lmap_t &node_map, amap_t &arc_map)
{
    node_pair_t node = stack.top();
    stack.pop();
    
    node_t split_node = nfa.addNode();
    node_map[split_node] = Split;
    // update_arc_map(nfa, node_map, arc_map, split_node, node.first);
    update_arc_map(nfa, node_map, arc_map, split_node, std::get<0>(node));

    node_t ghost_node = nfa.addNode();
    node_map[ghost_node] = Ghost;
    update_arc_map(nfa, node_map, arc_map, split_node, ghost_node);
    // update_arc_map(nfa, node_map, arc_map, node.second, ghost_node);
    update_arc_map(nfa, node_map, arc_map, std::get<1>(node), ghost_node);

    node_pair_t node_pair = std::make_tuple(split_node, ghost_node, (std::get<2>(node)+2));
    stack.push(node_pair);
}


void kleene_procedure(nfa_t &nfa, nfa_stack_t &stack, lmap_t &node_map, const uint8_t &k, amap_t &arc_map)
{
    node_pair_t subgraph = stack.top();
    stack.pop();
    size_t subgraph_node_count = std::get<2>(subgraph);

    // Create the Split Node which will be the entrance to the subgraph
    node_t split_node = nfa.addNode();
    node_map[split_node] = Split;
    // update_arc_map(nfa, node_map, arc_map, split_node, subgraph.first);
    update_arc_map(nfa, node_map, arc_map, split_node, std::get<0>(subgraph));

    // Create the Ghost node which will be the exit of the subgraph
    node_t ghost_node = nfa.addNode();
    ++subgraph_node_count;
    node_map[ghost_node] = Ghost;
    update_arc_map(nfa, node_map, arc_map, split_node, ghost_node); // * allows to skip the operand completely

    // node_t back_node = subgraph.second;
    node_t back_node = std::get<1>(subgraph);
    for(uint8_t i = 1; i < (k-1); ++i) // I iterate starting at 1 to represent how I already linearized one cycle
    {
        node_t inner_split = nfa.addNode();
        ++subgraph_node_count;
        node_map[inner_split] = Split;
        update_arc_map(nfa, node_map, arc_map, inner_split, ghost_node);

        node_pair_t new_subgraph;
        subgraph_node_count += copy_subgraph(subgraph, nfa, node_map, new_subgraph, arc_map);
        // Copy the subgraph before connecting it to anything downstream
        // so that the DFS in the copy routine doesn't go beyond the subgraph
        update_arc_map(nfa, node_map, arc_map, back_node, inner_split);
        // update_arc_map(nfa, node_map, arc_map, inner_split, new_subgraph.first);
        update_arc_map(nfa, node_map, arc_map, inner_split, std::get<0>(new_subgraph));

        if(i == (k-2)) // Is this right?
        {
            // update_arc_map(nfa, node_map, arc_map, new_subgraph.second, ghost_node);
            update_arc_map(nfa, node_map, arc_map, std::get<1>(new_subgraph), ghost_node);
            break;
        }
        // back_node = new_subgraph.second;
        back_node = std::get<1>(new_subgraph);
    }
    node_pair_t node_pair = std::make_tuple(split_node, ghost_node, subgraph_node_count);
    stack.push(node_pair);
}


void plus_procedure(nfa_t &nfa, nfa_stack_t &stack, lmap_t &node_map, const uint8_t &k, amap_t &arc_map)
{
    node_pair_t subgraph = stack.top();
    stack.pop();
    size_t subgraph_node_count = std::get<2>(subgraph);

    node_t ghost_node = nfa.addNode();
    ++subgraph_node_count;
    node_map[ghost_node] = Ghost;

    // node_t back_node = subgraph.second;
    node_t back_node = std::get<1>(subgraph);
    for(uint8_t i = 1; i < (k-1); ++i) // I iterate starting at 1 because one subgraph already exists in the graph
    {
        node_pair_t new_subgraph;
        node_t inner_split = nfa.addNode();
        ++subgraph_node_count;
        node_map[inner_split] = Split;
        subgraph_node_count += copy_subgraph(subgraph, nfa, node_map, new_subgraph, arc_map);
        update_arc_map(nfa, node_map, arc_map, back_node, inner_split);
        update_arc_map(nfa, node_map, arc_map, inner_split, ghost_node);
        // update_arc_map(nfa, node_map, arc_map, inner_split, new_subgraph.first);
        update_arc_map(nfa, node_map, arc_map, inner_split, std::get<0>(new_subgraph));
        if(i == (k-2)) // Is this right?
        {
            // update_arc_map(nfa, node_map, arc_map, new_subgraph.second, ghost_node);
            update_arc_map(nfa, node_map, arc_map, std::get<1>(new_subgraph), ghost_node);
            break;
        }
        // back_node = new_subgraph.second;
        back_node = std::get<1>(new_subgraph);
    }
    // node_pair_t node_pair = std::make_tuple(subgraph.first, ghost_node, subgraph_node_count);
    node_pair_t node_pair = std::make_tuple(std::get<0>(subgraph), ghost_node, subgraph_node_count);
    stack.push(node_pair);
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
    node_pair_t subgraph = stack.top();
    size_t extra = max-min;
    assert(max != 0); // This should be done with the USER in mind
    if(min == 0)
    {
        optional_procedure(nfa, stack, node_map, arc_map);
        --extra;
    }
    for(size_t i = 1; i < min; ++i)
    {
        node_pair_t new_subgraph;
        copy_subgraph(subgraph, nfa, node_map, new_subgraph, arc_map);
        concat_procedure(nfa, node_map, stack, arc_map, cats);
        stack.push(new_subgraph);
    }
    for(size_t i = 0; i < extra; ++i)
    {
        node_pair_t new_subgraph;
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
    node_t not_head = std::get<0>(stack.top());
    update_arc_map(nfa, node_map, arc_map, start_node, not_head);
    
    // and with an accepting node
    node_t match_node = nfa.addNode();
    node_map[match_node] = Match;
    // node_t tail_node = stack.top().second;
    node_t tail_node = std::get<1>(stack.top());
    update_arc_map(nfa, node_map, arc_map, tail_node, match_node);
    stack.pop();
    
    // Last but not least...
    stack.pop();
    return catsites;
}
