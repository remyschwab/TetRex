#include "construct_nfa.h"


void export_nfa_img(nfa_t &nfa)
{
    lemon::graphToEps(nfa, "nfa.eps").run();
}


node_pair_t copy_subgraph(node_pair_t &node_pair)
{
    //TODO: return a copy of whatever is between a source and target node
    // This is necessary for making an acyclic NFA with * and +
}


void default_procedure(nfa_t &nfa, nfa_stack_t &stack, lmap_t &node_map, const int &symbol)
{
    node_t node = nfa.addNode();
    node_map[node] = symbol;
    node_pair_t twins = std::make_pair(node, node);
    stack.push(twins);
}


void concat_procedure(nfa_t &nfa, nfa_stack_t &stack)
{
    node_pair_t node2 = stack.top();
    stack.pop();
    node_pair_t node1 = stack.top();
    stack.pop();
    nfa.addArc(node1.second, node2.first);
    node_pair_t node_pair = std::make_pair(node1.second, node2.first);
    stack.push(node_pair);
}


void union_procedure(nfa_t &nfa, nfa_stack_t &stack, lmap_t &node_map)
{
    node_pair_t node2 = stack.top();
    stack.pop();
    node_pair_t node1 = stack.top();
    stack.pop();

    node_t split_node = nfa.addNode();
    node_map[split_node] = SplitU;
    nfa.addArc(split_node, node1.first);
    nfa.addArc(split_node, node2.first);

    node_t ghost_node = nfa.addNode();
    node_map[ghost_node] = Ghost;
    nfa.addArc(node1.second, ghost_node);
    nfa.addArc(node2.second, ghost_node);

    node_pair_t node_pair = std::make_pair(split_node, ghost_node);
    stack.push(node_pair);
}


void optional_procedure(nfa_t &nfa, nfa_stack_t &stack, lmap_t &node_map)
{
    node_pair_t node = stack.top();
    stack.pop();

    node_t split_node = nfa.addNode();
    node_map[split_node] = SplitU;
    nfa.addArc(split_node, node.first);

    node_t ghost_node = nfa.addNode();
    node_map[ghost_node] = Ghost;
    nfa.addArc(split_node, ghost_node);
    nfa.addArc(node.second, ghost_node);

    node_pair_t node_pair = std::make_pair(split_node, ghost_node);
    stack.push(node_pair);
}


void kleene_procedure(nfa_t &nfa, nfa_stack_t &stack, lmap_t &node_map, const uint8_t &k)
{
    node_pair_t node = stack.top();
    stack.pop();
    const int symbol = node_map[node.first]; // This isn't a full solution!

    node_t split_node = nfa.addNode();
    node_map[split_node] = SplitK;

    node_t ghost_node = nfa.addNode();
    node_map[ghost_node] = Ghost;
    nfa.addArc(split_node, ghost_node);

    node_t back_node = split_node; // This functions as a pointer
    for(uint8_t i = 1; i < (k-1); ++i) // I iterate starting at 1 to represent how I already linearized one cycle
    {
        node_t new_node = nfa.addNode();
        node_map[new_node] = symbol;
        nfa.addArc(back_node, new_node);
        nfa.addArc(new_node, ghost_node);
        back_node = new_node;
    }
    node_pair_t node_pair = std::make_pair(split_node, ghost_node);
    stack.push(node_pair);
}


void plus_procedure(nfa_t &nfa, nfa_stack_t &stack, lmap_t &node_map, const uint8_t &k)
{
    node_pair_t node = stack.top();
    stack.pop();
    const int symbol = node_map[node.first]; // This isn't a full solution!

    node_t split_node = nfa.addNode();
    node_map[split_node] = SplitP;
    nfa.addArc(split_node, node.first);

    node_t ghost_node = nfa.addNode();
    node_map[ghost_node] = Ghost;
    nfa.addArc(node.second, ghost_node);

    node_t back_node = node.second;
    for(uint8_t i = 1; i < (k-1); ++i) // I iterate starting at 1 to represent how I already linearized one cycle
    {
        node_t new_node = nfa.addNode();
        node_map[new_node] = symbol;
        nfa.addArc(back_node, new_node);
        nfa.addArc(new_node, ghost_node);
        back_node = new_node;
    }
    node_pair_t node_pair = std::make_pair(split_node, ghost_node);
    stack.push(node_pair);
}


void construct_graph(const std::string &postfix, nfa_t &nfa, lmap_t &node_map, const uint8_t &k)
{
    nfa_stack_t stack;
    for(size_t i = 0; i < postfix.size(); i++)
    {
        int symbol = postfix[i];
        switch(symbol)
        {
            default: // Character
                default_procedure(nfa, stack, node_map, symbol);
                break;
            case '-': // Concat
                concat_procedure(nfa, stack);
                break;
            case '|': // Or
                union_procedure(nfa, stack, node_map);
                break;
            case '?': // Optional
                optional_procedure(nfa, stack, node_map);
                break;
            case '*': // Kleene
                kleene_procedure(nfa, stack, node_map, k);
                break;
            case '+': // One or More
                plus_procedure(nfa, stack, node_map, k);
                break;
        }
    }
    node_t match_node = nfa.addNode();
    node_map[match_node] = 256;
    node_t tail_node = stack.top().second;
    nfa.addArc(tail_node, match_node);
    stack.pop();
}