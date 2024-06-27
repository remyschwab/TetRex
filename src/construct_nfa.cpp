#include "construct_nfa.h"


// void export_nfa_img(nfa_t &nfa, std::string &title)
// {
//     lemon::graphToEps(nfa, "/Users/rschwab/Desktop/nfa.eps").title(title).copyright("(C) 2003-2009 LEMON Project").run();
// }

void paste_to_graph(nfa_t &NFA, lmap_t &node_map, const node_t &reference_node, node_t &paste_node)
{
    int symbol = node_map[reference_node];
    paste_node = NFA.addNode();
    node_map[paste_node] = symbol;
}

void copy_subgraph(node_pair_t &subgraph, nfa_t &NFA, lmap_t &node_map, node_pair_t &subgraph_copy, amap_t &arc_map)
{
    node_t new_node;
    paste_to_graph(NFA, node_map, subgraph.first, new_node);
    subgraph_copy.first = new_node;
    // If the operand is just a single character then just copy that one node
    if(subgraph.first == subgraph.second)
    {
        subgraph_copy.second = new_node;
        return;
    }
    // If the operand is a more complicated subgraph, then traverse with BFS copying nodes and arcs along the way
    nfa_t::NodeMap<node_t> reference_to_copy_map(NFA);
    robin_hood::unordered_set<int> copied_targets;
    lemon::Dfs<nfa_t> dfs(NFA);
    dfs.init();
    dfs.addSource(subgraph.first);
    reference_to_copy_map[subgraph.first] = new_node;
    while (!dfs.emptyQueue())
    {
        arc_t arc = dfs.processNextArc();
        node_t source = NFA.source(arc);
        node_t target = NFA.target(arc);
        node_t source_copy = reference_to_copy_map[source];
        if(copied_targets.find(NFA.id(target)) != copied_targets.end())
        {
            new_node = reference_to_copy_map[target];
            update_arc_map(NFA, node_map, arc_map, source_copy, new_node);
            continue;
        }
        paste_to_graph(NFA, node_map, target, new_node);
        update_arc_map(NFA, node_map, arc_map, source_copy, new_node);
        copied_targets.insert(NFA.id(target));
        reference_to_copy_map[target] = new_node;
    }
    subgraph_copy.second = new_node;
}

std::string generate_kmer_seq(uint64_t &kmer, uint8_t &k)
{
    robin_hood::unordered_map<uint64_t, char> nucleotides = {
        {0, 'A'},
        {1, 'C'},
        {2, 'T'},
        {3, 'G'}
    };
    std::string kmer_seq = "";
    size_t countdown = k;
    while(countdown != 0)
    {
        kmer_seq += nucleotides[(kmer & 0b11)];
        kmer = (kmer>>2);
        --countdown;
    }
    std::reverse(kmer_seq.begin(), kmer_seq.end());
    return kmer_seq;
}


void print_node_pointers(const amap_t &arc_map, nfa_t &nfa)
{
    for(auto [id, target_pair]: arc_map)
        std::cout << id << " --> " << nfa.id(target_pair.first) << "-" << nfa.id(target_pair.second) << std::endl;
}


void print_kgraph_arcs(const nfa_t &NFA)
{
    for(auto &&arc: NFA.arcs())
        std::cout << NFA.id(NFA.source(arc)) << " --> " <<NFA.id(NFA.target(arc)) << std::endl;
}


void print_node_ids(nfa_t &NFA, lmap_t &nmap)
{
    for(auto &&node: NFA.nodes())
    {
        int symbol = nmap[node];
        if(symbol < 256)
        {
            std::cout << NFA.id(node) << " " << static_cast<char>(nmap[node]) << std::endl;
        }
        else if(symbol == 256)
        {
            std::cout << NFA.id(node) << " ∆" << std::endl;
        }
        else if(symbol == 257)
        {
            std::cout << NFA.id(node) << " •" << std::endl;
        }
        else if(symbol == 258)
        {
            std::cout << NFA.id(node) << " Ø" << std::endl;
        }
    }
}


void update_arc_map(nfa_t &NFA, lmap_t &node_map, amap_t &arc_map, node_t &source, node_t &target)
{
    NFA.addArc(source, target); // Update the internal NFA arc map I guess
    int source_id = NFA.id(source);
    int symbol = node_map[source];
    if(symbol < 258) // If the node is not a split just place the same target in both slots
    {
        arc_map[source_id].first = target;
        arc_map[source_id].second = target; // I don't know if this is the best approach but fine for now
    }
    else // If it's a split
    {
        if(arc_map.find(source_id) == arc_map.end())
        {
            arc_map[source_id].first = target; // Place it in the first slot if we haven't seen it before
        }
        else
        {
            arc_map[source_id].second = target; // Place it in the second if we have
        }
    }
}


std::vector<int> run_top_sort(nfa_t &NFA)
{
    wmap_t list(NFA);
    lemon::topologicalSort(NFA, list);
    std::vector<int> priority_map;
    priority_map.resize(NFA.nodeNum());
    size_t rank = 1; // Start node is always at 0
    for(auto &&it: list.order)
    {
        priority_map[NFA.id(it)] = rank;
        rank++;
    }
    return priority_map;
}


void default_procedure(nfa_t &nfa, nfa_stack_t &stack, lmap_t &node_map, const int &symbol)
{
    node_t node = nfa.addNode();
    node_map[node] = symbol;
    node_pair_t twins = std::make_pair(node, node); // Source and target are the same node
    stack.push(twins);
}


void concat_procedure(nfa_t &nfa, lmap_t &node_map, nfa_stack_t &stack, amap_t &arc_map)
{
    node_pair_t node2 = stack.top();
    stack.pop();
    node_pair_t node1 = stack.top();
    stack.pop();
    update_arc_map(nfa, node_map, arc_map, node1.second, node2.first);
    node_pair_t node_pair = std::make_pair(node1.first, node2.second);
    stack.push(node_pair);
}


void union_procedure(nfa_t &nfa, nfa_stack_t &stack, lmap_t &node_map, amap_t &arc_map)
{
    node_pair_t node2 = stack.top();
    stack.pop();
    node_pair_t node1 = stack.top();
    stack.pop();

    node_t split_node = nfa.addNode();
    node_map[split_node] = Split;
    update_arc_map(nfa, node_map, arc_map, split_node, node1.first);
    update_arc_map(nfa, node_map, arc_map, split_node, node2.first);

    node_t ghost_node = nfa.addNode();
    node_map[ghost_node] = Ghost;
    update_arc_map(nfa, node_map, arc_map, node1.second, ghost_node);
    update_arc_map(nfa, node_map, arc_map, node2.second, ghost_node);

    node_pair_t node_pair = std::make_pair(split_node, ghost_node);
    stack.push(node_pair);
}


void optional_procedure(nfa_t &nfa, nfa_stack_t &stack, lmap_t &node_map, amap_t &arc_map)
{
    node_pair_t node = stack.top();
    stack.pop();

    node_t split_node = nfa.addNode();
    node_map[split_node] = Split;
    update_arc_map(nfa, node_map, arc_map, split_node, node.first);

    node_t ghost_node = nfa.addNode();
    node_map[ghost_node] = Ghost;
    update_arc_map(nfa, node_map, arc_map, split_node, ghost_node);
    update_arc_map(nfa, node_map, arc_map, node.second, ghost_node);

    node_pair_t node_pair = std::make_pair(split_node, ghost_node);
    stack.push(node_pair);
}


void kleene_procedure(nfa_t &nfa, nfa_stack_t &stack, lmap_t &node_map, const uint8_t &k, amap_t &arc_map)
{
    node_pair_t subgraph = stack.top();
    stack.pop();

    // Create the Split Node which will be the entrance to the subgraph
    node_t split_node = nfa.addNode();
    node_map[split_node] = Split;
    update_arc_map(nfa, node_map, arc_map, split_node, subgraph.first);

    // Create the Ghost node which will be the exit of the subgraph
    node_t ghost_node = nfa.addNode();
    node_map[ghost_node] = Ghost;
    update_arc_map(nfa, node_map, arc_map, split_node, ghost_node); // * allows to skip the operand completely


    node_t *back_node = &subgraph.second;
    for(uint8_t i = 1; i < (k-1); ++i) // I iterate starting at 1 to represent how I already linearized one cycle
    {
        node_t inner_split = nfa.addNode();
        node_map[inner_split] = Split;
        update_arc_map(nfa, node_map, arc_map, inner_split, ghost_node);

        node_pair_t new_subgraph;
        copy_subgraph(subgraph, nfa, node_map, new_subgraph, arc_map);
        // Copy the subgraph before connecting it to anything downstream
        // so that the DFS in the copy routine doesn't go beyond the subgraph
        update_arc_map(nfa, node_map, arc_map, *back_node, inner_split);
        update_arc_map(nfa, node_map, arc_map, inner_split, new_subgraph.first);

        if(i == (k-2)) // Is this right?
        {
            update_arc_map(nfa, node_map, arc_map, new_subgraph.second, ghost_node);
            break;
        }
        back_node = &new_subgraph.second;
    }
    node_pair_t node_pair = std::make_pair(split_node, ghost_node);
    stack.push(node_pair);
}


void plus_procedure(nfa_t &nfa, nfa_stack_t &stack, lmap_t &node_map, const uint8_t &k, amap_t &arc_map)
{
    node_pair_t node = stack.top();
    stack.pop();
    const int symbol = node_map[node.first]; // This isn't a full solution!

    node_t split_node = nfa.addNode();
    node_map[split_node] = Split;
    update_arc_map(nfa, node_map, arc_map, node.first, split_node);

    node_t ghost_node = nfa.addNode();
    node_map[ghost_node] = Ghost;
    update_arc_map(nfa, node_map, arc_map, split_node, ghost_node);

    node_t *back_node = &split_node;
    for(uint8_t i = 1; i < (k-1); ++i) // I iterate starting at 1 to represent how I already linearized one cycle
    {
        node_t new_node = nfa.addNode();
        node_map[new_node] = symbol;
        update_arc_map(nfa, node_map, arc_map, *back_node, new_node);
        if(i == (k-2))
        {
            update_arc_map(nfa, node_map, arc_map, new_node, ghost_node);
            break;
        }
        node_t inner_split = nfa.addNode();
        node_map[inner_split] = Split;
        update_arc_map(nfa, node_map, arc_map, new_node, inner_split);
        back_node = &inner_split;
    }
    node_pair_t node_pair = std::make_pair(node.first, ghost_node);
    stack.push(node_pair);
}


void construct_kgraph(const std::string &postfix, nfa_t &nfa, lmap_t &node_map, amap_t &arc_map, const uint8_t &k)
{
    nfa_stack_t stack;
    node_t start_node = nfa.addNode(); // I don't know why, but a buffer node is necessary for top sort...
    node_map[start_node] = Ghost;
    for(size_t i = 0; i < postfix.size(); i++)
    {
        int symbol = postfix[i];
        switch(symbol)
        {
            default: // Character
                default_procedure(nfa, stack, node_map, symbol);
                break;
            case '-': // Concat
                concat_procedure(nfa, node_map, stack, arc_map);
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
        }
    }
    // Cap the graph at both ends
    // With a start node
    node_t not_head = stack.top().first;
    update_arc_map(nfa, node_map, arc_map, start_node, not_head);
    
    // and with an accepting node
    node_t match_node = nfa.addNode();
    node_map[match_node] = Match;
    node_t tail_node = stack.top().second;
    update_arc_map(nfa, node_map, arc_map, tail_node, match_node);
    stack.pop();
    
    // Last but not least...
    stack.pop();
}
