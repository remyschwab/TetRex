#include "construct_reduced_nfa.h"


void copy_subgraph(const Subgraph& subgraph, nfa_t& NFA, lmap_t& node_map, Subgraph& subgraph_copy, amap_t& arc_map, buffer_t& buffer)
{
    using namespace lemon;
    using Digraph = SmartDigraph;

    if(subgraph.start == subgraph.end) // Copying is very simple if the subgraph is just a single node
    {
        node_t new_twin = NFA.addNode();
        node_map[new_twin] = buffer.top();
        subgraph_copy.start = new_twin;
        subgraph_copy.end = new_twin;
        subgraph_copy.copyMeta(subgraph);
        return;
    }

    // --- Keep snapshots of the original nodes/arcs (ids are dense & stable in SmartDigraph)
    std::vector<Digraph::Node> origNodes;
    for (Digraph::NodeIt n(NFA); n != INVALID; ++n) origNodes.push_back(n);

    std::vector<Digraph::Arc> origArcs;
    for (Digraph::ArcIt a(NFA); a != INVALID; ++a) origArcs.push_back(a);

    const int maxN = NFA.maxNodeId() + 1;
    const int maxA = NFA.maxArcId() + 1;

    node_t s = subgraph.start;
    node_t t = subgraph.end;

    // --- Reachability from s
    Dfs<Digraph> dfs_f(NFA);
    dfs_f.run(s);

    // --- Reachability to t via reverse graph
    ReverseDigraph<const Digraph> RG(NFA);
    Dfs<ReverseDigraph<const Digraph>> dfs_b(RG);
    dfs_b.run(t);

    // --- Mark nodes on any s->t path
    Digraph::NodeMap<bool> on_path(NFA, false);
    for (auto n : origNodes) {
        if (dfs_f.reached(n) && dfs_b.reached(n)) on_path[n] = true;
    }

    InplaceDuplicateResult res;
    res.old2newNode.assign(maxN, Digraph::Node(INVALID));
    res.old2newArc.assign(maxA,  Digraph::Arc(INVALID));

    // --- Duplicate nodes
    for (auto n : origNodes) {
        if (!on_path[n]) continue;
        auto nn = NFA.addNode();
        res.old2newNode[NFA.id(n)] = nn;
        node_map[nn] = node_map[n];
    }

    res.s_copy = res.old2newNode[NFA.id(s)];
    res.t_copy = res.old2newNode[NFA.id(t)];

    // --- Duplicate arcs whose endpoints are both on_path
    for (auto a : origArcs) {
        auto u = NFA.source(a);
        auto v = NFA.target(a);
        if (!on_path[u] || !on_path[v]) continue;

        auto uu = res.old2newNode[NFA.id(u)];
        auto vv = res.old2newNode[NFA.id(v)];
        arc_t aa = update_arc_map(NFA, node_map, arc_map, uu, vv);
        res.old2newArc[NFA.id(a)] = aa;
    }
    subgraph_copy.start = res.s_copy;
    subgraph_copy.end = res.t_copy;
    subgraph_copy.copyMeta(subgraph);
}


bool redundancy_test(buffer_t &buffer)
{
    if(buffer.size() == 1) return false;
    int sym2 = buffer.top();
    buffer.pop();
    int sym1 = buffer.top();
    buffer.push(sym2);
    bool redundant = (sym1 == sym2);
    return redundant;
}


bool twin_test(const Subgraph &node_pair)
{
    return node_pair.start == node_pair.end;
}


Subgraph add_node(nfa_t &nfa, lmap_t &node_map, const int &symbol)
{
    node_t node = nfa.addNode();
    node_map[node] = symbol;
    Subgraph twins{node, node, 0u, 1u, {1u}, Default}; // Source and target are the same node
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
    Subgraph twins = {node, node, 0u, 1u, {1u}, Default}; // Source and target are the same node
    stack.push(twins); // Basically just push a pair of dummy twins onto the stack
}


void concat_procedure(nfa_t &nfa, lmap_t &node_map, nfa_stack_t &stack, amap_t &arc_map, buffer_t &buffer, catsites_t& cats)
{
    Subgraph subgraph2 = stack.top();
    stack.pop();
    Subgraph subgraph1 = stack.top();
    stack.pop();
    if(twin_test(subgraph2)) twin_procedure(subgraph2, buffer, nfa, stack, node_map);
    if(twin_test(subgraph1)) twin_procedure(subgraph1, buffer, nfa, stack, node_map);
    arc_t cat_arc = update_arc_map(nfa, node_map, arc_map, subgraph1.end, subgraph2.start);
    Subgraph subgraph{subgraph1.start, subgraph2.end};
    subgraph.concatInfo(subgraph1, subgraph2);
    detect_bad_graphs(subgraph1, subgraph2, subgraph, nfa, cat_arc, cats);
    stack.push(subgraph);
}


void union_procedure(nfa_t &nfa, nfa_stack_t &stack, lmap_t &node_map, amap_t &arc_map, buffer_t &buffer)
{
    Subgraph subgraph2 = stack.top();
    stack.pop();
    Subgraph subgraph1 = stack.top();
    stack.pop();
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

    Subgraph new_subgraph{split_node, ghost_node};
    new_subgraph.unionInfo(subgraph1, subgraph2);
    stack.push(new_subgraph);
}


void optional_procedure(nfa_t &nfa, nfa_stack_t &stack, lmap_t &node_map, amap_t &arc_map, buffer_t &buffer)
{
    Subgraph subgraph = stack.top();
    stack.pop();
    if(twin_test(subgraph)) twin_procedure(subgraph, buffer, nfa, stack, node_map);

    node_t split_node = nfa.addNode();
    node_map[split_node] = Split;
    update_arc_map(nfa, node_map, arc_map, split_node, subgraph.start);

    node_t ghost_node = nfa.addNode();
    node_map[ghost_node] = Ghost;
    update_arc_map(nfa, node_map, arc_map, split_node, ghost_node);
    update_arc_map(nfa, node_map, arc_map, subgraph.end, ghost_node);

    Subgraph new_subgraph{split_node, ghost_node};
    new_subgraph.optionInfo(subgraph);
    stack.push(new_subgraph);
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
        copy_subgraph(subgraph, nfa, node_map, new_subgraph, arc_map, buffer);
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
    Subgraph stack_subgraph{split_node, ghost_node};
    stack_subgraph.kleeneInfo(subgraph, k);
    stack.push(stack_subgraph);
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
        copy_subgraph(subgraph, nfa, node_map, new_subgraph, arc_map, buffer);
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


bool quant_procedure(nfa_t &nfa, nfa_stack_t &stack, lmap_t &node_map, const uint8_t &k, amap_t &arc_map, const size_t min, const size_t max, catsites_t& cats, buffer_t &buffer)
{
    bool skip = false;
    if(min == 0) // Quants like {0,4}
    {
        kleene_procedure(nfa, stack, node_map, (max+1), arc_map, buffer);
        if(stack.size() != 1)
        {
            concat_procedure(nfa, node_map, stack, arc_map, buffer, cats);
            skip = true;
        } 
        return skip;
    }
    Subgraph subgraph = stack.top();
    int sym = buffer.top();
    bool twins = twin_test(subgraph);
    if(stack.size() != 1) // If there's one item on the stack then the subgraph being copied is the first subgraph in the regex
    {
        concat_procedure(nfa, node_map, stack, arc_map, buffer, cats);
        if(twin_test(subgraph)) default_procedure(buffer, sym, stack);
        skip = true; // If it's not first then skip the next concat operator
    }
    size_t extra = (max == 0) ? 0 : (max-min);
    for(size_t i = 1; i < min; ++i)
    {
        Subgraph new_subgraph;
        new_subgraph.copyMeta(subgraph);
        // copy_subgraph(subgraph, nfa, node_map, new_subgraph, arc_map, buffer);
        // stack.push(new_subgraph);
        concat_procedure(nfa, node_map, stack, arc_map, buffer, cats);
    }
    for(size_t i = 0; i < extra; ++i)
    {
        Subgraph new_subgraph;
        copy_subgraph(subgraph, nfa, node_map, new_subgraph, arc_map, buffer);
        stack.push(new_subgraph);
        optional_procedure(nfa, stack, node_map, arc_map, buffer);
        concat_procedure(nfa, node_map, stack, arc_map, buffer, cats);
    }
    return skip;
}


catsites_t construct_reduced_kgraph(const std::string &postfix, nfa_t &nfa, lmap_t &node_map, amap_t &arc_map, const uint8_t &k)
{
    nfa_stack_t stack;
    buffer_t buffer;
    node_t start_node = nfa.addNode(); // I don't know why, but a buffer node is necessary for top sort...
    node_map[start_node] = Ghost;
    std::pair<size_t, size_t> min_max;
    catsites_t catsites;
    bool skip = false;
    for(size_t i = 0; i < postfix.size(); i++)
    {
        int symbol = postfix[i];
        if(std::isdigit(symbol)) continue;
        // DBG(static_cast<char>(symbol));
        switch(symbol)
        {
            default: // Character
                // default_procedure(nfa, stack, node_map, symbol);
                default_procedure(buffer, symbol, stack);
                break;
            case '-': // Concat
                if(skip)
                {
                    skip = false;
                    continue;
                }
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
            case '{': // Start of Quantifier
                min_max = parse_quant(postfix, i);
                if(min_max == OPT_QUANT) // A special case
                {
                    optional_procedure(nfa, stack, node_map, arc_map, buffer);
                    break;
                }
                skip = quant_procedure(nfa, stack, node_map, k, arc_map, min_max.first, min_max.second, catsites, buffer);
                break;
            case '}': // End of Quantifier
            case ',':
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
    // stack.pop();
    assert(stack.size() == 0);
    return catsites;
}
