#include "construction_tools.h"


void paste_to_graph(nfa_t &NFA, lmap_t &node_map, const node_t &reference_node, node_t &paste_node)
{
    int symbol = node_map[reference_node];
    paste_node = NFA.addNode();
    node_map[paste_node] = symbol;
}

size_t copy_subgraph(Subgraph &subgraph, nfa_t &NFA, lmap_t &node_map, Subgraph &subgraph_copy, amap_t &arc_map)
{ // Lord help anyone who ever needs to debug this method
    node_t new_node;
    size_t added_node_count = 0;
    // paste_to_graph(NFA, node_map, subgraph.first, new_node);
    paste_to_graph(NFA, node_map, subgraph.start, new_node);
    ++added_node_count;
    // subgraph_copy.first = new_node;
    subgraph_copy.start = new_node;
    // If the operand is just a single character then just copy that one node
    // if(subgraph.first == subgraph.second)
    if(subgraph.start == subgraph.end)
    {
        // subgraph_copy.second = new_node;
        subgraph_copy.end = new_node;
        return added_node_count;
    }
    // If the operand is a more complicated subgraph, then traverse with DFS copying nodes and arcs along the way
    nfa_t::NodeMap<node_t> reference_to_copy_map(NFA); // A mapping of nodes in the old subgraph to the new one
    robin_hood::unordered_set<int> copied_targets; // It's possible that paths converge to the same (ghost) node and we don't want to copy it twice
    
    lemon::Dfs<nfa_t> dfs(NFA);
    dfs.init();
    // dfs.addSource(subgraph.first);
    dfs.addSource(subgraph.start);
    
    // reference_to_copy_map[subgraph.first] = new_node;
    reference_to_copy_map[subgraph.start] = new_node;
    while(!dfs.emptyQueue())
    {
        arc_t arc = dfs.processNextArc();
        node_t source = NFA.source(arc);
        // if(source == subgraph.second) break;
        if(source == subgraph.end) break; // i don't totally get why this works
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
    subgraph_copy.end = new_node;
    return added_node_count;
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


void print_graph(nfa_t &NFA, lmap_t &nmap, const catsites_t& cats, const bool& augment)
{
    std::fstream f;
    const std::string filename = "kgraph_visualizer.gv";
    f.open(filename, std::ios::out);
    f << "digraph kGraph\n{\n\trankdir=\"LR\";\n";
    // Collect and style all the nodes
    nfa_t::ArcMap<bool> filter(NFA, true);
    if(augment) for(auto && cat: cats) filter[cat.arc_] = false;
    lemon::Bfs<nfa_t>  bfs(NFA);
    bfs.run(NFA.nodeFromId(0));
    for(nfa_t::NodeIt n(NFA); n != lemon::INVALID; ++n)
    {
        const int id = NFA.id(n);
        if(id == 0)
        {
            f << "\t" << NFA.id(n) << " [shape=point label=\"\"];\n";
            continue;
        }
        if(nmap[n] == Split)
        {
            f << "\t" << NFA.id(n) << " [label=\"" << "Ø\"];\n";
            continue;
        }
        if(nmap[n] == Ghost)
        {
            f << "\t" << NFA.id(n) << " [label=\"•\"];\n";
            continue;
        }
        if(nmap[n] == Match)
        {
            f << "\t" << NFA.id(n) << " [shape=doublecircle label=\"\"];\n";
            continue;
        }
        if(nmap[n] == Gap)
        {
            f << "\t" << NFA.id(n) << " [label=\"GAP\"];\n";
            continue;
        }
        f << "\t" << NFA.id(n) << " [label=\"" << (char)nmap[n] << "\"];\n";
    }
    // Traverse the Graph and get all the transitions
    lemon::Dfs<nfa_t> dfs(NFA);
    dfs.init();
    dfs.addSource(NFA.nodeFromId(0));
    while (!dfs.emptyQueue())
    {
        nfa_t::Arc arc = dfs.processNextArc();
        if(filter[arc]) f << "\t" << NFA.id(NFA.source(arc)) << "->" << NFA.id(NFA.target(arc)) << ";" << std::endl;
    }
    f << "}";
    f.close();
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


arc_t update_arc_map(nfa_t &NFA, lmap_t &node_map, amap_t &arc_map, const node_t &source, const node_t &target)
{
    arc_t new_arc = NFA.addArc(source, target); // Update the internal NFA arc map I guess
    int source_id = NFA.id(source);
    int symbol = node_map[source];
    if(symbol < 258 || symbol == Gap) // If the node is not a split just place the same target in both slots
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
    return new_arc;
}
