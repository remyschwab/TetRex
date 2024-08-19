#include "construction_tools.h"


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
    // If the operand is a more complicated subgraph, then traverse with DFS copying nodes and arcs along the way
    nfa_t::NodeMap<node_t> reference_to_copy_map(NFA); // A mapping of nodes in the old subgraph to the new one
    robin_hood::unordered_set<int> copied_targets; // It's possible that paths converge to the same (ghost) node and we don't want to copy it twice
    
    lemon::Dfs<nfa_t> dfs(NFA);
    dfs.init();
    dfs.addSource(subgraph.first);
    
    reference_to_copy_map[subgraph.first] = new_node;
    while(!dfs.emptyQueue())
    {
        arc_t arc = dfs.processNextArc();
        node_t source = NFA.source(arc);
        if(source == subgraph.second) break; // i don't totally get why this works
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


void print_graph(nfa_t &NFA, lmap_t &nmap)
{
    std::fstream f;
    const std::string filename = "kgraph_visualizer.gv";
    
    // f.open(filename, std::ios::out);
    // f << "digraph kGraph"<<"\n";
    // f << "{"<<"\n";
    // f << "\trankdir=\"LR\"\n";
    
    lemon::Dfs<nfa_t> dfs(NFA);
    dfs.init();
    dfs.addSource(NFA.nodeFromId(0));
    std::vector<lemon::ListPath<nfa_t>> paths;
    lemon::ListPath<nfa_t> currPath;
    nfa_t::Node prevNode = NFA.nodeFromId(0);

    while (!dfs.emptyQueue())
    {
        nfa_t::Arc arc = dfs.processNextArc();

        if (NFA.source(arc) == prevNode)
        {
            currPath.addBack(arc);
            prevNode = NFA.target(arc);
        }
        else
        {
            prevNode = NFA.target(arc);
            currPath = {};
            currPath.addBack(arc);
        }
        if (nmap[NFA.target(arc)] == Match)
        {
            paths.push_back(currPath);
        }
    }

    for (auto path : paths)
    {
        for (int i = 0; i < path.length(); i++)
        {
            std::cout << NFA.id(NFA.source(path.nth(i))) << " -> " << NFA.id(NFA.target(path.nth(i))) << ", ";
        }
        std::cout << std::endl;
    }
    
    // dfs.addSource(start);
    // while(!dfs.emptyQueue())
    // {
    //     arc_t arc = dfs.processNextArc();
    //     node_t source = NFA.source(arc);
    //     if(NFA.id(source) == 0)
    //     {
    //         f << "\t• -> ";
    //     }
    //     node_t target = NFA.target(arc);
    //     if(nmap[target] == Match)
    //     {
    //         f <<  "∆;\n\t";
    //         continue;
    //     }
    //     f << NFA.id(target) << " -> ";
    // }
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
