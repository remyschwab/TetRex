#include "construction_tools.h"



struct InplaceDuplicateResult
{
    std::vector<lemon::SmartDigraph::Node> old2newNode; // size = old maxNodeId+1
    std::vector<lemon::SmartDigraph::Arc>  old2newArc;  // size = old maxArcId+1
    lemon::SmartDigraph::Node s_copy{lemon::INVALID};
    lemon::SmartDigraph::Node t_copy{lemon::INVALID};
};

void copy_subgraph(const Subgraph& subgraph, nfa_t& NFA, lmap_t& node_map, Subgraph& subgraph_copy, amap_t& arc_map)
{
    using namespace lemon;
    using Digraph = SmartDigraph;

    if(subgraph.start == subgraph.end) // Copying is very simple if the subgraph is just a single node
    {
        node_t new_twin = NFA.addNode();
        node_map[new_twin] = node_map[subgraph.start];
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

std::pair<size_t, size_t> parse_quant(const std::string& postfix, size_t quant_start)
{ // Parse quantifiers that look like {3,4} or {4} BUT NOT {3,}
    std::pair<size_t, size_t> min_max = {0,0};
    size_t comma = postfix.find(',', quant_start);
    size_t end = postfix.find('}', quant_start);
    // Check to see if there is a comma at all or if it might belong to another quant
    if(comma == std::string::npos || comma > end)
    {
        min_max.first = std::stoi(postfix.substr(quant_start+1, end - quant_start));
        return min_max;
    }
    min_max.first = std::stoi(postfix.substr(quant_start+1, comma - quant_start));
    min_max.second = std::stoi(postfix.substr(comma + 1, end - comma - 1));
    return min_max;
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
