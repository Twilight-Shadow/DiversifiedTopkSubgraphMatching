#ifndef BASICGRAPH_H
#define BASICGRAPH_H

#include "configuration/Types.h"
#include "MatchGraph.h"
#include <vector>
#include <map>
#include <cstring>
#include <string>
#include <bitset>
#include <unordered_set>

class GraphOP {
public:
    Node** head;
    Node* edge;
    int node_count, edge_count;
    int cnt;
    ui* degree;
public:
    GraphOP();
    GraphOP(pair<ui, ui> max_index);
    ~GraphOP();
    void build(std::vector<pair<VertexID, VertexID> >& edges);
    void add_edge(VertexID u, VertexID v);
    BFS_TREE* pre_process(VertexID start_vertex, Graph* query_graph);
    void Match_BFS(std::vector<VertexID>& start_vertex, BFS_TREE* Qroot, std::vector<ui>& labels, Graph* query_graph, int data_graph_index, GraphOP* QG, std::unordered_set<VertexID>* Graph_Matrix);
    void Match_DFS(BFS_TREE* root, GraphOP* QG, BFS_TREE** table, std::unordered_set<VertexID>* Graph_Matrix, std::map<VertexID, VertexID>& Matchres);
};

#endif // BASICGRAPH_H