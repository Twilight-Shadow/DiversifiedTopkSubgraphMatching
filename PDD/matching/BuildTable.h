#ifndef BUILDTABLE_H
#define BUILDTABLE_H

#include "graph/MatchGraph.h"
#include <vector>

class BuildTable {
public:
    static void buildTables(const Graph* data_graph, const Graph* query_graph, ui** candidates, ui* candidates_count,
                            Edges*** edge_matrix);

    static void printTableCardinality(const Graph* query_graph, Edges*** edge_matrix);
    static void printTableCardinality(const Graph *query_graph, TreeNode *tree, ui *order,
                                         std::vector<std::unordered_map<VertexID, std::vector<VertexID >>> &TE_Candidates,
                                         std::vector<std::vector<std::unordered_map<VertexID, std::vector<VertexID>>>> &NTE_Candidates);
    static void printTableInfo(const Graph* query_graph, Edges*** edge_matrix);
    static void printTableInfo(VertexID begin_vertex, VertexID end_vertex, Edges*** edge_matrix);
    static size_t computeMemoryCostInBytes(const Graph* query_graph, ui* candidates_count, Edges*** edge_matrix);
    static size_t computeMemoryCostInBytes(const Graph *query_graph, ui *candidates_count, ui *order, TreeNode *tree,
                                               std::vector<std::unordered_map<VertexID, std::vector<VertexID >>> &TE_Candidates,
                                               std::vector<std::vector<std::unordered_map<VertexID, std::vector<VertexID>>>> &NTE_Candidates);
};

#endif //BUILDTABLE_H