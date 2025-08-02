#ifndef INTER_MATCH_H
#define INTER_MATCH_H

#include <map>
#include <set>
#include <vector>
#include <unordered_set>
#include <bitset>
#include <queue>
#include <unordered_map>

#include "graph/MatchGraph.h"
#include "configuration/Config.h"

namespace InterPartitionFilter {
    bool Filter(const Graph *data_graph, const Graph *query_graph, ui **&candidates, ui *&candidates_count,
                            ui *&order, TreeNode *&tree,   std::vector<std::unordered_map<VertexID, std::vector<VertexID>>> &TE_Candidates,
                            std::vector<std::vector<std::unordered_map<VertexID, std::vector<VertexID>>>> &NTE_Candidates, 
                            VertexID root, std::unordered_map<uint32_t, std::vector<VertexID> >& getReplicabyLabel);

    void computeCandidateWithNLF(const Graph *data_graph, const Graph *query_graph, VertexID query_vertex,
                                        ui &count, ui *buffer = NULL);

    void computeRootCandidateWithNLF(const Graph *data_graph, const Graph *query_graph, VertexID query_vertex, ui &count, 
                                            std::unordered_map<uint32_t, std::vector<VertexID> >& getReplicabyLabel, ui *buffer = NULL);

    void sortCandidates(ui** candidates, ui* candidates_count, ui num);

    void allocateBuffer(const Graph* data_graph, const Graph* query_graph, ui** &candidates, ui* &candidates_count);

    void compactCandidates(ui** &candidates, ui* &candidates_count, ui query_vertex_num);

    bool isCandidateSetValid(ui** &candidates, ui* &candidates_count, ui query_vertex_num);
};

namespace GenerateFilteringPlan {
    void generateFilterPlan(const Graph *data_graph, const Graph *query_graph, TreeNode *&tree,
                                       VertexID *&order, VertexID root);

    VertexID selectStartVertex(const Graph *data_graph, const Graph *query_graph);
};

namespace InterPartitionQueryPlan {
    void generateQueryPlan(const Graph* query_graph, TreeNode *tree, ui *bfs_order, ui *&order, ui *&pivot);
};

namespace InterPartitionEnumeration {
    size_t enumerate(const Graph *data_graph, const Graph *query_graph, TreeNode *tree, ui **candidates,
                                      ui *candidates_count,
                                      std::vector<std::unordered_map<VertexID, std::vector<VertexID>>> &TE_Candidates,
                                      std::vector<std::vector<std::unordered_map<VertexID, std::vector<VertexID>>>> &NTE_Candidates,
                                      ui *order, size_t &output_limit_num, size_t &call_count,
                                      const std::unordered_map<VertexID, std::pair<std::set<VertexID>, std::set<VertexID>>>& order_constraints, 
                                      pair<VertexID, VertexID>* ret_embedding);

    void generateValidCandidates(ui depth, ui *embedding, ui *idx_count, ui **valid_candidates, ui *order,
                                            ui *&temp_buffer, TreeNode *tree,
                                            std::vector<std::unordered_map<VertexID, std::vector<VertexID>>> &TE_Candidates,
                                            std::vector<std::vector<std::unordered_map<VertexID, std::vector<VertexID>>>> &NTE_Candidates);
    
    void pruneCandidatesBySymmetryBreaking(ui depth, ui *embedding, ui* order,
                                                  ui *idx_count, ui **valid_candidates,
                                                  const std::unordered_map<VertexID, std::pair<std::set<VertexID>, std::set<VertexID>>>& ordered_constraints);

    void computeAncestor(const Graph *query_graph, VertexID *order, std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> &ancestors);
    
    void AncestorWithSymmetryBreaking(const Graph *query_graph,
                                             std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> &ancestors,
                                             const std::unordered_map<VertexID, std::pair<std::set<VertexID>, std::set<VertexID>>>& ordered_constraints);
};

#endif // INTER_MATCH_H