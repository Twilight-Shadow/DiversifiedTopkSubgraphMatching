#ifndef INTRA_MATCH_H
#define INTRA_MATCH_H

#include <vector>
#include <queue>
#include <bitset>
#include <map>

#include "graph/MatchGraph.h"

class IntraPartitionFilter {
public:
    static bool NLFFilter(const Graph* data_graph, const Graph* query_graph, ui** &candidates, ui* &candidates_count);
    static bool Filter(const Graph *data_graph, const Graph *query_graph, ui **&candidates, ui *&candidates_count);
    static void computeCandidateWithNLF(const Graph *data_graph, const Graph *query_graph, VertexID query_vertex,
                                        ui &count, ui *buffer = NULL);
    static void generateCandidates(const Graph *data_graph, const Graph *query_graph, VertexID query_vertex,
                                      VertexID *pivot_vertices, ui pivot_vertices_count, VertexID **candidates,
                                      ui *candidates_count, ui *flag, ui *updated_flag);
    static void sortCandidates(ui** candidates, ui* candidates_count, ui num);
private:
    static void allocateBuffer(const Graph* data_graph, const Graph* query_graph, ui** &candidates, ui* &candidates_count);
    static bool verifyExactTwigIso(const Graph *data_graph, const Graph *query_graph, ui data_vertex, ui query_vertex,
                                   bool **valid_candidates, int *left_to_right_offset, int *left_to_right_edges,
                                   int *left_to_right_match, int *right_to_left_match, int* match_visited,
                                   int* match_queue, int* match_previous);
    static void compactCandidates(ui** &candidates, ui* &candidates_count, ui query_vertex_num);
    static bool isCandidateSetValid(ui** &candidates, ui* &candidates_count, ui query_vertex_num);
};

class IntraPartitionQueryPlan {
public:
    static void generateQueryPlan(const Graph *data_graph, const Graph *query_graph, ui *candidates_count,
                                         ui *&order, ui *&pivot);
    static void checkQueryPlanCorrectness(const Graph* query_graph, ui* order, ui* pivot);
private:
    static VertexID selectStartVertex(const Graph *query_graph, ui *candidates_count);
    static void updateValidVertices(const Graph* query_graph, VertexID query_vertex, std::vector<bool>& visited, std::vector<bool>& adjacent);
};

class IntraPartitionEnumeration {
public:
    static size_t LFTJ(const Graph *data_graph, const Graph *query_graph, Edges ***edge_matrix, ui **candidates, ui *candidates_count,
                           ui *order, size_t output_limit_num, size_t &call_count, pair<VertexID, VertexID>* ret_embedding);
private:
    static void generateBN(const Graph *query_graph, ui *order, ui **&bn, ui *&bn_count);
    static void allocateBuffer(const Graph *query_graph, const Graph *data_graph, ui *candidates_count, ui *&idx,
                                   ui *&idx_count, ui *&embedding, ui *&idx_embedding, ui *&temp_buffer,
                                   ui **&valid_candidate_idx, bool *&visited_vertices);
    static void releaseBuffer(ui query_vertices_num, ui *idx, ui *idx_count, ui *embedding, ui *idx_embedding,
                                  ui *temp_buffer, ui **valid_candidate_idx, bool *visited_vertices, ui **bn, ui *bn_count);
    static void generateValidCandidateIndex(ui depth, ui *idx_embedding, ui *idx_count, ui **valid_candidate_index,
                                                Edges ***edge_matrix, ui **bn, ui *bn_cnt, ui *order, ui *&temp_buffer);
};

#endif //INTRA_MATCH_H