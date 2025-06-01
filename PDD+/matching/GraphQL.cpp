#include <memory.h>
#include <vector>
#include <algorithm>
#include <ctime>
#include <cstdlib>
#include <random>
#include <limits>
#include <cassert>
#include <cstring>
#include <fstream>

#include "GraphQL.h"
#include "utility/Operation.h"
#include "utility/ComputeSetIntersection.h"

#define INVALID_VERTEX_ID 100000000

using std::pair;
using std::make_pair;

bool FilterVertices::NLFFilter(const Graph *data_graph, const Graph *query_graph, ui **&candidates, ui *&candidates_count) {
    allocateBuffer(data_graph, query_graph, candidates, candidates_count);

    for (ui i = 0; i < query_graph->getVerticesCount(); ++i) {
        VertexID query_vertex = i;
        computeCandidateWithNLF(data_graph, query_graph, query_vertex, candidates_count[query_vertex], candidates[query_vertex]);

        if (candidates_count[query_vertex] == 0) {
            return false;
        }
    }

    return true;
}

bool FilterVertices::GQLFilter(const Graph *data_graph, const Graph *query_graph, ui **&candidates, ui *&candidates_count) {
    // Local refinement.
    if (!NLFFilter(data_graph, query_graph, candidates, candidates_count))
        return false;

    // Allocate buffer.
    ui query_vertex_num = query_graph->getVerticesCount();
    ui data_vertex_num = data_graph->getVerticesCount();

    bool** valid_candidates = new bool*[query_vertex_num];
    for (ui i = 0; i < query_vertex_num; ++i) {
        valid_candidates[i] = new bool[data_vertex_num];
        memset(valid_candidates[i], 0, sizeof(bool) * data_vertex_num);
    }

    ui query_graph_max_degree = query_graph->getGraphMaxDegree();
    ui data_graph_max_degree = data_graph->getGraphMaxDegree();

    int* left_to_right_offset = new int[query_graph_max_degree + 1];
    int* left_to_right_edges = new int[query_graph_max_degree * data_graph_max_degree];
    int* left_to_right_match = new int[query_graph_max_degree];
    int* right_to_left_match = new int[data_graph_max_degree];
    int* match_visited = new int[data_graph_max_degree + 1];
    int* match_queue = new int[query_vertex_num];
    int* match_previous = new int[data_graph_max_degree + 1];

    // Record valid candidate vertices for each query vertex.
    for (ui i = 0; i < query_vertex_num; ++i) {
        VertexID query_vertex = i;
        for (ui j = 0; j < candidates_count[query_vertex]; ++j) {
            VertexID data_vertex = candidates[query_vertex][j];
            valid_candidates[query_vertex][data_vertex] = true;
        }
    }

    // Global refinement.
    for (ui l = 0; l < 2; ++l) {
        for (ui i = 0; i < query_vertex_num; ++i) {
            VertexID query_vertex = i;
            for (ui j = 0; j < candidates_count[query_vertex]; ++j) {
                VertexID data_vertex = candidates[query_vertex][j];

                if (data_vertex == INVALID_VERTEX_ID)
                    continue;

                if (!verifyExactTwigIso(data_graph, query_graph, data_vertex, query_vertex, valid_candidates,
                                        left_to_right_offset, left_to_right_edges, left_to_right_match,
                                        right_to_left_match, match_visited, match_queue, match_previous)) {
                    candidates[query_vertex][j] = INVALID_VERTEX_ID;
                    valid_candidates[query_vertex][data_vertex] = false;
                }
            }
        }
    }

    // Compact candidates.
    compactCandidates(candidates, candidates_count, query_vertex_num);

    // Release memory.
    for (ui i = 0; i < query_vertex_num; ++i) {
        delete[] valid_candidates[i];
    }
    delete[] valid_candidates;
    delete[] left_to_right_offset;
    delete[] left_to_right_edges;
    delete[] left_to_right_match;
    delete[] right_to_left_match;
    delete[] match_visited;
    delete[] match_queue;
    delete[] match_previous;

    return isCandidateSetValid(candidates, candidates_count, query_vertex_num);
}

void FilterVertices::sortCandidates(ui **candidates, ui *candidates_count, ui num) {
    for (ui i = 0; i < num; ++i) {
        std::sort(candidates[i], candidates[i] + candidates_count[i]);
    }
}

void FilterVertices::allocateBuffer(const Graph *data_graph, const Graph *query_graph, ui **&candidates,
                                    ui *&candidates_count) {
    ui query_vertex_num = query_graph->getVerticesCount();
    ui candidates_max_num = data_graph->getGraphMaxLabelFrequency();

    candidates_count = new ui[query_vertex_num];
    memset(candidates_count, 0, sizeof(ui) * query_vertex_num);

    candidates = new ui*[query_vertex_num];

    for (ui i = 0; i < query_vertex_num; ++i) {
        candidates[i] = new ui[candidates_max_num];
    }
}

void
FilterVertices::computeCandidateWithNLF(const Graph *data_graph, const Graph *query_graph, VertexID query_vertex,
                                               ui &count, ui *buffer) {
    LabelID label = query_graph->getVertexLabel(query_vertex);
    ui degree = query_graph->getVertexDegree(query_vertex);
    const std::unordered_map<LabelID, ui>* query_vertex_nlf = query_graph->getVertexNLF(query_vertex);
    ui data_vertex_num;
    const ui* data_vertices = data_graph->getVerticesByLabel(label, data_vertex_num);
    count = 0;
    for (ui j = 0; j < data_vertex_num; ++j) {
        ui data_vertex = data_vertices[j];
        if (data_graph->getVertexDegree(data_vertex) >= degree) {
            // NFL check
            const std::unordered_map<LabelID, ui>* data_vertex_nlf = data_graph->getVertexNLF(data_vertex);

            if (data_vertex_nlf->size() >= query_vertex_nlf->size()) {
                bool is_valid = true;

                for (auto element : *query_vertex_nlf) {
                    auto iter = data_vertex_nlf->find(element.first);
                    if (iter == data_vertex_nlf->end() || iter->second < element.second) {
                        is_valid = false;
                        break;
                    }
                }

                if (is_valid) {
                    if (buffer != NULL) {
                        buffer[count] = data_vertex;
                    }
                    count += 1;
                }
            }
        }
    }
}

bool FilterVertices::verifyExactTwigIso(const Graph *data_graph, const Graph *query_graph, ui data_vertex, ui query_vertex,
                                   bool **valid_candidates, int *left_to_right_offset, int *left_to_right_edges,
                                   int *left_to_right_match, int *right_to_left_match, int* match_visited,
                                   int* match_queue, int* match_previous) {
    // Construct the bipartite graph between N(query_vertex) and N(data_vertex)
    ui left_partition_size;
    ui right_partition_size;
    const VertexID* query_vertex_neighbors = query_graph->getVertexNeighbors(query_vertex, left_partition_size);
    const VertexID* data_vertex_neighbors = data_graph->getVertexNeighbors(data_vertex, right_partition_size);

    ui edge_count = 0;
    for (int i = 0; i < left_partition_size; ++i) {
        VertexID query_vertex_neighbor = query_vertex_neighbors[i];
        left_to_right_offset[i] = edge_count;

        for (int j = 0; j < right_partition_size; ++j) {
            VertexID data_vertex_neighbor = data_vertex_neighbors[j];

            if (valid_candidates[query_vertex_neighbor][data_vertex_neighbor]) {
                left_to_right_edges[edge_count++] = j;
            }
        }
    }
    left_to_right_offset[left_partition_size] = edge_count;

    memset(left_to_right_match, -1, left_partition_size * sizeof(int));
    memset(right_to_left_match, -1, right_partition_size * sizeof(int));

    GraphOperations::match_bfs(left_to_right_offset, left_to_right_edges, left_to_right_match, right_to_left_match,
                               match_visited, match_queue, match_previous, left_partition_size, right_partition_size);
    for (int i = 0; i < left_partition_size; ++i) {
        if (left_to_right_match[i] == -1)
            return false;
    }

    return true;
}

void FilterVertices::compactCandidates(ui **&candidates, ui *&candidates_count, ui query_vertex_num) {
    // std::mt19937 gen(std::random_device{}());
    for (ui i = 0; i < query_vertex_num; ++i) {
        VertexID query_vertex = i;
        ui next_position = 0;
        for (ui j = 0; j < candidates_count[query_vertex]; ++j) {
            VertexID data_vertex = candidates[query_vertex][j];

            if (data_vertex != INVALID_VERTEX_ID) {
                candidates[query_vertex][next_position++] = data_vertex;
            }
        }

        candidates_count[query_vertex] = next_position;
        // if (next_position > 1) {
        //     std::shuffle(candidates[query_vertex], candidates[query_vertex] + next_position, gen);
        // }
    }
}

bool FilterVertices::isCandidateSetValid(ui **&candidates, ui *&candidates_count, ui query_vertex_num) {
    for (ui i = 0; i < query_vertex_num; ++i) {
        if (candidates_count[i] == 0)
            return false;
    }
    return true;
}

void GenerateQueryPlan::generateGQLQueryPlan(const Graph *data_graph, const Graph *query_graph, ui *candidates_count,
                                             ui *&order, ui *&pivot) {
     /**
      * Select the vertex v such that (1) v is adjacent to the selected vertices; and (2) v has the minimum number of candidates.
      */
     std::vector<bool> visited_vertices(query_graph->getVerticesCount(), false);
     std::vector<bool> adjacent_vertices(query_graph->getVerticesCount(), false);
     order = new ui[query_graph->getVerticesCount()];
     pivot = new ui[query_graph->getVerticesCount()];
    //  srand(time(0));
     VertexID start_vertex = selectGQLStartVertex(query_graph, candidates_count);
    //  VertexID start_vertex = rand() % (query_graph->getVerticesCount());
    //  std::cout << "........" << start_vertex << "........\n";
    //  std::cout << "......." << query_graph->getVerticesCount() << ".........\n";
     order[0] = start_vertex;
     updateValidVertices(query_graph, start_vertex, visited_vertices, adjacent_vertices);

     for (ui i = 1; i < query_graph->getVerticesCount(); ++i) {
          VertexID next_vertex;
          ui min_value = data_graph->getVerticesCount() + 1;
          for (ui j = 0; j < query_graph->getVerticesCount(); ++j) {
               VertexID cur_vertex = j;

               if (!visited_vertices[cur_vertex] && adjacent_vertices[cur_vertex]) {
                    if (candidates_count[cur_vertex] < min_value) {
                         min_value = candidates_count[cur_vertex];
                         next_vertex = cur_vertex;
                    }
                    else if (candidates_count[cur_vertex] == min_value && query_graph->getVertexDegree(cur_vertex) > query_graph->getVertexDegree(next_vertex)) {
                         next_vertex = cur_vertex;
                    }
               }
          }
          updateValidVertices(query_graph, next_vertex, visited_vertices, adjacent_vertices);
          order[i] = next_vertex;
     }

     // Pick a pivot randomly.
     for (ui i = 1; i < query_graph->getVerticesCount(); ++i) {
         VertexID u = order[i];
         for (ui j = 0; j < i; ++j) {
             VertexID cur_vertex = order[j];
             if (query_graph->checkEdgeExistence(u, cur_vertex)) {
                 pivot[i] = cur_vertex;
                 break;
             }
         }
     }
}

VertexID GenerateQueryPlan::selectGQLStartVertex(const Graph *query_graph, ui *candidates_count) {
    /**
     * Select the vertex with the minimum number of candidates as the start vertex.
     * Tie Handling:
     *  1. degree
     *  2. label id
     */

     ui start_vertex = 0;

     for (ui i = 1; i < query_graph->getVerticesCount(); ++i) {
          VertexID cur_vertex = i;

          if (candidates_count[cur_vertex] < candidates_count[start_vertex]) {
               start_vertex = cur_vertex;
          }
          else if (candidates_count[cur_vertex] == candidates_count[start_vertex]
                   && query_graph->getVertexDegree(cur_vertex) > query_graph->getVertexDegree(start_vertex)) {
               start_vertex = cur_vertex;
          }
     }

     return start_vertex;
}

void GenerateQueryPlan::updateValidVertices(const Graph *query_graph, VertexID query_vertex, std::vector<bool> &visited,
                                            std::vector<bool> &adjacent) {
     visited[query_vertex] = true;
     ui nbr_cnt;
     const ui* nbrs = query_graph->getVertexNeighbors(query_vertex, nbr_cnt);

     for (ui i = 0; i < nbr_cnt; ++i) {
          ui nbr = nbrs[i];
          adjacent[nbr] = true;
     }
}

void GenerateQueryPlan::checkQueryPlanCorrectness(const Graph *query_graph, ui *order, ui *pivot) {
    ui query_vertices_num = query_graph->getVerticesCount();
    std::vector<bool> visited_vertices(query_vertices_num, false);
    // Check whether each query vertex is in the order.
    for (ui i = 0; i < query_vertices_num; ++i) {
        VertexID vertex = order[i];
        assert(vertex < query_vertices_num && vertex >= 0);

        visited_vertices[vertex] = true;
    }

    for (ui i = 0; i < query_vertices_num; ++i) {
        VertexID vertex = i;
        assert(visited_vertices[vertex]);
    }

    // Check whether the order is connected.

    std::fill(visited_vertices.begin(), visited_vertices.end(), false);
    visited_vertices[order[0]] = true;
    for (ui i = 1; i < query_vertices_num; ++i) {
        VertexID vertex = order[i];
        VertexID pivot_vertex = pivot[i];
        assert(query_graph->checkEdgeExistence(vertex, pivot_vertex));
        assert(visited_vertices[pivot_vertex]);
        visited_vertices[vertex] = true;
    }
}

size_t
EvaluateQuery::LFTJ(const Graph *data_graph, const Graph *query_graph, Edges ***edge_matrix, ui **candidates,
                    ui *candidates_count,
                    ui *order, size_t output_limit_num, size_t &call_count, int data_graph_index, VertexID* remap_info) {

    // Generate bn.
    ui **bn;
    ui *bn_count;
    generateBN(query_graph, order, bn, bn_count);

    // Allocate the memory buffer.
    ui *idx;
    ui *idx_count;
    ui *embedding;
    ui *idx_embedding;
    ui *temp_buffer;
    ui **valid_candidate_idx;
    bool *visited_vertices;
    // bool prt = 0;
    allocateBuffer(data_graph, query_graph, candidates_count, idx, idx_count, embedding, idx_embedding,
                   temp_buffer, valid_candidate_idx, visited_vertices);

    size_t embedding_cnt = 0;
    int cur_depth = 0;
    int max_depth = query_graph->getVerticesCount();
    VertexID start_vertex = order[0];

    idx[cur_depth] = 0;
    idx_count[cur_depth] = candidates_count[start_vertex];

    for (ui i = 0; i < idx_count[cur_depth]; ++i) {
        valid_candidate_idx[cur_depth][i] = i;
    }
    // bool* used = new bool[query_graph->getVerticesCount() + 1];
    // for (int i = 0; i < query_graph->getVerticesCount(); ++i) {
    //     used[i] = false;
    // }
    // for (int i = 0; i < query_graph->getVerticesCount(); ++i) {
    //     std::cout << order[i] << '\n';
    // }

    std::ofstream ofile("./MatchRes/Matching" + std::to_string(data_graph_index) + ".txt", std::ios::out);
    std::vector<pair<VertexID, VertexID> > MStack;
    while (true) {
        while (idx[cur_depth] < idx_count[cur_depth]) {
            ui valid_idx = valid_candidate_idx[cur_depth][idx[cur_depth]];
            VertexID u = order[cur_depth];
            VertexID v = candidates[u][valid_idx];

            if (visited_vertices[v]) {
                idx[cur_depth] += 1;
                continue;
            }

            embedding[u] = v;
            MStack.push_back(make_pair(u, v));
            // ofile << "Matched vertices: " << u << ' ' << remap_info[v] << '\n';
            // if(!prt) {
            //     prt = 1;
            //     std::cout << "Query anchor: " << u << ", Data anchor: " << v << '\n';
            // }
            idx_embedding[u] = valid_idx;
            visited_vertices[v] = true;
            idx[cur_depth] += 1;

            if (cur_depth == max_depth - 1) {
                embedding_cnt += 1;
                visited_vertices[v] = false;
                // ofile << "----------------cnt_" << embedding_cnt << "----------------\n";
                if (embedding_cnt >= output_limit_num) {
                    goto EXIT;
                }
            } else {
                call_count += 1;
                cur_depth += 1;

                idx[cur_depth] = 0;
                generateValidCandidateIndex(cur_depth, idx_embedding, idx_count, valid_candidate_idx, edge_matrix, bn,
                                            bn_count, order, temp_buffer);
            }
        }
        if(!MStack.empty()) MStack.pop_back();
        cur_depth -= 1;
        if (cur_depth < 0)
            break;
        else {
            VertexID u = order[cur_depth];
            visited_vertices[embedding[u]] = false;
        }
    }

    // Release the buffer.
    EXIT:
    for (auto i = MStack.begin(); i != MStack.end(); ++i) {
        ofile << "Matched vertices: " << i->first << ' ' << remap_info[i->second] << '\n';
    }
    ofile.close();
    MStack.clear();
    releaseBuffer(max_depth, idx, idx_count, embedding, idx_embedding, temp_buffer, valid_candidate_idx,
                  visited_vertices,
                  bn, bn_count);

    return embedding_cnt;
}

void EvaluateQuery::generateBN(const Graph *query_graph, ui *order, ui **&bn, ui *&bn_count) {
    ui query_vertices_num = query_graph->getVerticesCount();
    bn_count = new ui[query_vertices_num];
    std::fill(bn_count, bn_count + query_vertices_num, 0);
    bn = new ui *[query_vertices_num];
    for (ui i = 0; i < query_vertices_num; ++i) {
        bn[i] = new ui[query_vertices_num];
    }

    std::vector<bool> visited_vertices(query_vertices_num, false);
    visited_vertices[order[0]] = true;
    for (ui i = 1; i < query_vertices_num; ++i) {
        VertexID vertex = order[i];

        ui nbrs_cnt;
        const ui *nbrs = query_graph->getVertexNeighbors(vertex, nbrs_cnt);
        for (ui j = 0; j < nbrs_cnt; ++j) {
            VertexID nbr = nbrs[j];

            if (visited_vertices[nbr]) {
                bn[i][bn_count[i]++] = nbr;
            }
        }

        visited_vertices[vertex] = true;
    }
}

void EvaluateQuery::allocateBuffer(const Graph *data_graph, const Graph *query_graph, ui *candidates_count, ui *&idx,
                              ui *&idx_count, ui *&embedding, ui *&idx_embedding, ui *&temp_buffer,
                              ui **&valid_candidate_idx, bool *&visited_vertices) {
    ui query_vertices_num = query_graph->getVerticesCount();
    ui data_vertices_num = data_graph->getVerticesCount();
    ui max_candidates_num = candidates_count[0];

    for (ui i = 1; i < query_vertices_num; ++i) {
        VertexID cur_vertex = i;
        ui cur_candidate_num = candidates_count[cur_vertex];

        if (cur_candidate_num > max_candidates_num) {
            max_candidates_num = cur_candidate_num;
        }
    }

    idx = new ui[query_vertices_num];
    idx_count = new ui[query_vertices_num];
    embedding = new ui[query_vertices_num];
    idx_embedding = new ui[query_vertices_num];
    visited_vertices = new bool[data_vertices_num];
    temp_buffer = new ui[max_candidates_num];
    valid_candidate_idx = new ui *[query_vertices_num];
    for (ui i = 0; i < query_vertices_num; ++i) {
        valid_candidate_idx[i] = new ui[max_candidates_num];
    }

    std::fill(visited_vertices, visited_vertices + data_vertices_num, false);
}

void EvaluateQuery::releaseBuffer(ui query_vertices_num, ui *idx, ui *idx_count, ui *embedding, ui *idx_embedding,
                                  ui *temp_buffer, ui **valid_candidate_idx, bool *visited_vertices, ui **bn,
                                  ui *bn_count) {
    delete[] idx;
    delete[] idx_count;
    delete[] embedding;
    delete[] idx_embedding;
    delete[] visited_vertices;
    delete[] bn_count;
    delete[] temp_buffer;
    for (ui i = 0; i < query_vertices_num; ++i) {
        delete[] valid_candidate_idx[i];
        delete[] bn[i];
    }

    delete[] valid_candidate_idx;
    delete[] bn;
}

void EvaluateQuery::generateValidCandidateIndex(ui depth, ui *idx_embedding, ui *idx_count, ui **valid_candidate_index,
                                                Edges ***edge_matrix, ui **bn, ui *bn_cnt, ui *order,
                                                ui *&temp_buffer) {
    VertexID u = order[depth];
    VertexID previous_bn = bn[depth][0];
    ui previous_index_id = idx_embedding[previous_bn];
    ui valid_candidates_count = 0;

    Edges& previous_edge = *edge_matrix[previous_bn][u];

    valid_candidates_count = previous_edge.offset_[previous_index_id + 1] - previous_edge.offset_[previous_index_id];
    ui* previous_candidates = previous_edge.edge_ + previous_edge.offset_[previous_index_id];

    memcpy(valid_candidate_index[depth], previous_candidates, valid_candidates_count * sizeof(ui));

    ui temp_count;
    for (ui i = 1; i < bn_cnt[depth]; ++i) {
        VertexID current_bn = bn[depth][i];
        Edges& current_edge = *edge_matrix[current_bn][u];
        ui current_index_id = idx_embedding[current_bn];

        ui current_candidates_count = current_edge.offset_[current_index_id + 1] - current_edge.offset_[current_index_id];
        ui* current_candidates = current_edge.edge_ + current_edge.offset_[current_index_id];

        ComputeSetIntersection::ComputeCandidates(current_candidates, current_candidates_count, valid_candidate_index[depth], valid_candidates_count,
                        temp_buffer, temp_count);

        std::swap(temp_buffer, valid_candidate_index[depth]);
        valid_candidates_count = temp_count;
    }

    idx_count[depth] = valid_candidates_count;
}