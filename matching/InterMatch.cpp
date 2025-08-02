#include <vector>
#include <algorithm>
#include <chrono>
#include <memory.h>

#include "InterMatch.h"
#include "utility/Operation.h"
#include "utility/ComputeSetIntersection.h"
#define INVALID_VERTEX_ID 2147483647

bool InterPartitionFilter::Filter(const Graph *data_graph, const Graph *query_graph, ui **&candidates, ui *&candidates_count,
                            ui *&order, TreeNode *&tree,  std::vector<std::unordered_map<VertexID, std::vector<VertexID >>> &TE_Candidates,
                            std::vector<std::vector<std::unordered_map<VertexID, std::vector<VertexID>>>> &NTE_Candidates, 
                            VertexID root, std::unordered_map<uint32_t, std::vector<VertexID> >& getReplicabyLabel) {
    GenerateFilteringPlan::generateFilterPlan(data_graph, query_graph, tree, order, root);

    allocateBuffer(data_graph, query_graph, candidates, candidates_count);

    ui query_vertices_count = query_graph->getVerticesCount();
    ui data_vertices_count = data_graph->getVerticesCount();
    
    computeRootCandidateWithNLF(data_graph, query_graph, root, candidates_count[root], getReplicabyLabel, candidates[root]);

    if (candidates_count[root] == 0) return false;
    
    std::vector<ui> updated_flag(data_vertices_count);
    std::vector<ui> flag(data_vertices_count);
    std::fill(flag.begin(), flag.end(), 0);
    std::vector<bool> visited_query_vertex(query_vertices_count);
    std::fill(visited_query_vertex.begin(), visited_query_vertex.end(), false);

    visited_query_vertex[root] = true;

    TE_Candidates.resize(query_vertices_count);

    for (ui i = 1; i < query_vertices_count; ++i) {
        VertexID u = order[i];
        VertexID u_p = tree[u].parent_;

        ui u_label = query_graph->getVertexLabel(u);
        ui u_degree = query_graph->getVertexDegree(u);
        const std::unordered_map<LabelID, ui>* u_nlf = query_graph->getVertexNLF(u);
        candidates_count[u] = 0;

        visited_query_vertex[u] = true;
        VertexID* frontiers = candidates[u_p];
        ui frontiers_count = candidates_count[u_p];

        for (ui j = 0; j < frontiers_count; ++j) {
            VertexID v_f = frontiers[j];

            if (v_f == INVALID_VERTEX_ID)
                continue;

            ui nbrs_cnt;
            const VertexID* nbrs = data_graph->getVertexNeighbors(v_f, nbrs_cnt);

            auto iter_pair = TE_Candidates[u].emplace(v_f, std::vector<VertexID>());
            for (ui k = 0; k < nbrs_cnt; ++k) {
                VertexID v = nbrs[k];

                if (data_graph->getVertexLabel(v) == u_label && data_graph->getVertexDegree(v) >= u_degree) {
                    const std::unordered_map<LabelID, ui>* v_nlf = data_graph->getVertexNLF(v);

                    if (v_nlf->size() >= u_nlf->size()) {
                        bool is_valid = true;

                        for (auto element : *u_nlf) {
                            auto iter = v_nlf->find(element.first);
                            if (iter == v_nlf->end() || iter->second < element.second) {
                                is_valid = false;
                                break;
                            }
                        }

                        if (is_valid) {
                            iter_pair.first->second.push_back(v);
                            if (flag[v] == 0) {
                                flag[v] = 1;
                                candidates[u][candidates_count[u]++] = v;
                            }
                        }
                    }
                }
            }

            if (iter_pair.first->second.empty()) {
                frontiers[j] = INVALID_VERTEX_ID;
                for (ui k = 0; k < tree[u_p].children_count_; ++k) {
                    VertexID u_c = tree[u_p].children_[k];
                    if (visited_query_vertex[u_c]) {
                        TE_Candidates[u_c].erase(v_f);
                    }
                }
            }
        }

        if (candidates_count[u] == 0)
            return false;

        for (ui j = 0; j < candidates_count[u]; ++j) {
            VertexID v = candidates[u][j];
            flag[v] = 0;
        }
    }

    
    NTE_Candidates.resize(query_vertices_count);
    for (auto& value : NTE_Candidates) {
        value.resize(query_vertices_count);
    }

    for (ui i = 1; i < query_vertices_count; ++i) {
        VertexID u = order[i];
        TreeNode &u_node = tree[u];

        ui u_label = query_graph->getVertexLabel(u);
        ui u_degree = query_graph->getVertexDegree(u);
        const std::unordered_map<LabelID, ui> *u_nlf = query_graph->getVertexNLF(u);
        for (ui l = 0; l < u_node.bn_count_; ++l) {
            VertexID u_p = u_node.bn_[l];
            VertexID *frontiers = candidates[u_p];
            ui frontiers_count = candidates_count[u_p];

            for (ui j = 0; j < frontiers_count; ++j) {
                VertexID v_f = frontiers[j];

                if (v_f == INVALID_VERTEX_ID)
                    continue;

                ui nbrs_cnt;
                const VertexID *nbrs = data_graph->getVertexNeighbors(v_f, nbrs_cnt);

                auto iter_pair = NTE_Candidates[u][u_p].emplace(v_f, std::vector<VertexID>());
                for (ui k = 0; k < nbrs_cnt; ++k) {
                    VertexID v = nbrs[k];

                    if (data_graph->getVertexLabel(v) == u_label && data_graph->getVertexDegree(v) >= u_degree) {
                        const std::unordered_map<LabelID, ui> *v_nlf = data_graph->getVertexNLF(v);

                        if (v_nlf->size() >= u_nlf->size()) {
                            bool is_valid = true;

                            for (auto element : *u_nlf) {
                                auto iter = v_nlf->find(element.first);
                                if (iter == v_nlf->end() || iter->second < element.second) {
                                    is_valid = false;
                                    break;
                                }
                            }

                            if (is_valid) {
                                iter_pair.first->second.push_back(v);
                            }
                        }
                    }
                }

                if (iter_pair.first->second.empty()) {
                    frontiers[j] = INVALID_VERTEX_ID;
                    for (ui k = 0; k < tree[u_p].children_count_; ++k) {
                        VertexID u_c = tree[u_p].children_[k];
                        TE_Candidates[u_c].erase(v_f);
                    }
                }
            }
        }
    }

    
    std::vector<std::vector<ui>> cardinality(query_vertices_count);
    for (ui i = 0; i < query_vertices_count; ++i) {
        cardinality[i].resize(candidates_count[i], 1);
    }

    std::vector<ui> local_cardinality(data_vertices_count);
    std::fill(local_cardinality.begin(), local_cardinality.end(), 0);

    for (int i = query_vertices_count - 1; i >= 0; --i) {
        VertexID u = order[i];
        TreeNode& u_node = tree[u];

        ui flag_num = 0;
        ui updated_flag_count = 0;

        
        for (ui j = 0; j < candidates_count[u]; ++j) {
            VertexID v = candidates[u][j];

            if (v == INVALID_VERTEX_ID)
                continue;

            if (flag[v] == flag_num) {
                flag[v] += 1;
                updated_flag[updated_flag_count++] = v;
            }
        }

        for (ui j = 0; j < u_node.bn_count_; ++j) {
            VertexID u_bn = u_node.bn_[j];
            flag_num += 1;
            for (auto iter = NTE_Candidates[u][u_bn].begin(); iter != NTE_Candidates[u][u_bn].end(); ++iter) {
                for (auto v : iter->second) {
                    if (flag[v] == flag_num) {
                        flag[v] += 1;
                    }
                }
            }
        }

        flag_num += 1;

        
        for (ui j = 0; j < candidates_count[u]; ++j) {
            VertexID v = candidates[u][j];
            if (v != INVALID_VERTEX_ID && flag[v] == flag_num) {
                local_cardinality[v] = cardinality[u][j];
            }
            else {
                cardinality[u][j] = 0;
            }
        }

        VertexID u_p = u_node.parent_;
        VertexID* frontiers = candidates[u_p];
        ui frontiers_count = candidates_count[u_p];

        
        for (ui j = 0; j < frontiers_count; ++j) {
            VertexID v_f = frontiers[j];

            if (v_f == INVALID_VERTEX_ID) {
                cardinality[u_p][j] = 0;
                continue;
            }

            ui temp_score = 0;
            for (auto iter = TE_Candidates[u][v_f].begin(); iter != TE_Candidates[u][v_f].end();) {
                VertexID v = *iter;
                temp_score += local_cardinality[v];
                if (local_cardinality[v] == 0) {
                    iter = TE_Candidates[u][v_f].erase(iter);
                    for (ui k = 0; k < u_node.children_count_; ++k) {
                        VertexID u_c = u_node.children_[k];
                        TE_Candidates[u_c].erase(v);
                    }

                    for (ui k = 0; k < u_node.fn_count_; ++k) {
                        VertexID u_c = u_node.fn_[k];
                        NTE_Candidates[u_c][u].erase(v);
                    }
                }
                else {
                    ++iter;
                }
            }

            cardinality[u_p][j] *= temp_score;
        }

        
        for (ui j = 0; j < updated_flag_count; ++j) {
            flag[updated_flag[j]] = 0;
            local_cardinality[updated_flag[j]] = 0;
        }
    }

    compactCandidates(candidates, candidates_count, query_vertices_count);
    sortCandidates(candidates, candidates_count, query_vertices_count);


    for (ui i = 0; i < query_vertices_count; ++i) {
        if (candidates_count[i] == 0) {
            return false;
        }
    }

    for (ui i = 1; i < query_vertices_count; ++i) {
        VertexID u = order[i];
        TreeNode& u_node = tree[u];

        
        {
            VertexID u_p = u_node.parent_;
            auto iter = TE_Candidates[u].begin();
            while (iter != TE_Candidates[u].end()) {
                VertexID v_f = iter->first;
                if (!std::binary_search(candidates[u_p], candidates[u_p] + candidates_count[u_p], v_f)) {
                    iter = TE_Candidates[u].erase(iter);
                }
                else {
                    std::sort(iter->second.begin(), iter->second.end());
                    iter++;
                }
            }
        }

        
        {
            for (ui j = 0; j < u_node.bn_count_; ++j) {
                VertexID u_p = u_node.bn_[j];
                auto iter = NTE_Candidates[u][u_p].end();
                while (iter != NTE_Candidates[u][u_p].end()) {
                    VertexID v_f = iter->first;
                    if (!std::binary_search(candidates[u_p], candidates[u_p] + candidates_count[u_p], v_f)) {
                        iter = NTE_Candidates[u][u_p].erase(iter);
                    }
                    else {
                        std::sort(iter->second.begin(), iter->second.end());
                        iter++;
                    }
                }
            }
        }
    }
    return true;
}

void
InterPartitionFilter::computeCandidateWithNLF(const Graph *data_graph, const Graph *query_graph, VertexID query_vertex,
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

void
InterPartitionFilter::computeRootCandidateWithNLF(const Graph *data_graph, const Graph *query_graph, VertexID query_vertex, ui &count,
                                                    std::unordered_map<uint32_t, std::vector<VertexID> >& getReplicabyLabel, ui *buffer) {
    LabelID label = query_graph->getVertexLabel(query_vertex);
    ui degree = query_graph->getVertexDegree(query_vertex);
    const std::unordered_map<LabelID, ui>* query_vertex_nlf = query_graph->getVertexNLF(query_vertex);
    std::vector<VertexID> data_vertices = getReplicabyLabel[label];
    count = 0;
    for (ui j = 0; j < data_vertices.size(); ++j) {
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

void InterPartitionFilter::sortCandidates(ui **candidates, ui *candidates_count, ui num) {
    for (ui i = 0; i < num; ++i) {
        std::sort(candidates[i], candidates[i] + candidates_count[i]);
    }
}

void InterPartitionFilter::allocateBuffer(const Graph *data_graph, const Graph *query_graph, ui **&candidates,
                                    ui *&candidates_count) {
    ui query_vertices_num = query_graph->getVerticesCount();
    ui candidates_max_num = data_graph->getGraphMaxLabelFrequency();

    candidates_count = new ui[query_vertices_num];
    memset(candidates_count, 0, sizeof(ui) * query_vertices_num);

    candidates = new ui*[query_vertices_num];

    for (ui i = 0; i < query_vertices_num; ++i) {
        candidates[i] = new ui[candidates_max_num];
    }
}

void InterPartitionFilter::compactCandidates(ui **&candidates, ui *&candidates_count, ui query_vertices_num) {
    for (ui i = 0; i < query_vertices_num; ++i) {
        VertexID query_vertex = i;
        ui next_position = 0;
        for (ui j = 0; j < candidates_count[query_vertex]; ++j) {
            VertexID data_vertex = candidates[query_vertex][j];

            if (data_vertex != INVALID_VERTEX_ID) {
                candidates[query_vertex][next_position++] = data_vertex;
            }
        }

        candidates_count[query_vertex] = next_position;
    }
}

bool InterPartitionFilter::isCandidateSetValid(ui **&candidates, ui *&candidates_count, ui query_vertex_num) {
    for (ui i = 0; i < query_vertex_num; ++i) {
        if (candidates_count[i] == 0)
            return false;
    }
    return true;
}

void GenerateFilteringPlan::generateFilterPlan(const Graph *data_graph, const Graph *query_graph, TreeNode *&tree,
                                                   VertexID *&order, VertexID root) {
    VertexID start_vertex = root;
    GraphOperations::bfsTraversal(query_graph, start_vertex, tree, order);

    ui query_vertices_num = query_graph->getVerticesCount();
    std::vector<ui> order_index(query_vertices_num);
    for (ui i = 0; i < query_vertices_num; ++i) {
        VertexID query_vertex = order[i];
        order_index[query_vertex] = i;
    }

    for (ui i = 0; i < query_vertices_num; ++i) {
        VertexID u = order[i];
        tree[u].under_level_count_ = 0;
        tree[u].bn_count_ = 0;
        tree[u].fn_count_ = 0;

        ui u_nbrs_count;
        const VertexID* u_nbrs = query_graph->getVertexNeighbors(u, u_nbrs_count);
        for (ui j = 0; j < u_nbrs_count; ++j) {
            VertexID u_nbr = u_nbrs[j];
            if (u_nbr != tree[u].parent_ && order_index[u_nbr] < order_index[u]) {
                tree[u].bn_[tree[u].bn_count_++] = u_nbr;
                tree[u_nbr].fn_[tree[u_nbr].fn_count_++] = u;
            }
        }
    }
}

VertexID GenerateFilteringPlan::selectStartVertex(const Graph *data_graph, const Graph *query_graph) {
    double min_score = data_graph->getVerticesCount();
    VertexID start_vertex = 0;

    for (ui i = 0; i < query_graph->getVerticesCount(); ++i) {
        ui degree = query_graph->getVertexDegree(i);
        ui count = 0;
        InterPartitionFilter::computeCandidateWithNLF(data_graph, query_graph, i, count);
        double cur_score = count / (double)degree;
        if (cur_score < min_score) {
            min_score = cur_score;
            start_vertex = i;
        }
    }

    return start_vertex;
}

void InterPartitionQueryPlan::generateQueryPlan(const Graph *query_graph, TreeNode *tree, ui *bfs_order, ui *&order,
                                              ui *&pivot) {
    ui query_vertices_num = query_graph->getVerticesCount();
    order = new ui[query_vertices_num];
    pivot = new ui[query_vertices_num];

    for (ui i = 0; i < query_vertices_num; ++i) {
        order[i] = bfs_order[i];
    }

    for (ui i = 1; i < query_vertices_num; ++i) {
        pivot[i] = tree[order[i]].parent_;
    }
}

size_t
InterPartitionEnumeration::enumerate(const Graph *data_graph, const Graph *query_graph, TreeNode *tree, ui **candidates,
                                ui *candidates_count,
                                std::vector<std::unordered_map<VertexID, std::vector<VertexID>>> &TE_Candidates,
                                std::vector<std::vector<std::unordered_map<VertexID, std::vector<VertexID>>>> &NTE_Candidates,
                                ui *order, size_t &output_limit_num, size_t &call_count,
                                const std::unordered_map<VertexID, std::pair<std::set<VertexID>, std::set<VertexID>>>& ordered_constraints,
                                pair<VertexID, VertexID>* ret_embedding) {

    ui max_depth = query_graph->getVerticesCount();
    ui data_vertices_count = data_graph->getVerticesCount();
    ui max_valid_candidates_count = 0;
    for (ui i = 0; i < max_depth; ++i) {
        if (candidates_count[i] > max_valid_candidates_count) {
            max_valid_candidates_count = candidates_count[i];
        }
    }
    
    ui *idx = new ui[max_depth];
    ui *idx_count = new ui[max_depth];
    ui *embedding = new ui[max_depth];
    ui *temp_buffer = new ui[max_valid_candidates_count];
    ui **valid_candidates = new ui *[max_depth];
    for (ui i = 0; i < max_depth; ++i) {
        valid_candidates[i] = new ui[max_valid_candidates_count];
    }
    bool *visited_vertices = new bool[data_vertices_count];
    std::fill(visited_vertices, visited_vertices + data_vertices_count, false);

    
    size_t embedding_cnt = 0;
    int cur_depth = 0;
    VertexID start_vertex = order[0];

    idx[cur_depth] = 0;
    idx_count[cur_depth] = candidates_count[start_vertex];

    for (ui i = 0; i < idx_count[cur_depth]; ++i) {
        valid_candidates[cur_depth][i] = candidates[start_vertex][i];
    }

    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> ancestors;
    computeAncestor(query_graph, order, ancestors);
    AncestorWithSymmetryBreaking(query_graph, ancestors, ordered_constraints);
    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> vec_failing_set(max_depth);
    std::unordered_map<VertexID, VertexID> reverse_embedding;
    reverse_embedding.reserve(MAXIMUM_QUERY_GRAPH_SIZE * 2);

    while (true) {
        while (idx[cur_depth] < idx_count[cur_depth]) {
            VertexID u = order[cur_depth];
            VertexID v = valid_candidates[cur_depth][idx[cur_depth]];
            idx[cur_depth] += 1;

            if (visited_vertices[v]) {
                vec_failing_set[cur_depth] = ancestors[u];
                vec_failing_set[cur_depth] |= ancestors[reverse_embedding[v]];
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                continue;
            }

            embedding[u] = v;
            visited_vertices[v] = true;

            reverse_embedding[v] = u;
            if (cur_depth == max_depth - 1) {
#ifdef OUTPUT_RESULT
                if (embedding_cnt < output_limit_num) {
                    for (int i = 0; i < max_depth; ++i) {
                        VertexID query_vertex = order[i];
                        VertexID data_vertex = embedding[i];
                        ret_embedding[embedding_cnt * max_depth + i] = make_pair(query_vertex, data_vertex);
                    }
                }
#endif
                embedding_cnt += 1;
                visited_vertices[v] = false;
                reverse_embedding.erase(embedding[u]);
                vec_failing_set[cur_depth].set();
                vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                if (embedding_cnt >= output_limit_num) {
                    goto EXIT;
                }
            } else {
                call_count += 1;
                cur_depth += 1;
                idx[cur_depth] = 0;
                generateValidCandidates(cur_depth, embedding, idx_count, valid_candidates, order, temp_buffer, tree,
                                        TE_Candidates,
                                        NTE_Candidates);
                pruneCandidatesBySymmetryBreaking(cur_depth, embedding, order,
                                                  idx_count, valid_candidates, ordered_constraints);
                if (idx_count[cur_depth] == 0) {
                    vec_failing_set[cur_depth - 1] = ancestors[order[cur_depth]];
                } else {
                    vec_failing_set[cur_depth - 1].reset();
                }
            }
        }

        cur_depth -= 1;
        if (cur_depth < 0)
            break;
        else {
            VertexID u = order[cur_depth];
            reverse_embedding.erase(embedding[u]);
            if (cur_depth != 0) {
                if (!vec_failing_set[cur_depth].test(u)) {
                    vec_failing_set[cur_depth - 1] = vec_failing_set[cur_depth];
                    idx[cur_depth] = idx_count[cur_depth];
                } else {
                    vec_failing_set[cur_depth - 1] |= vec_failing_set[cur_depth];
                }
            }
            visited_vertices[embedding[u]] = false;
        }
    }

    
    EXIT:
    delete[] idx;
    delete[] idx_count;
    delete[] embedding;
    delete[] temp_buffer;
    delete[] visited_vertices;
    for (ui i = 0; i < max_depth; ++i) {
        delete[] valid_candidates[i];
    }
    delete[] valid_candidates;

    return embedding_cnt;
}

void InterPartitionEnumeration::generateValidCandidates(ui depth, ui *embedding, ui *idx_count, ui **valid_candidates, ui *order,
                                            ui *&temp_buffer, TreeNode *tree,
                                            std::vector<std::unordered_map<VertexID, std::vector<VertexID>>> &TE_Candidates,
                                            std::vector<std::vector<std::unordered_map<VertexID, std::vector<VertexID>>>> &NTE_Candidates) {

    VertexID u = order[depth];
    idx_count[depth] = 0;
    ui valid_candidates_count = 0;
    {
        VertexID u_p = tree[u].parent_;
        VertexID v_p = embedding[u_p];

        auto iter = TE_Candidates[u].find(v_p);
        if (iter == TE_Candidates[u].end() || iter->second.empty()) {
            return;
        }

        valid_candidates_count = iter->second.size();
        VertexID *v_p_nbrs = iter->second.data();

        for (ui i = 0; i < valid_candidates_count; ++i) {
            valid_candidates[depth][i] = v_p_nbrs[i];
        }
    }
    ui temp_count;
    for (ui i = 0; i < tree[u].bn_count_; ++i) {
        VertexID u_p = tree[u].bn_[i];
        VertexID v_p = embedding[u_p];

        auto iter = NTE_Candidates[u][u_p].find(v_p);
        if (iter == NTE_Candidates[u][u_p].end() || iter->second.empty()) {
            return;
        }

        ui current_candidates_count = iter->second.size();
        ui *current_candidates = iter->second.data();

        ComputeSetIntersection::ComputeCandidates(current_candidates, current_candidates_count,
                                                  valid_candidates[depth], valid_candidates_count,
                                                  temp_buffer, temp_count);

        std::swap(temp_buffer, valid_candidates[depth]);
        valid_candidates_count = temp_count;
    }

    idx_count[depth] = valid_candidates_count;
}

void 
InterPartitionEnumeration::pruneCandidatesBySymmetryBreaking(ui depth, ui *embedding, ui* order,
                                                  ui *idx_count, ui **valid_candidates,
                                                  const std::unordered_map<VertexID, std::pair<std::set<VertexID>, std::set<VertexID>>>& ordered_constraints) {
    VertexID u = order[depth];

    auto it = ordered_constraints.find(u);
    if (it != ordered_constraints.end()) {
        
        VertexID min_id = 0;
        VertexID max_id = -1;

        const auto& u_cons_prior = it->second.first;
        const auto& u_cons_inferior = it->second.second;

        for (auto& u_p: u_cons_prior) {
            if (embedding[u_p] > min_id) min_id = embedding[u_p]; 
        }

        for (auto& u_i: u_cons_inferior) {
            if (embedding[u_i] < max_id) max_id = embedding[u_i]; 
        }

        
        int valid_cnt = 0;
        for (int i = 0; i < idx_count[depth]; i++) {
            VertexID v = valid_candidates[depth][i];
            if (min_id <= v && v < max_id) {
                valid_candidates[depth][valid_cnt++] = v;
            }
        }
        idx_count[depth] = valid_cnt;
    }
}

void InterPartitionEnumeration::computeAncestor(const Graph *query_graph, VertexID *order,
                                    std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> &ancestors) {
    ui query_vertices_num = query_graph->getVerticesCount();
    ancestors.resize(query_vertices_num);

    
    for (ui i = 0; i < query_vertices_num; ++i) {
        VertexID u = order[i];
        ancestors[u].set(u);
        for (ui j = 0; j < i; ++j) {
            VertexID u_bn = order[j];
            if (query_graph->checkEdgeExistence(u, u_bn)) {
                ancestors[u] |= ancestors[u_bn];
            }
        }
    }
}

void 
InterPartitionEnumeration::AncestorWithSymmetryBreaking(const Graph *query_graph,
                             std::vector<std::bitset<MAXIMUM_QUERY_GRAPH_SIZE>> &ancestors,
                             const std::unordered_map<VertexID, std::pair<std::set<VertexID>, std::set<VertexID>>>& ordered_constraints) {
    ui query_vertices_num = query_graph->getVerticesCount();

    for (auto it = ordered_constraints.begin(); it != ordered_constraints.end(); it++) {
        VertexID u = it->first;
        const auto& u_cons_prior = it->second.first;
        const auto& u_cons_inferior = it->second.second;

        
        
        for (auto& u_p: u_cons_prior) {      
            ancestors[u] |= ancestors[u_p];
        }

        for (auto& u_i: u_cons_inferior) {
            ancestors[u] |= ancestors[u_i];
        }
    }
}