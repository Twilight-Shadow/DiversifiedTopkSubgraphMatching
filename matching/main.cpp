#include <chrono>
#include <future>
#include <thread>
#include <fstream>
#include <random>
#include <dirent.h>
#include <sys/types.h>
#include <cstdlib>
#include <cstdio>
#include <unistd.h>
#include <vector>
#include <map>
#include <set>
#include <atomic>
#include <mutex>
#include <condition_variable>
#include <unordered_set>
#include <queue>
#include <pthread.h>
#include <semaphore.h>
#include <random>
#include <ctime>

#include "CommandList.h"
#include "graph/MatchGraph.h"
#include "IntraMatch.h"
#include "InterMatch.h"
#include "BuildTable.h"
#include "ctpl/ctpl_stl.h"
#include "configuration/Config.h"
#include "utility/analyze_symmetry.h"

using std::pair;
using std::make_pair;

uint32_t embedding_cnt;
bool _terminate_;
std::mutex emb_lock;
bool* inter_terminate_flag;

uint32_t searchIntraPartition(int id, Graph* query_graph, Graph* partition_graph, pair<VertexID, VertexID>* ret_embedding, 
                            uint32_t partition_index, std::unordered_map<VertexID, VertexID>& part2data) {
    if (_terminate_) return 0;
    int query_size = query_graph->getVerticesCount();
    
    // filter & sort candidates
    ui** candidates = NULL;
    ui* candidates_count = NULL;
    IntraPartitionFilter::Filter(partition_graph, query_graph, candidates, candidates_count);
    IntraPartitionFilter::sortCandidates(candidates, candidates_count, query_size);

    // build edge matrix for candidates
    Edges ***edge_matrix = NULL;
    edge_matrix = new Edges **[query_size];
    for (int i = 0; i < query_size; ++i) {
        edge_matrix[i] = new Edges *[query_size];
    }
    BuildTable::buildTables(partition_graph, query_graph, candidates, candidates_count, edge_matrix);

    // generate query plan
    ui* matching_order = NULL;
    ui* pivots = NULL;
    IntraPartitionQueryPlan::generateQueryPlan(partition_graph, query_graph, candidates_count, matching_order, pivots);

    // intra-partition matching
    size_t output_limit = 1, call_count = 0;
    uint32_t intra_embedding_cnt = IntraPartitionEnumeration::LFTJ(partition_graph, query_graph, edge_matrix, candidates, candidates_count,
                                                                    matching_order, output_limit, call_count, ret_embedding);
    
    if (intra_embedding_cnt >= 1) {
        for (int i = 0; i < query_size; ++i) {
            ret_embedding[i].second = part2data[ret_embedding[i].second];
        }
    }

    // release intra matching memory                                                
    delete[] candidates_count;
    delete[] matching_order;
    delete[] pivots;
    for (int i = 0; i < query_size; ++i) {
        delete[] candidates[i];
    }
    delete[] candidates;
    if (edge_matrix != NULL) {
        for (int i = 0; i < query_size; ++i) {
            for (int j = 0; j < query_size; ++j) {
                delete edge_matrix[i][j];
            }
            delete[] edge_matrix[i];
        }
        delete[] edge_matrix;
    }
    part2data.clear();
    return intra_embedding_cnt;
}

pair<uint32_t, pair<VertexID, VertexID>*> 
searchInterPartition(int id, VertexID root_ID, Graph* query_graph, Graph* data_graph, uint32_t partition_index,
                        std::unordered_map<uint32_t, std::vector<VertexID> >& getReplicabyLabel,
                        std::vector<std::set<std::pair<VertexID, VertexID>>>& permutations,
                        std::vector<std::pair<VertexID, VertexID>>& constraints,
                        std::unordered_map<VertexID, std::pair<std::set<VertexID>, std::set<VertexID>>>& full_constraints) {
    if (_terminate_ || inter_terminate_flag[partition_index]) {
        return make_pair(0, (pair<VertexID, VertexID>*)NULL);
    }

    int query_size = query_graph->getVerticesCount();
    ui** candidate = NULL;
    ui* candidate_count = NULL;
    TreeNode* tree = NULL;
    ui* order = NULL;
    std::vector<std::unordered_map<VertexID, std::vector<VertexID>>> TE_Candidates;
    std::vector<std::vector<std::unordered_map<VertexID, std::vector<VertexID>>>> NTE_Candidates;
    InterPartitionFilter::Filter(data_graph, query_graph, candidate, candidate_count, order, tree,
                                    TE_Candidates, NTE_Candidates, root_ID, getReplicabyLabel);

    ui* matching_order = NULL;
    ui* pivots = NULL;
    InterPartitionQueryPlan::generateQueryPlan(query_graph, tree, order, matching_order, pivots);

    std::unordered_map<VertexID, std::pair<std::set<VertexID>, std::set<VertexID>>> ordered_constraints;
    ANALYZE_SYMMETRY::make_ordered_constraints(matching_order, query_size, full_constraints, ordered_constraints);

    size_t output_limit = 1, call_count = 0;
    pair<VertexID, VertexID>* try_match = new pair<VertexID, VertexID>[query_size];
    uint32_t try_match_cnt = 0;
    for (int i = 0; i < query_size; ++i) {
        if (candidate_count[i] == 0) goto EXIT_INTER;
    }
    try_match_cnt = InterPartitionEnumeration::enumerate(
        data_graph, query_graph, tree, candidate, candidate_count, TE_Candidates, NTE_Candidates,
        order, output_limit, call_count, ordered_constraints, try_match
    );

EXIT_INTER:
    delete[] candidate_count;
    delete[] tree;
    delete[] order;
    delete[] matching_order;
    delete[] pivots;
    TE_Candidates.clear();
    NTE_Candidates.clear();
    ordered_constraints.clear();
    for (int i = 0; i < query_size; ++i) {
        delete[] candidate[i];
    }
    delete[] candidate;
    return make_pair(try_match_cnt, try_match);
}

void ParallelSubgraphMatching(ctpl::thread_pool& pool, uint32_t k, Graph* query_graph, Graph* data_graph, Graph* partition_graph, int partition_index, 
                                std::unordered_map<uint32_t, std::vector<VertexID>>& getReplicabyLabel, 
                                std::unordered_map<VertexID, VertexID>& part2data, pair<VertexID, VertexID>* embedding, 
                                std::vector<std::set<std::pair<VertexID, VertexID>>>& permutations, std::vector<std::pair<VertexID, VertexID>>& constraints,
                                std::unordered_map<VertexID, std::pair<std::set<VertexID>, std::set<VertexID>>>& full_constraints) {
    
    if (_terminate_ || inter_terminate_flag[partition_index]) return;
    // load data
    int query_size = query_graph->getVerticesCount();
    uint32_t intra_embedding_cnt = 0, inter_embedding_cnt = 0;
    pair<VertexID, VertexID>* ret_embedding = new pair<VertexID, VertexID>[query_size];

    // intra-partition matching
    std::future<uint32_t> intra_future = pool.push(
        searchIntraPartition, query_graph, partition_graph, ret_embedding, partition_index, std::ref(part2data)
    );
    intra_embedding_cnt = intra_future.get();

    // inter-partition matching
    if (intra_embedding_cnt <= 0) {
        if (_terminate_ || inter_terminate_flag[partition_index]) goto EXIT_MATCHING;
        std::mutex part_lock;
        std::vector<std::future<pair<uint32_t, pair<VertexID, VertexID>*>>> inter_futures;

        // inter-partition matching
        for (int i = 0; i < query_size; ++i) {
            inter_futures.push_back(
                pool.push(
                    searchInterPartition,
                    i,
                    query_graph,
                    data_graph,
                    partition_index,
                    std::ref(getReplicabyLabel),
                    std::ref(permutations),
                    std::ref(constraints),
                    std::ref(full_constraints)
                )
            );
        }
        for (auto& future : inter_futures) {
            auto result = future.get();
            if (result.first > 0) {
                inter_embedding_cnt = result.first;
                delete[] ret_embedding;
                ret_embedding = result.second;
                inter_terminate_flag[partition_index] = 1;
                break;
            }
        }
    }

EXIT_MATCHING:
    emb_lock.lock();
    if (!_terminate_ && (intra_embedding_cnt || inter_embedding_cnt)){
#ifdef OUTPUT_RESULT
        for (int i = 0; i < query_size; ++i) {
            embedding[embedding_cnt * query_size + i] = ret_embedding[i];
        }
        delete[] ret_embedding;
#endif
        embedding_cnt += 1;
        if (embedding_cnt >= k) _terminate_ = true;
    }
    emb_lock.unlock();
    return;
}

template <typename T>
void LoadMatrix(std::string input_data_path, int n, int m, T** mat) {
    std::ifstream input_file(input_data_path, std::ios::in);
    if (!input_file.is_open()) {
        std::cerr << input_data_path + ": no such file or directory.\n";
        return;
    }
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            input_file >> mat[i][j];
        }
    }
    input_file.close();
    std::cout << "Finish load file " + input_data_path << ".\n";
    return;
}

void LoadNodeReplica(std::string input_data_path, std::vector<VertexID>& replica) {
    std::ifstream input_file(input_data_path, std::ios::in);
    replica.clear();
    VertexID node = 0;
    while (input_file >> node) {
        replica.push_back(node);
    }
    return;
}

void SelectKFarthestPartition(uint32_t k, std::vector<VertexID>& candidates, std::vector<VertexID>& selected, \
                                uint32_t** distance_matrix, bool* alloc_bitmap, uint32_t n_partition) {
    std::vector<uint32_t> min_dist(n_partition, 0x7fffffff);
    uint32_t init = 0;
    if (candidates.size() == 0) {
        uint32_t root = rand() % n_partition;
        candidates.push_back(root);
        selected.push_back(root);
        alloc_bitmap[root] = 1;
        init = 1;
    }
    for (int i = 0; i < n_partition; ++i) {
        if (alloc_bitmap[i] == 1) {
            min_dist[i] = 0;
            continue;
        }
        for (auto j = candidates.begin(); j != candidates.end(); ++j) {
            min_dist[i] = std::min(min_dist[i], distance_matrix[i][*j]);
        }
    }
    for (int i = init; i < k; ++i) {
        uint32_t max_dist = 0;
        int choice = -1;
        for (int j = 0; j < n_partition; ++j) {
            if (alloc_bitmap[j]) continue;
            if (min_dist[j] > max_dist) {
                max_dist = min_dist[j];
                choice = j;
            }
        }
        if (choice == -1) break;
        candidates.push_back(choice);
        selected.push_back(choice);
        alloc_bitmap[choice] = 1;
        // update
        for (int j = 0; j < n_partition; ++j) {
            if (!alloc_bitmap[j]) {
                min_dist[j] = std::min(min_dist[j], distance_matrix[choice][j]);
            }
        }
    }
    // if (selected.size() > 0) {
    //     std::cout << "Finish selecting " << selected.size() << " partitions: \n";
    //     for (auto i = selected.begin(); i != selected.end(); ++i) {
    //         std::cout << *i << ' ';
    //     }
    //     std::cout << '\n';
    // }
    return;
}

int main(int argc, char** argv){
    MatchingCommand command(argc, argv);
    std::string query_graph_path = command.getQueryGraphFilePath();
    std::string partitioned_data_graph_path = command.getPartitionedDataGraphDirectoryPath();
    std::string replicated_node_path = command.getReplicatedNodeDirectory();
    std::string data_graph_path = command.getDataGraphFilePath();
    std::string distance_matrix_path = command.getDistanceFilePath();
    std::string groundtruth_path = command.getGroundTruthFilePath();
    uint32_t n_thread = std::stoi(command.getThreadCount());
    uint32_t k = std::stoi(command.getK());
    int n_partition = std::stoi(command.getPartitionNum());

    std::cout << "Query Graph: " << query_graph_path << '\n';
    std::cout << "Data Graph: " << data_graph_path << '\n';
    std::cout << "Parititoned Data Graph: " << partitioned_data_graph_path << '\n';
    std::cout << "Replicated Nodes Info: " << replicated_node_path << '\n';
    std::cout << "Distance Matrix: " << distance_matrix_path << '\n';
    std::cout << "Ground Truth: " << groundtruth_path << '\n';
    std::cout << "Number of Partition: " << n_partition << '\n';
    std::cout << "Top-k value: " << k << '\n';
    std::cout << "Parallel Thread Number: " << n_thread << '\n';

    ctpl::thread_pool pool(n_thread);
    std::vector<std::thread> partition_thread;
    std::srand(static_cast<unsigned int>(std::time(nullptr)));
    embedding_cnt = 0;
    _terminate_ = false;

    inter_terminate_flag = new bool[n_partition];
    for (int i = 0; i < n_partition; ++i) {
        inter_terminate_flag[i] = false;
    }

    // Load Data
    uint32_t** distance_matrix = new uint32_t*[n_partition];
    for (int i = 0; i < n_partition; ++i) {
        distance_matrix[i] = new uint32_t[n_partition];
    }
    LoadMatrix(distance_matrix_path, n_partition, n_partition, distance_matrix);

    // load query graph
    Graph* query_graph = new Graph(true);
    query_graph->loadGraphFromFile(query_graph_path);
    query_graph->buildCoreTable();
    query_graph->printGraphMetaData();

    // load complete data graph
    Graph* data_graph = new Graph(true);
    data_graph->loadGraphFromFile(data_graph_path);
    data_graph->printGraphMetaData();

    // load partition graph
    Graph** partition_graphs = new Graph*[n_partition];
    std::vector<std::unordered_map<VertexID, VertexID>> part2data(n_partition);
    for (int i = 0; i < n_partition; ++i) {
        partition_graphs[i] = new Graph(true);
        partition_graphs[i]->loadGraphFromFile(partitioned_data_graph_path + std::to_string(i) + ".graph", part2data[i]);
    }

    // load node replica
    std::vector<std::unordered_map<uint32_t, std::vector<VertexID>>> replica_label_group(n_partition);
    for (int i = 0; i < n_partition; ++i) {
        std::vector<VertexID> replica;
        LoadNodeReplica(replicated_node_path + std::to_string(i) + ".txt", replica);
        for (VertexID j : replica) {
            replica_label_group[i][data_graph->getVertexLabel(j)].push_back(j);
        }
        replica.clear();
    }

    std::vector<std::set<std::pair<VertexID, VertexID>>> permutations;
    std::vector<std::pair<VertexID, VertexID>> constraints;
    std::unordered_map<VertexID, std::pair<std::set<VertexID>, std::set<VertexID>>> full_constraints;
    std::unordered_map<VertexID, std::set<VertexID>> cosets = ANALYZE_SYMMETRY::analyze_symmetry(query_graph, permutations);
    ANALYZE_SYMMETRY::make_constraints(cosets, constraints);
    ANALYZE_SYMMETRY::make_full_constraints(constraints, full_constraints);
    pair<VertexID, VertexID>* embedding = new pair<VertexID, VertexID>[query_graph->getVerticesCount() * k];

#ifdef PDD          // run algorithm PDD, switch in configuration/Config.h
    bool* alloc_bitmap = new bool[n_partition];
    for (int i = 0; i < n_partition; ++i) {
        alloc_bitmap[i] = false;
    }
    std::vector<uint32_t> candidates, selected;

    auto start = std::chrono::high_resolution_clock::now();

BACKTRACKING:
    selected.clear();
    SelectKFarthestPartition(k, candidates, selected, distance_matrix, alloc_bitmap, n_partition);
    if (selected.size() == 0) goto EXIT_PDD;
    for (auto i = selected.begin(); i != selected.end(); ++i) {
        partition_thread.emplace_back(
            ParallelSubgraphMatching, 
            std::ref(pool),
            k,
            query_graph, data_graph,
            partition_graphs[*i],
            *i,
            std::ref(replica_label_group[*i]),
            std::ref(part2data[*i]),
            embedding,
            std::ref(permutations),
            std::ref(constraints),
            std::ref(full_constraints)
        );
    }
    for (auto& t : partition_thread) {
        if (t.joinable()) t.join();
    }
    partition_thread.clear();
    if (!_terminate_) goto BACKTRACKING;
EXIT_PDD:

    auto end = std::chrono::high_resolution_clock::now();

    candidates.clear();
    selected.clear();
    delete[] alloc_bitmap;
#endif

#ifdef PDD_PLUS     // run algorithm PDD+, switch in configuration/Config.h
    // Load extra data
    bool** ground_truth = new bool*[1];
    ground_truth[0] = new bool[n_partition];
    LoadMatrix(groundtruth_path, 1, n_partition, ground_truth);

    // build PDG
    std::vector<pair<uint32_t, uint32_t> > PDG[n_partition];
    std::vector<uint32_t> available_partition;
    ui* PDG_deg = new ui[n_partition];
    for (int i = 0; i < n_partition; ++i) {
        PDG_deg[i] = 0;
    }
    for (int i = 0; i < n_partition; ++i) {
        for (int j = 0; j < n_partition; ++j) {
            if (distance_matrix[i][j] <= 1 || ground_truth[0][i] == 0 || ground_truth[0][j] == 0) continue;
            PDG[i].push_back(make_pair(j, distance_matrix[i][j]));
            PDG[j].push_back(make_pair(i, distance_matrix[i][j]));
            ++PDG_deg[i];
            ++PDG_deg[j];
        }
        if (ground_truth[0][i]) available_partition.push_back(i);
    }
    std::cout << "Finish building PDG.\n";
    std::cout << "Available Partition Number: " << available_partition.size() << '\n';

    // select candidates in positive partitions
    bool* alloc_bitmap = new bool[n_partition];
    for (int i = 0; i < n_partition; ++i) {
        alloc_bitmap[i] = false;
    }

    auto start = std::chrono::high_resolution_clock::now();

    VertexID isolated_part = 0;
    uint32_t iso_max_deg = 0;
    for (auto i = available_partition.begin(); i != available_partition.end(); ++i) {
        if (PDG_deg[*i] > iso_max_deg) {
            isolated_part = *i;
            iso_max_deg = PDG_deg[*i];
        }
    }
    alloc_bitmap[isolated_part] = 1;
    std::vector<uint32_t> min_dist(n_partition, 0x7fffffff);
    for (int i = 0; i < n_partition; ++i) {
        min_dist[i] = distance_matrix[isolated_part][i];
    }
    partition_thread.emplace_back(
        ParallelSubgraphMatching, 
        std::ref(pool),
        k,
        query_graph, data_graph,
        partition_graphs[isolated_part],
        isolated_part,
        std::ref(replica_label_group[isolated_part]),
        std::ref(part2data[isolated_part]),
        embedding,
        std::ref(permutations),
        std::ref(constraints),
        std::ref(full_constraints)
    );

    while (true) {
        uint32_t max_dist = 0, max_deg = 0;
        int choice = -1;
        for (uint32_t j : available_partition) {
            if (alloc_bitmap[j]) continue;
            if (min_dist[j] > max_dist) {
                max_dist = min_dist[j];
                max_deg = PDG_deg[j];
                choice = j;
            }
            else if (min_dist[j] == max_dist) {
                if (PDG_deg[j] > max_deg) {
                    max_dist = min_dist[j];
                    max_deg = PDG_deg[j];
                    choice = j;
                }
            }
        }
        if (choice == -1) break;
        alloc_bitmap[choice] = 1;
        partition_thread.emplace_back(
            ParallelSubgraphMatching, 
            std::ref(pool),
            k,
            query_graph, data_graph,
            partition_graphs[choice],
            choice,
            std::ref(replica_label_group[choice]),
            std::ref(part2data[choice]),
            embedding,
            std::ref(permutations),
            std::ref(constraints),
            std::ref(full_constraints)
        );

        // update
        for (uint32_t j : available_partition) {
            if (!alloc_bitmap[j]) {
                min_dist[j] = std::min(min_dist[j], distance_matrix[choice][j]);
            }
        }
    }
    for (auto& t : partition_thread) {
        if (t.joinable()) t.join();
    }
    std::vector<uint32_t> candidates, selected;
    if (_terminate_) goto EXIT_PDD_PLUS;

    // if not enough, retry with PDD
    puts("\nRetry with PDD...");
    partition_thread.clear();
    for (auto j : available_partition) candidates.push_back(j);
BACKTRACKING:
    selected.clear();
    SelectKFarthestPartition(k, candidates, selected, distance_matrix, alloc_bitmap, n_partition);
    if (selected.size() == 0) goto EXIT_PDD_PLUS;
    for (auto i = selected.begin(); i != selected.end(); ++i) {
        partition_thread.emplace_back(
            ParallelSubgraphMatching, 
            std::ref(pool),
            k,
            query_graph, data_graph,
            partition_graphs[*i],
            *i,
            std::ref(replica_label_group[*i]),
            std::ref(part2data[*i]),
            embedding,
            std::ref(permutations),
            std::ref(constraints),
            std::ref(full_constraints)
        );
    }
    for (auto& t : partition_thread) {
        if (t.joinable()) t.join();
    }
    partition_thread.clear();
    if (!_terminate_) goto BACKTRACKING;
EXIT_PDD_PLUS:

    auto end = std::chrono::high_resolution_clock::now();

    candidates.clear();
    selected.clear();
    min_dist.clear();
    delete alloc_bitmap;
#endif

#ifdef OUTPUT_RESULT
    std::ofstream output_file("embedding.bin", std::ios::out | std::ios::binary);
    uint32_t query_size = query_graph->getVerticesCount();
    output_file.write(reinterpret_cast<const char*>(&embedding_cnt), sizeof(uint32_t));
    output_file.write(reinterpret_cast<const char*>(&query_size), sizeof(uint32_t));
    output_file.write(reinterpret_cast<const char*>(embedding), embedding_cnt * query_size * sizeof(std::pair<uint32_t, uint32_t>));
    output_file.close();
    std::cout << "Embeddings saved to embedding.bin\n";
    delete[] embedding;
#endif
    std::cout << "Embedding Count: " << embedding_cnt << '\n';
    std::chrono::duration<double> duration = end - start;
    std::cout << "Total time: " << duration.count() << " s\n";

    delete query_graph;
    for (int i = 0; i < n_partition; ++i) {
        delete partition_graphs[i];
    }
    delete data_graph;
    delete[] inter_terminate_flag;
    for (int i = 0; i < n_partition; ++i) {
        delete[] distance_matrix[i];
    }
    delete[] distance_matrix;
    return 0;
}