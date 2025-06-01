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
#include "ctpl/ctpl_stl.h"

#include "CommandList.h"
#include "graph/MatchGraph.h"
#include "graph/BasicGraph.h"
#include "GraphQL.h"
#include "BuildTable.h"

using std::pair;
using std::make_pair;

std::atomic<bool>* stopFlag;
std::mutex* resultMutex;
std::vector<pair<VertexID, VertexID> >* MatchResult;
std::mutex timeMutex;
std::mutex finishMutex;
std::vector<pair<std::chrono::time_point<std::chrono::high_resolution_clock>, std::chrono::time_point<std::chrono::high_resolution_clock> > > time_span;
int finished;
int _k;
bool* status;
int running_count;
std::mutex rMutex;

void ParallelSubgraphMatching(int id, std::string input_query_graph_file, std::string input_data_graph_path, int data_graph_index, 
                                GraphOP* G, pair<ui, ui> max_index, std::vector<ui>& label, int t_index, std::unordered_set<VertexID>* Graph_Matrix) {
    // std::cout << label.size() << '\n';
    rMutex.lock();
    ++running_count;
    // std::cout << "adding..." << data_graph_index << '\n';
    rMutex.unlock();
    Graph* query_graph = new Graph(true);
    query_graph->loadGraphFromFile(input_query_graph_file);

    query_graph->buildCoreTable();
    std::cout << "Finish Loading Query Graph" << data_graph_index << '\n';

    Graph* data_graph = new Graph(true);
    data_graph->loadGraphFromFile(input_data_graph_path + "/data_graph/" + std::to_string(data_graph_index) + ".graph");
    std::cout << "Finish Loading Data Graph" << data_graph_index << '\n';

    int qsize = query_graph->getVerticesCount();

    // auto start = std::chrono::high_resolution_clock::now();

    ui** candidates = NULL;
    ui* candidates_count = NULL;
    // std::vector<std::unordered_map<VertexID, std::vector<VertexID >>> TE_Candidates;
    // std::vector<std::vector<std::unordered_map<VertexID, std::vector<VertexID>>>> NTE_Candidates;
    FilterVertices::GQLFilter(data_graph, query_graph, candidates, candidates_count);
    FilterVertices::sortCandidates(candidates, candidates_count, (ui)qsize);
    std::cout << "Finish Finding Candidates" << data_graph_index << '\n';

    Edges ***edge_matrix = NULL;
    edge_matrix = new Edges **[qsize];
    for (ui i = 0; i < qsize; ++i) {
        edge_matrix[i] = new Edges *[qsize];
    }
    BuildTable::buildTables(data_graph, query_graph, candidates, candidates_count, edge_matrix);
    std::cout << "Finish Building Edge Matrix" << data_graph_index << '\n';

    ui* matching_order = NULL;
    ui* pivots = NULL;
    // size_t order_num = 100;
    GenerateQueryPlan::generateGQLQueryPlan(data_graph, query_graph, candidates_count, matching_order, pivots);
    GenerateQueryPlan::checkQueryPlanCorrectness(query_graph, matching_order, pivots);
    // std::cout << "Finish Generating Query Plan\n";

    size_t output_limit = 1;
    size_t embedding_count = 0;
    size_t call_count = 0;
    // size_t time_limit = 20050509;
    std::ofstream file("./MatchRes/MatchResult" + std::to_string(data_graph_index) + ".txt", std::ios::out);
    VertexID u, v;
    // std::map<VertexID, VertexID> map_info;
    // std::map<VertexID, VertexID> remap_info;
    VertexID* map_info = new VertexID[max_index.first + 2];
    VertexID* remap_info = new VertexID[max_index.first + 2];
    // for (int i = 0; i <= max_index.first; ++i) {
    //     map_info[i] = 0x7fffffff;
    //     remap_info[i] = 0x7fffffff;
    // }
    std::fstream ffile;
    ffile.open(input_data_graph_path + "/map_info/" + std::to_string(data_graph_index) + ".txt", std::ios::in);
    while (ffile >> u >> v) { // origin -- map
        map_info[u] = v;
        remap_info[v] = u;
    }
    ffile.close();
    ffile.clear();

    auto start = std::chrono::high_resolution_clock::now();

    embedding_count = EvaluateQuery::LFTJ(data_graph, query_graph, edge_matrix, candidates, candidates_count,
                                              matching_order, output_limit, call_count, data_graph_index, remap_info);

    auto end = std::chrono::high_resolution_clock::now();
    // std::cout << "time_basic_match: " << time_in_ns << '\n';
    timeMutex.lock();
    time_span.push_back(make_pair(start, end));
    timeMutex.unlock();

    if (embedding_count >= 1) {
        file << "Query matches " + std::to_string(data_graph_index) + ".graph" + " without overflow." << '\n';
        finishMutex.lock();
        ++finished;
        std::cout << "fin" << data_graph_index << ' ' << finished << '\n';
        finishMutex.unlock();
        status[data_graph_index] = 1;
    }
    else {
        start = std::chrono::high_resolution_clock::now();
        std::vector<VertexID> partition;
        ffile.open(input_data_graph_path + "/partition_point/" + std::to_string(data_graph_index) + ".txt", std::ios::in);
        VertexID u, v, w;
        while (ffile >> u) partition.push_back(u);
        ffile.close();
        ffile.clear();

        GraphOP* Query = new GraphOP(make_pair(qsize, query_graph->getEdgesCount()));
        std::vector<pair<VertexID, VertexID> > query_edges;
        ffile.open(input_query_graph_file);
        char t;
        while (ffile >> t) {
            if (t == 'v') ffile >> u >> v >> w;
            else ffile >> u >> v;
            if (t == 'e') query_edges.push_back(make_pair(u, v));
        }
        ffile.close();
        ffile.clear();
        Query->build(query_edges);

        end = std::chrono::high_resolution_clock::now();
        // time_in_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
        // std::cout << "time_load: " << time_in_ns << '\n';

        start = std::chrono::high_resolution_clock::now();

        BFS_TREE** roots = new BFS_TREE*[qsize];
        for (int i = 0; i < qsize; ++i) {
            roots[i] = Query->pre_process(i, query_graph);
        }

        end = std::chrono::high_resolution_clock::now();
        // time_in_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
        // std::cout << "time_pre: " << time_in_ns << '\n';

        // std::queue<BFS_TREE*> Q;
        // bool vis[8] = {0, 0, 0, 0, 0, 0, 0, 0};
        // Q.push(roots[2]);
        // vis[2] = 1;
        // while (!Q.empty()) {
        //     BFS_TREE* u = Q.front();
        //     Q.pop();
        //     std::cout << u->index[0] << '\n';
        //     for (int i = 0; i < u->tree_node[0]->neigh.size(); ++i) std::cout << u->tree_node[0]->neigh[i] << ' ';
        //     puts("");
        //     for (int i = 0; i < u->child.size(); ++i) {
        //         std::cout << u->child[i]->index[0] << "u ";
        //         if (vis[u->child[i]->index[0]]) continue;
        //         Q.push(u->child[i]);
        //         vis[u->child[i]->index[0]] = 1;
        //     }
        //     puts("");
        // }

        start = std::chrono::high_resolution_clock::now();

        std::vector<VertexID>* start_vertex = new std::vector<VertexID>[qsize];

        for (auto i = partition.begin(); i != partition.end(); ++i) {
            for (int j = 0; j < qsize; ++j) {
                if (label[*i] == query_graph->getVertexLabel(j) && G->degree[*i] >= Query->degree[j]) {
                    start_vertex[j].push_back(*i);
                }
            }
        }

        // std::vector<pair<VertexID, VertexID> > MatchResult;
        // for (int i = 0; i < query_graph->getVerticesCount(); ++i) {
        //     MatchResult = G->Match_BFS(start_vertex[i], roots[i], label, query_graph);
        //     if (MatchResult.size() >= query_graph->getVerticesCount()) break;
        // }

        std::thread* BFS_Thread = new std::thread[qsize];
        for (int i = 0; i < qsize; ++i) {
            BFS_Thread[i] = std::thread(&GraphOP::Match_BFS, G, std::ref(start_vertex[i]), roots[i], std::ref(label), query_graph, t_index, Query, Graph_Matrix);
            // std::cout << start_vertex[i].size() << '\n';
        }

        for (int i = 0; i < qsize; ++i) {
            if (BFS_Thread[i].joinable()) BFS_Thread[i].join();
        }

        end = std::chrono::high_resolution_clock::now();
        // time_in_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
        // std::cout << "time_match: " << time_in_ns << '\n';
        // MatchResult = G->Match_BFS(start_vertex[2], roots[2], label, query_graph);

        if (MatchResult[t_index].size() >= qsize) {
            file << "Query matches " + std::to_string(data_graph_index) + ".graph" + " with overflow." << '\n';
            std::ofstream ofile("./MatchRes/Matching" + std::to_string(data_graph_index) + ".txt", std::ios::out);
            for (auto i = MatchResult[t_index].begin(); i != MatchResult[t_index].end(); ++i) {
                ofile << "Matched vertices: " << i->first << ' ' << i->second << '\n';
            }
            ofile.close();
            finishMutex.lock();
            ++finished;
            finishMutex.unlock();
            status[data_graph_index] = 1;
        }
        else {
            file << "Query can not match " + std::to_string(data_graph_index) + ".graph" << '\n';
        }

        partition.clear();
        query_edges.clear();
        delete Query;
        delete[] roots;
        for (int i = 0; i < qsize; ++i) {
            start_vertex[i].clear();
        }
        delete[] start_vertex;
        delete[] BFS_Thread;
    }
    file.close();
    delete[] candidates_count;
    delete[] matching_order;
    delete[] pivots;
    for (ui i = 0; i < qsize; ++i) {
        delete[] candidates[i];
    }
    delete[] candidates;
    if (edge_matrix != NULL) {
        for (ui i = 0; i < qsize; ++i) {
            for (ui j = 0; j < qsize; ++j) {
                delete edge_matrix[i][j];
            }
            delete[] edge_matrix[i];
        }
        delete[] edge_matrix;
    }
    delete query_graph;
    delete data_graph;
    rMutex.lock();
    --running_count;
    // std::cout << "minusing..." << data_graph_index << '\n';
    rMutex.unlock();
    return;
}

void LoadLabel(std::string input_data_path, std::vector<ui>& label) {
    std::ifstream infile(input_data_path + "/all_labels.txt", std::ios::in);
    if (!infile.is_open()) {
        std::cerr << "Can not open graph file\n";
        exit(-1);
    }
    label.clear();
    ui u, v;
    while (infile >> u >> v) {
        label.push_back(v);
    }
    infile.close();
    return;
}

void LoadGraph(std::string input_data_path, std::vector<pair<VertexID, VertexID> >& edges, pair<ui, ui>& max_index) {
    std::ifstream infile(input_data_path + "/WG.txt");
    if (!infile.is_open()) {
        std::cerr << "Can not open graph file\n";
        exit(-1);
    }
    VertexID u, v;
    while (infile >> u >> v) {
        edges.push_back(make_pair(u, v));
        ++max_index.second;
        max_index.first = std::max(max_index.first, std::max(u, v));
    }
    infile.close();
    return;
}

int GetQueryIndex(std::string filename) {
    int i = filename.length() - 7;
    std::string num;
    for ( ; i; --i) {
        if (filename[i] < '0' || filename[i] > '9') break;
        num.push_back(filename[i]);
    }
    std::reverse(num.begin(), num.end());
    return std::stoi(num.c_str());
}

int main(int argc, char** argv){
    MatchingCommand command(argc, argv);
    std::string input_query_graph_file = command.getQueryGraphFilePath();
    std::string input_data_graph_path = command.getDataGraphFilePath();
    // std::string whole_graph_path = command.getWholeGraphPath();
    int k = std::stoi(command.getk());
    int np = std::stoi(command.getPartitionNum());
    _k = k;
    // int nq = std::stoi(command.getQueryNum());
    
    stopFlag = new std::atomic<bool>[2 * np];
    resultMutex = new std::mutex[2 * np];
    MatchResult = new std::vector<pair<VertexID, VertexID> >[2 * np];
    status = new bool[np];
    for (int i = 0; i < np; ++i) {
        status[i] = 0;
    }

    std::pair<ui, ui> max_index = make_pair(0, 0);
    std::vector<pair<VertexID, VertexID> > edges;
    LoadGraph(input_data_graph_path, edges, max_index);
    // std::cout << max_index.first << ' ' << max_index.second << '\n';
    // std::bitset<1200000>* Graph_Matrix = new std::bitset<1200000>[max_index.first + 1];
    std::unordered_set<VertexID>* Graph_Matrix = new std::unordered_set<VertexID>[max_index.first + 1];
    for (auto i = edges.begin(); i != edges.end(); ++i) {
        Graph_Matrix[i->first].insert(i->second);
        Graph_Matrix[i->second].insert(i->first);
    }

    GraphOP* G = new GraphOP(max_index);
    G->build(edges);

    std::vector<ui> label;
    LoadLabel(input_data_graph_path, label);

    for (int i = 0; i < np; ++i) stopFlag[i] = false;
    ctpl::thread_pool pool(1);
    std::vector<pair<int, int> > Partition_Graph[np];
    std::ifstream file(input_data_graph_path + "/Distance_Matrix.txt", std::ios::in);
    int** distance = new int*[np];
    for (int i = 0; i < np; ++i) {
        distance[i] = new int[np];
        for (int j = 0; j < np; ++j) {
            file >> distance[i][j];
        }
    }
    file.close();
    file.clear();

    int query_index = GetQueryIndex(input_query_graph_file);
    int* truth = new int[np];
    file.open(input_data_graph_path + "/ground_truth.txt", std::ios::in);
    for (int i = 0; i <= query_index; ++i) {
        for (int j = 0; j < np; ++j) {
            file >> truth[j];
        }
    }
    file.close();
    file.clear();

    // std::cout << query_index << '\n';
    // for (int i = 0; i < np; ++i) std::cout << truth[i] << ' ';
    // puts("");

    file.open(input_data_graph_path + "/Partition_Matrix.txt", std::ios::in);
    int* deg = new int[np];
    std::vector<pair<int, int> > Nodes;
    for (int i = 0; i < np; ++i) deg[i] = 0;
    for (int i = 0; i < np; ++i) {
        for (int j = 0, e; j < np; ++j) {
            file >> e;
            if (!e && truth[i] && truth[j]) {
                ++deg[i];
                Partition_Graph[i].push_back(make_pair(j, distance[i][j]));
            }
        }
    }

    int* used = new int[np];
    int* vis = new int[np];
    for (int i = 0; i < np; ++i) {
        if (truth[i]) Nodes.push_back(make_pair(deg[i], i));
        used[i] = vis[i] = 0;
    }

    std::sort(Nodes.begin(), Nodes.end(), std::greater<pair<int, int> >());
    std::priority_queue<clique> Q;
    if (Nodes.size()) {
        Q.push(clique{Nodes.begin()->first, 0, Nodes.begin()->second});
        used[Nodes.begin()->second] = 1;
    }
    int pushed = 0;

    while (!Q.empty() && finished < k) {
        int u = Q.top().idx;
        Q.pop();
        // std::cout << u << "\n";
        pool.push(ParallelSubgraphMatching, input_query_graph_file, input_data_graph_path, u, G, max_index, std::ref(label), pushed, Graph_Matrix);
        ++pushed;
        for (auto i = Partition_Graph[u].begin(); i != Partition_Graph[u].end(); ++i) {
            if (i->second > 1 && !used[i->first]) {
                Q.push(clique{deg[i->first], i->second, i->first});
                used[i->first] = 1;
            }
            // else used[i->first] = 2;
        }
    }
    // pool.stop(true);
    // sleep(5);
    // std::cout << pushed << ' ' << finished << ' ' << running_count << '\n';
    while (finished < k && running_count) {
        sleep(1);
        // std::cout << finished << ' ' << running_count << "u\n";
    }
    if (finished >= k) {
        pool.clear_queue();
        finishMutex.lock();
        ++finished;
        finishMutex.unlock();
    }

    if (finished < k) {
        if (Nodes.size()) Q.push(clique{Nodes.begin()->first, 0, Nodes.begin()->second});
        while (!Q.empty() && finished < k) {
            int u = Q.top().idx;
            Q.pop();
            if (used[u] != 1) {
                pool.push(ParallelSubgraphMatching, input_query_graph_file, input_data_graph_path, u, G, max_index, std::ref(label), pushed, Graph_Matrix);
                // std::cout << u << '\n';
                ++pushed;
            }
            for (auto i = Partition_Graph[u].begin(); i != Partition_Graph[u].end(); ++i) {
                if (vis[i->first] != 1) {
                    Q.push(clique{deg[i->first], i->second, i->first});
                    vis[i->first] = 1;
                }
            }
        }
    }
    // std::cout << pushed << ' ' << finished << ' ' << running_count << '\n';
    while (finished < k && running_count) {
        sleep(1);
        // std::cout << finished << ' ' << running_count << "v\n";
    }
    if (finished >= k) {
        pool.clear_queue();
        finishMutex.lock();
        ++finished;
        finishMutex.unlock();
    }

    if (finished < k) {
        for (int i = 0; i < np && finished < k; ++i) {
            if (used[i] != 1 && vis[i] != 1) {
                pool.push(ParallelSubgraphMatching, input_query_graph_file, input_data_graph_path, i, G, max_index, std::ref(label), pushed, Graph_Matrix);
                ++pushed;
            }
        }
    }

    while (finished < k && running_count) {
        sleep(1);
        // std::cout << finished << ' ' << running_count << "w\n";
    }
    std::cout << pushed << ' ' << finished << ' ' << running_count << '\n';
    if (finished >= k) {
        pool.clear_queue();
        finishMutex.lock();
        ++finished;
        finishMutex.unlock();
    }
    pool.stop(true);

    std::ofstream ofile("./timeres.txt", std::ios::out);
    std::vector<pair<std::chrono::time_point<std::chrono::high_resolution_clock>, std::chrono::time_point<std::chrono::high_resolution_clock> > > real_time_span;
    std::sort(time_span.begin(), time_span.end());
    for (auto i = time_span.begin(); i != time_span.end(); ++i) {
        if (real_time_span.empty()) {
            real_time_span.push_back(*i);
        }
        else {
            auto& lst = real_time_span.back();
            if (i->first < lst.second) lst.second = i->second;
            else real_time_span.push_back(*i);
        }
    }
    double time_match = 0;
    for (auto i = real_time_span.begin(); i != real_time_span.end(); ++i) {
        time_match += std::chrono::duration_cast<std::chrono::nanoseconds>(i->second - i->first).count();
    }
    ofile << time_match << '\n';
    ofile.close();
    ofile.clear();
    
    std::vector<int> MatchList;
    for (int i = 0; i < np; ++i) {
        if (status[i] == 1) MatchList.push_back(i);
    }
    ofile.open("./MatchList.txt", std::ios::out);
    ofile << MatchList.size() << '\n';
    for (auto i = MatchList.begin(); i != MatchList.end(); ++i) {
        ofile << *i << ' ';
    }
    ofile << '\n';
    ofile.close();

    delete G;
    for (int i = 0; i < max_index.first; ++i) {
        Graph_Matrix[i].clear();
    }
    delete[] Graph_Matrix;
    delete[] stopFlag;
    delete[] resultMutex;
    delete[] MatchResult;
    delete[] status;
    for (int i = 0; i < np; ++i) delete[] distance[i];
    delete[] distance;
    delete[] deg;
    delete[] truth;
    return 0;
}