#include "BasicGraph.h"
#include <queue>
#include <iostream>
#include <fstream>
#include <set>
#include <unordered_set>
#include <immintrin.h>
#include <bitset>

using std::pair;
using std::make_pair;

GraphOP::GraphOP() {
    head = nullptr;
    edge = nullptr;
    degree = nullptr;
}

GraphOP::GraphOP(pair<ui, ui> max_index) {
    node_count = max_index.first + 1;
    edge_count = max_index.second;
    head = new Node*[node_count + 1];
    edge = new Node[(edge_count + 1) << 1];
    degree = new ui[node_count + 1];
    for (int i = 0; i < node_count; ++i) {
        head[i] = nullptr;
        degree[i] = 0;
    }
    cnt = 0;
}

GraphOP::~GraphOP() {
    delete[] degree;
    delete[] head;
    delete[] edge;
}

void GraphOP::add_edge(VertexID u, VertexID v) {
    ++cnt;
    edge[cnt].next_ = head[u];
    head[u] = &edge[cnt];
    head[u]->to = v;
    return;
}

void GraphOP::build(std::vector<pair<VertexID, VertexID> >& edges) {
    for (auto i = edges.begin(); i != edges.end(); ++i) {
        add_edge(i->first, i->second);
        add_edge(i->second, i->first);
        ++degree[i->first];
        ++degree[i->second];
    }
    return;
}

BFS_TREE* GraphOP::pre_process(VertexID start_vertex, Graph* query_graph) {
    BFS_TREE* root = new BFS_TREE;
    bool* visited = new bool[node_count + 1];
    for (int i = 0; i < node_count; ++i) {
        visited[i] = false;
    }
    SNode* sroot = new SNode;
    std::queue<BFS_TREE*> Q;
    Q.push(root);
    visited[start_vertex] = true;
    sroot->degree = degree[start_vertex];
    sroot->label = query_graph->getVertexLabel(start_vertex);
    sroot->neigh.resize(node_count + 1, 0);
    root->tree_node.push_back(sroot);
    root->index.push_back(start_vertex);
    BFS_TREE* now_node = root;
    while (!Q.empty()) {
        BFS_TREE* u = Q.front();
        Q.pop();
        // if (start_vertex == 3) std::cout << u->index[0] << '\n';
        for (Node* i = head[u->index[0]]; i != nullptr; i = i->next_) {
            VertexID v = i->to;
            if (visited[v]) {
                if (u->tree_node[0]->neigh[v] == 0) ++u->tree_node[0]->neigh_cnt;
                u->tree_node[0]->neigh[v] = 1;
            }
            else {
                // if (start_vertex == 3) std::cout << u->index[0] << ' ' << v << '\n';
                now_node = new BFS_TREE;
                sroot = new SNode;
                sroot->degree = degree[v];
                sroot->label = query_graph->getVertexLabel(v);
                sroot->neigh.resize(node_count + 1, 0);
                now_node->index.push_back(v);
                now_node->tree_node.push_back(sroot);
                u->child.push_back(now_node);
                Q.push(now_node);
                visited[v] = true;
            }
        }
    }
    delete[] visited;
    return root;
}

void GraphOP::Match_BFS(std::vector<VertexID>& start_vertex, BFS_TREE* Qroot, std::vector<ui>& labels, Graph* query_graph, int data_graph_index, GraphOP* QG, std::unordered_set<VertexID>* Graph_Matrix) {

    auto start = std::chrono::high_resolution_clock::now();
    auto end = std::chrono::high_resolution_clock::now();
    double time_in_ns;

    // bool* visited = new bool[node_count + 1];
    const size_t MAX_NODE = 30000000;
    // const size_t MAX_DEGREE = 2000;
    std::bitset<MAX_NODE> visited;
    int qsize = query_graph->getVerticesCount();
    // std::vector<BFS_TREE*>* Matched = new std::vector<BFS_TREE*>[node_count + 1];
    std::map<VertexID, std::vector<BFS_TREE*> > Matched;
    // std::unordered_set<BFS_TREE*>* Matched = new std::unordered_set<BFS_TREE*>[node_count + 1];
    // std::set<BFS_TREE*>* Matched = new std::set<BFS_TREE*>[node_count + 1];
    BFS_TREE** table = new BFS_TREE*[qsize + 1];
    for (int i = 0; i < qsize; ++i) {
        table[i] = nullptr;
    }
    // for (int i = 0; i < node_count; ++i) {
    //     visited[i] = false;
    // }
    // memset(visited, 0, sizeof(visited));
    // puts("Finish Initialize");
    BFS_TREE* root = new BFS_TREE;
    std::queue<BFS_TREE*> Q;
    SNode* sroot = nullptr;
    for (auto i = start_vertex.begin(); i != start_vertex.end(); ++i) {
        visited[*i] = 1;
        // sroot = new SNode;
        // sroot->degree = degree[*i];
        // sroot->label = labels[*i];
        // root->tree_node.push_back(sroot);
        root->index.push_back(*i);
        root->Qre = Qroot;
        root->vaild[*i] = true;
        Matched[*i].push_back(root);
        // Matched[*i].insert(root);
    }
    Q.push(root);
    BFS_TREE* now_node = nullptr;
    table[Qroot->index[0]] = root;
    std::map<VertexID, VertexID> MatchResult_own;
    std::vector<pair<VertexID, BFS_TREE*> > MatchSave;
    // std::vector<pair<pair<BFS_TREE*, int>, int> > stack_;
    bool* MatchL = new bool [qsize];
    // puts("Finish Start Vertices");
    while (!Q.empty()) {
        if (stopFlag[data_graph_index].load()) goto End;
        BFS_TREE* u = Q.front();
        Q.pop();
        // std::cout << Q.size() << ' ' << u << '\n';
        for (int t = 0; t < u->index.size(); ++t) {
            if (finished > _k) return;
            // std::cout << u->index[t] << '\n';
            // if (!u->vaild[u->index[t]]) continue;
            std::vector<bool> u_neigh;
            std::vector<VertexID> expand_able;
            u_neigh.resize(qsize, 0);
            for (Node* i = head[u->index[t]]; i != nullptr; i = i->next_) { // check if u is vaild
                VertexID v = i->to;
                // std::cout << u->index[t] << ' ' << v << '\n';
                if (visited[v]) {
                    for (int j = 0; j < Matched[v].size(); ++j) {
                        // std::cout << u->Qre->index[0] << "R\n";
                        if (Matched[v][j]->vaild[v]) u_neigh[Matched[v][j]->Qre->index[0]] = 1;
                    }
                    // for (auto j = Matched[v].begin(); j != Matched[v].end(); ++j) {
                    //     u_neigh[(*j)->Qre->index[0]] = 1; 
                    // }
                }
                else expand_able.push_back(v);
            }
            // puts("Finish Check");
            // std::cout << count << ' ' << u->Qre->tree_node[0]->neigh_cnt << '\n';
            ui count = 0;
            for (int i = 0; i < qsize; ++i) {
                if (u_neigh[i] == 1 && u->Qre->tree_node[0]->neigh[i] == 1) ++count;
            }
            if (count < u->Qre->tree_node[0]->neigh_cnt) {
                u->vaild[u->index[t]] = false;
                // Matched[u->index[t]].erase(u);
                continue;
            }
            // puts("Finish RE");
            
            for (int i = 0; i < qsize; ++i) MatchL[i] = false;
            MatchSave.clear();

            for (auto i = expand_able.begin(); i != expand_able.end(); ++i) { // expand
                // VertexID v = i->to;
                VertexID v = *i;
                // std::cout << u->Qre->index[0] << '\n';
                for (int j = 0; j < u->Qre->child.size(); ++j) {
                    // std::cout << u->Qre->child[j]->index[0] << "r\n";
                    if (labels[v] == u->Qre->child[j]->tree_node[0]->label && degree[v] >= u->Qre->child[j]->tree_node[0]->degree) {
                        // std::cout << v << "r\n";
                        MatchL[u->Qre->child[j]->index[0]] = true;
                        MatchSave.push_back(make_pair(v, u->Qre->child[j]));
                    }
                }
            }
            int total = 0;
            for (int i = 0; i < qsize; ++i) {
                total += (MatchL[i] == true);
            }
            if (total < u->Qre->child.size()) {
                u->vaild[u->index[t]] = false;
                continue;
            }
            else {
                for (auto i = MatchSave.begin(); i != MatchSave.end(); ++i) {
                    if (table[i->second->index[0]] == nullptr) {
                        BFS_TREE* new_node = new BFS_TREE;
                        table[i->second->index[0]] = new_node;
                        new_node->Qre = i->second;
                        Q.push(new_node);
                        u->child.push_back(new_node);
                    }
                    now_node = table[i->second->index[0]];
                    // sroot = new SNode;
                    // sroot->degree = degree[v];
                    // sroot->label = labels[v];
                    // now_node->tree_node.push_back(sroot);
                    now_node->index.push_back(i->first);
                    now_node->vaild[i->first] = true;
                    Matched[i->first].push_back(now_node);
                    // Matched[v].insert(now_node);
                    // u->kids.push_back(make_pair(u->index[t], i->first));
                    visited[i->first] = 1;
                }
            }
        }
        // puts("");
    }
    
    end = std::chrono::high_resolution_clock::now();
    time_in_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    // std::cout << "time_BFS_match: " << time_in_ns << '\n';
    start = std::chrono::high_resolution_clock::now();

    // puts("Finish BFS");
    // std::swap(table[2]->index[0], table[2]->index[1]);
    for (int i = 0; i < qsize; ++i) {
        if (table[i] == nullptr) return;
        // for (auto j = table[i]->index.begin(); j != table[i]->index.end(); ++j) {
        //     // if (table[i]->vaild[*j]) std::cout << *j << ' ';
        //     std::cout << *j << ' '; 
        // }
        // puts("");
    }
    
    if (stopFlag[data_graph_index].load() || finished > _k) goto End;

    Match_DFS(root, QG, table, Graph_Matrix, MatchResult_own);

    if (MatchResult_own.size() < qsize) goto End;
    if (!stopFlag[data_graph_index].exchange(true)) {
        std::lock_guard<std::mutex> lock(resultMutex[data_graph_index]);
        // MatchResult[data_graph_index].assign(MatchResult_own.begin(), MatchResult_own.end());
        for (auto i = MatchResult_own.begin(); i != MatchResult_own.end(); ++i) {
            MatchResult[data_graph_index].push_back(*i);
        }
    }

    end = std::chrono::high_resolution_clock::now();
    time_in_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
    // std::cout << "time_DFS_match: " << time_in_ns << '\n';
    timeMutex.lock();
    time_span.push_back(make_pair(start, end));
    timeMutex.unlock();
    // puts("PPS");
    End:
    std::map<BFS_TREE*, bool> del;
    // delete[] visited;
    for (auto i = Matched.begin(); i != Matched.end(); ++i) {
        for (auto j = i->second.begin(); j != i->second.end(); ++j) {
            BFS_TREE* now = *j;
            if (del.find(now) != del.end()) continue;
            if (now != nullptr) {
                delete now;
                del[now] = 1;
            }
        }
        i->second.clear();
    }
    delete[] table;
    delete[] MatchL;
    Matched.clear();
    MatchResult_own.clear();
    del.clear();
    return;
}

void GraphOP::Match_DFS(BFS_TREE* root, GraphOP* QG, BFS_TREE** table, std::unordered_set<VertexID>* Graph_Matrix, std::map<VertexID, VertexID>& Matchres) {

    bool* ensured = new bool[QG->node_count];
    int* finished = new int[QG->node_count];
    int* now_child = new int[QG->node_count];
    for (int i = 0; i < QG->node_count; ++i) {
        ensured[i] = false;
        finished[i] = 0;
        now_child[i] = 0;
    }
    std::vector<BFS_TREE*> node_stack;

    node_stack.push_back(root);
    while (!node_stack.empty()) {
        // if (Matchres.size() >= QG->node_count) break;
        BFS_TREE* u = node_stack.back();
        // std::cout << u->Qre->index[0] << '\n';
        // for (auto i = Matchres.begin(); i != Matchres.end(); ++i) {
        //     std::cout << i->first << ' ' << i->second << '\n';
        // }
        // puts("_____________________________________________________");
        if (::finished > ::_k) return;

        if (ensured[u->Qre->index[0]]) {
            if (now_child[u->Qre->index[0]] > 0 && ensured[now_child[u->Qre->index[0]] - 1] == false) {
                ensured[u->Qre->index[0]] = false;
                now_child[u->Qre->index[0]] = 0;
                continue;
            }
            if (now_child[u->Qre->index[0]] < u->child.size()) {
                node_stack.push_back(u->child[now_child[u->Qre->index[0]]]);
                ++now_child[u->Qre->index[0]];
            }
            else {
                node_stack.pop_back();
                now_child[u->Qre->index[0]] = 0;
                finished[u->Qre->index[0]] = 0;
            }
            continue;
        }
        bool success = false;
        for (int i = finished[u->Qre->index[0]]; i < u->index.size(); ++i) {
            if (!u->vaild[u->index[i]]) continue;
            bool tag = true;
            for (Node* j = QG->head[u->Qre->index[0]]; j != nullptr; j = j->next_) {
                VertexID v = j->to;
                if (ensured[v]) {
                    if (Graph_Matrix[u->index[i]].find(Matchres[v]) == Graph_Matrix[u->index[i]].end()) {
                        tag = false;
                        break;
                    }
                }
                else {
                    bool rtag = false;
                    for (auto k = table[v]->index.begin(); k != table[v]->index.end(); ++k) {
                        if (Graph_Matrix[u->index[i]].find(*k) != Graph_Matrix[u->index[i]].end()) {
                            rtag = true;
                            break;
                        }
                    }
                    if (!rtag) {
                        tag = false;
                        break;
                    }
                }
            }
            success |= tag;
            if (!tag) continue;
            Matchres.insert(make_pair(u->Qre->index[0], u->index[i]));
            ensured[u->Qre->index[0]] = true;
            finished[u->Qre->index[0]] = i + 1;
            break;
        }
        if (!success) {
            node_stack.pop_back();
            finished[u->Qre->index[0]] = 0;
            now_child[u->Qre->index[0]] = 0;
            ensured[u->Qre->index[0]] = false;
            Matchres.erase(u->Qre->index[0]);
        }
    }
    delete[] ensured;
    delete[] finished;
    delete[] now_child;
    return;
}
