#include <iostream>
#include <fstream>
#include <vector>
#include <cstdint>
#include <utility>
#include <queue>

using std::pair;
using std::make_pair;

int n, m, cnt;
bool vis[4000000];
struct node{
    node* nxt;
    int to;
    node() {
        nxt = NULL;
        to = 0;
    }
    ~node() {
        nxt = NULL;
        to = 0;
    }
};

void add_edge(node** head, node* edge, int x, int y) {
    ++cnt;
    edge[cnt].nxt = head[x];
    head[x] = &edge[cnt];
    head[x]->to = y;
    return;
}

void Dijkstra(node** head, int* dis, int s) {
    std::priority_queue<pair<int, int> > Q;
    for (int i = 0; i < n; ++i) {
        dis[i] = 0x3f3f3f3f;
        vis[i] = 0;
    }
    Q.push(make_pair(0, s));
    dis[s] = 0;
    while(!Q.empty()){
        int u = Q.top().second;
        Q.pop();
        if(vis[u]) continue;
        vis[u] = 1;
        for (node *i = head[u]; i != NULL; i = i->nxt) {
            int v = i->to;
            if (dis[v] > dis[u] + 1) {
                dis[v] = dis[u] + 1;
                Q.push(make_pair(-dis[v], v));
            }
        }
    }
    return;
}

int main(int argc, char** argv) {   // dataset graph, embedding file
    std::ifstream input_file(argv[2], std::ios::in | std::ios::binary);
    if (!input_file.is_open()) {
        std::cerr << "Error: Could not open embedding.bin for reading.\n";
        return 1;
    }

    uint32_t embedding_cnt;
    input_file.read(reinterpret_cast<char*>(&embedding_cnt), sizeof(uint32_t));
    if (!input_file) {
        std::cerr << "Error: Failed to read embedding_cnt.\n";
        input_file.close();
        return 1;
    }

    uint32_t query_size;
    input_file.read(reinterpret_cast<char*>(&query_size), sizeof(uint32_t));
    if (!input_file) {
        std::cerr << "Error: Failed to read query_size.\n";
        input_file.close();
        return 1;
    }

    size_t total_pairs = static_cast<size_t>(embedding_cnt) * static_cast<size_t>(query_size);
    std::vector<std::pair<uint32_t, uint32_t>> embedding(total_pairs);

    input_file.read(
        reinterpret_cast<char*>(embedding.data()),
        total_pairs * sizeof(std::pair<uint32_t, uint32_t>)
    );
    if (!input_file) {
        std::cerr << "Error: Failed to read embedding data.\n";
        input_file.close();
        return 1;
    }
    input_file.close();

    // std::cout << "embedding_cnt: " << embedding_cnt << "\n";
    // std::cout << "query_size: " << query_size << "\n";
    // std::cout << "Total pairs read: " << total_pairs << "\n";

    // std::cout << "\nElements of embedding:\n";
    // for (size_t i = 0; i < embedding_cnt; ++i) {
    //     printf("Embedding %d:\n", i);
    //     for (size_t j = i * query_size; j < (i + 1) * query_size; ++j) {
    //         printf("(%d, %d) ", embedding[j].first, embedding[j].second);
    //     }
    //     printf("\n");
    // }

    node* edge = new node[34000000];
    node** head = new node*[4000000];

    input_file.open(argv[1], std::ios::in);
    char type;
    int u, v;
    while (input_file >> type) {
        if (type == 't') input_file >> u >> v;
        if (type == 'v') {
            input_file >> u >> v >> u;
            ++n;
        }
        if (type == 'e') {
            input_file >> u >> v;
            add_edge(head, edge, u, v);
            add_edge(head, edge, v, u);
            ++m;
        }
    }

    int** dis = new int*[query_size];
    for (int i = 0; i < query_size; ++i) {
        dis[i] = new int[n];
    }

    double total_dis = 0.0;
    for (int i = 0; i < embedding_cnt; ++i) {
        for (int  j = i * query_size; j < (i + 1) * query_size; ++j) {
            Dijkstra(head, dis[j - i * query_size], embedding[j].second);
        }
        for (int j = i + 1; j < embedding_cnt; ++j) {
            int sum = 0x7fffffff;
            for (int  k = i * query_size; k < (i + 1) * query_size; ++k) {
                for (int  l = j * query_size; l < (j + 1) * query_size; ++l) {
                    sum = std::min(sum, dis[k - i * query_size][embedding[l].second]);
                }
            }
            total_dis += sum;
        }
    }

    std::cout << "DISTANCE: " << total_dis / (((embedding_cnt - 1) * embedding_cnt) / 2) << '\n';

    delete[] edge;
    delete[] head;
    for (int i = 0; i < query_size; ++i) {
        delete[] dis[i];
    }
    delete[] dis;
    return 0;
}