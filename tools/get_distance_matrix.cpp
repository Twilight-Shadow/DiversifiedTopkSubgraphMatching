#include <bits/stdc++.h>
#include "ctpl/ctpl_stl.h"
#include <immintrin.h>
using namespace std;

size_t simd_intersect(const uint32_t* A, size_t lenA, const uint32_t* B, size_t lenB, uint32_t* out) {
    size_t i = 0, j = 0, count = 0;

    while (i < lenA && j < lenB) {
        // Process 4 elements from A and 4 from B
        if (i + 4 <= lenA && j + 4 <= lenB) {
            __m128i vA = _mm_loadu_si128((__m128i*)&A[i]);
            __m128i vB = _mm_loadu_si128((__m128i*)&B[j]);

            for (int k = 0; k < 4; ++k) {
                __m128i bcast = _mm_set1_epi32(A[i + k]);
                __m128i cmp = _mm_cmpeq_epi32(bcast, vB);
                int mask = _mm_movemask_ps(_mm_castsi128_ps(cmp));

                if (mask != 0) {
                    out[count++] = A[i + k];
                }
            }

            // Move forward
            if (A[i + 3] < B[j + 3]) {
                i += 4;
            } else if (A[i + 3] > B[j + 3]) {
                j += 4;
            } else {
                i += 4;
                j += 4;
            }
        } else {
            // Scalar fallback
            if (A[i] < B[j]) {
                ++i;
            } else if (A[i] > B[j]) {
                ++j;
            } else {
                out[count++] = A[i];
                ++i;
                ++j;
            }
        }
    }

    return count;
}

std::vector<int> dijkstra_heap(int id, int n, uint32_t** graph, int source) {
    std::vector<int> dist(n, 0x7fffffff);
    std::vector<bool> visited(n, false);

    // Min-heap: pair<distance, node>
    using pii = std::pair<int, int>;
    std::priority_queue<pii, std::vector<pii>, std::greater<pii>> pq;

    dist[source] = 0;
    pq.emplace(0, source);

    while (!pq.empty()) {
        auto [d, u] = pq.top(); pq.pop();
        if (visited[u]) continue;
        visited[u] = true;

        for (int v = 0; v < n; ++v) {
            if (graph[u][v] != 0x7fffffff && !visited[v]) {
                if (dist[v] > dist[u] + graph[u][v]) {
                    dist[v] = dist[u] + graph[u][v];
                    pq.emplace(dist[v], v);
                }
            }
        }
    }

    return dist;
}

std::string pad(int num) {
    std::ostringstream oss;
    oss << std::setw(5) << std::setfill('0') << num;
    return oss.str();
}

int main(int argc, char** argv) {   // n_part, directory, target
    ctpl::thread_pool pool(32);
    std::ifstream file;
    int n = std::stoi(argv[1]);
    string dir(argv[2]);

    std::vector<std::vector<uint32_t> > p(n);

    for (int i = 0; i < n; ++i) {
        file.open(dir + "/" + pad(i), std::ios::in);
        uint32_t u, v;
        std::set<uint32_t> s;
        s.clear();
        while (file >> u >> v) {
            s.insert(u);
            s.insert(v);
        }
        for (uint32_t j : s) {
            p[i].push_back(j);
        }
        file.close();
    }

    uint32_t** pm = new uint32_t*[n];
    for (int i = 0; i < n; ++i) {
        pm[i] = new uint32_t[n];
        for (int j = 0; j < n; ++j) {
            pm[i][j] = 0x7fffffff;
        }
        pm[i][i] = 0;
    }

    std::vector<uint32_t> result(20000);
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (simd_intersect(p[i].data(), p[i].size(), p[j].data(), p[j].size(), result.data()) >= 1) {
                pm[i][j] = pm[j][i] = 1;
            }
        }
    }
    // for (int i = 0; i < n; ++i) {
    //     std::cout << p[i].size() << '\n';
    //     // for (int j = 0; j < n; ++j) {
    //     //     std::cout << pm[i][j] << ' ';
    //     // }
    //     // std::cout << '\n';
    // }

    std::vector<std::future<std::vector<int>>> f(n);

    for (int i = 0; i < n; ++i) {
        f[i] = pool.push(dijkstra_heap, n, pm, i);
    }
    pool.stop(true);

    string target(argv[3]);
    std::ofstream out(target, std::ios::out);
    for (int i = 0; i < n; ++i) {
        std::vector<int> d = f[i].get();
        for (int j : d) out << j << ' ';
        out << '\n';
    }
    out.close();
    return 0;
}