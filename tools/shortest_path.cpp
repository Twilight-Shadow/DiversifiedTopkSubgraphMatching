#include <bits/stdc++.h>
using namespace std;
long long f[1005][1005];
int main() {
    // Partition Matrix stores the adjacency matrix of partitions
    std::ifstream file("./Partition_Matrix.txt", ios::in);
    for (int i = 0; i < 1000; ++i) {
        for(int j = 0; j < 1000; ++j) {
            file >> f[i][j];
            if (!f[i][j]) f[i][j] = 0x7fffffff;
        }
    }
    file.close();
    for (int i = 0; i < 1000; ++i) f[i][i] = 0;
    for (int k = 0; k < 1000; ++k) {
        for (int i = 0; i < 1000; ++i) {
            for (int j = 0; j < 1000; ++j) {
                f[i][j] = std::min(f[i][j], f[i][k] + f[k][j]);
            }
        }
    }
    std::ofstream ofile("./Distance_Matrix.txt", ios::out);
    for (int i = 0; i < 1000; ++i) {
        for (int j = 0; j < 1000; ++j) {
            ofile << f[i][j] << ' ';
        }
        ofile << '\n';
    }
    ofile.close();
    return 0;
}