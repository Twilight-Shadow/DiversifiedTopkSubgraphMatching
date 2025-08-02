#include <iostream>
#include <fstream>
#include <vector>
#include <cstdint>
#include <utility>

int main() {
    std::ifstream input_file("embedding.bin", std::ios::in | std::ios::binary);
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

    std::cout << "embedding_cnt: " << embedding_cnt << "\n";
    std::cout << "query_size: " << query_size << "\n";
    std::cout << "Total pairs read: " << total_pairs << "\n";

    std::cout << "\nElements of embedding:\n";
    for (size_t i = 0; i < embedding_cnt; ++i) {
        printf("Embedding %d:\n", i);
        for (size_t j = i * query_size; j < (i + 1) * query_size; ++j) {
            printf("(%d, %d) ", embedding[j].first, embedding[j].second);
        }
        printf("\n");
    }
    return 0;
}