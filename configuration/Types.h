#ifndef TYPES_H
#define TYPES_H

#include <cstdint>
#include <cstdlib>
#include <vector>
#include <map>
#include <atomic>
#include <mutex>
#include <condition_variable>
#include <chrono>
#include <future>
#include <thread>
#include <pthread.h>
#include <semaphore.h>

typedef unsigned int ui;
typedef uint32_t VertexID;
typedef ui LabelID;

using std::pair;
using std::make_pair;

extern uint32_t embedding_cnt;
extern bool _terminate_;
extern std::mutex emb_lock;
extern bool* inter_terminate_flag;

enum MatchingIndexType {
    VertexCentric = 0,
    EdgeCentric = 1
};

class TreeNode {
public:
    VertexID id_;
    VertexID parent_;
    ui level_;
    ui under_level_count_;
    ui children_count_;
    ui bn_count_;
    ui fn_count_;
    VertexID* under_level_;
    VertexID* children_;
    VertexID* bn_;
    VertexID* fn_;
    size_t estimated_embeddings_num_;
public:
    TreeNode() {
        id_ = 0;
        under_level_ = NULL;
        bn_ = NULL;
        fn_ = NULL;
        children_ = NULL;
        parent_ = 0;
        level_ = 0;
        under_level_count_ = 0;
        children_count_ = 0;
        bn_count_ = 0;
        fn_count_ = 0;
        estimated_embeddings_num_ = 0;
    }

    ~TreeNode() {
        delete[] under_level_;
        delete[] bn_;
        delete[] fn_;
        delete[] children_;
    }

    void initialize(const ui size) {
        under_level_ = new VertexID[size];
        bn_ = new VertexID[size];
        fn_ = new VertexID[size];
        children_ = new VertexID[size];
    }
};

class Edges {
public:
    ui* offset_;
    ui* edge_;
    ui vertex_count_;
    ui edge_count_;
    ui max_degree_;
public:
    Edges() {
        offset_ = NULL;
        edge_ = NULL;
        vertex_count_ = 0;
        edge_count_ = 0;
        max_degree_ = 0;
    }

    ~Edges() {
        delete[] offset_;
        delete[] edge_;
    }
};

class Node {
public:
    VertexID to;
    Node* next_;
public:
    Node() {
        to = 0;
        next_ = nullptr;
    }
    
    ~Node() {
        to = 0;
        next_ = nullptr;
    }
};

#endif // TYPES_H