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
#include "ctpl/ctpl_stl.h"

typedef unsigned int ui;
typedef uint32_t VertexID;
typedef ui LabelID;

using std::pair;
using std::make_pair;

extern std::atomic<bool>* stopFlag;
extern std::mutex* resultMutex;
extern std::vector<pair<VertexID, VertexID> >* MatchResult;
extern std::mutex timeMutex;
extern std::mutex finishMutex;
extern std::vector<pair<std::chrono::time_point<std::chrono::high_resolution_clock>, std::chrono::time_point<std::chrono::high_resolution_clock> > > time_span;
extern int finished;
extern int _k;
extern bool* status;
extern int running_count;
extern std::mutex rMutex;

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

class SNode {
public:
    ui label, degree;
    std::vector<VertexID> neigh;
    ui neigh_cnt;
public:
    SNode() {
        neigh_cnt = 0;
        neigh.clear();
    }
    ~SNode() {
        neigh.clear();
    }
};

class BFS_TREE {
public:
    std::vector<VertexID> index;
    std::vector<SNode*> tree_node;
    std::vector<BFS_TREE*> child;
    std::map<VertexID, bool> vaild;
    // std::vector<pair<VertexID, VertexID> > kids;
    // std::set<VertexID> tempM;
    BFS_TREE* Qre;
    BFS_TREE() {
        Qre = nullptr;
    }
    ~BFS_TREE() {
        index.clear();
        tree_node.clear();
        child.clear();
        vaild.clear();
        for (auto i = tree_node.begin(); i != tree_node.end(); ++i) {
            delete *i;
        }
    }
};

struct clique {
    int deg, len, idx;
    bool operator < (const clique& x) const {
        if (x.len == len) return deg < x.deg;
        else return len < x.len;
    }
};

#endif // TYPES_H