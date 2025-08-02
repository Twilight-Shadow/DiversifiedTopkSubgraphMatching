# Diversified Top-k Subgraph Matching

## STEP-1: Graph Partitioning

Any edge partitioning method can be used for this step, we recommend to use [Distributed NE](https://github.com/masatoshihanai/DistributedNE) or [SHEEP](https://github.com/dmargo/sheep) .

## STEP-2: Partition Prediction

The source code of partition prediction algorithm can be found in [SPMiner](https://github.com/snap-stanford/neural-subgraph-learning-GNN), we only use the subgraph matching part of the experimental branch, because only this branch takes node features into account. Our implementation can be found in `prediction/` .

### Environment Setup

You can refer to `prediction/requirements.txt` to install dependencies.

### Dataset

For training, validation and testing, we all use the synthetic dataset that is used in SPMiner. The model only needs to be trained once, and it can be used in any other dataset in our experiment.

### Training

Train the encoder by command: `python3 -m subgraph_matching.train --node_anchored`. The model will be saved to `ckpt/model.pt`. You can analyze the trained encoder via `python3 -m subgraph_matching.test --node_anchored`.

### Testing

The module `python3 -m subgraph_matching.alignment.py [--query_path=...] [--target_path=...]` provides a utility to obtain all pairs of corresponding matching scores, given a pickle file of the query and target graphs in networkx format. Run the module without these arguments for an example using random graphs. 

## STEP-3: Subgraph Matching

### Environment Setup

`GCC 11.4.0` and `CMake` are needed for compilation. To efficiently manage threads, we use [CTPL](https://github.com/vit-vit/ctpl) library, which has been already involved in the directory.

### Build

```bash
mkdir build && cd build
cmake ..
make -j4
```

### Execuate

For a quick start, you can run the following command in bash:

```bash
./build/matching/SubgraphMatching.out -q test/query.graph -d test/partitions/ -n 4 -k 3 -g test/test.graph -p 3 -m test/distance_matrix.txt -t test/ground_truth.txt -r test/node_replica/
```

which finds top-3 diversified matches of `query.graph` in `test.graph` using algorithm PDD+.

### Data Preparation

**Query graph**, **data graph** and **partition graph** are organized in CSR format:

``` 
t <vertex count> <edge count>
v <vertex ID> <vertex label> <vertex degree>
...
e <src vertex> <dst vertex>
...
```

All partition graphs should be in one folder and be named from `0.graph` to `[partition_number - 1].graph`. All graphs used in our experiment are undirected graph with node label.

Since we use edge partitioning method for graph partitioning, the vertices that are replicated between partitions should be recorded for inter-partition subgraph matching. **If a vertex is in more than one partition simultaneously, then it is a replicated node.** For each partition, the index of replicated nodes should be stored in `[partition_index].txt`, and these txt files should be in one folder.

A **distance matrix** stores the pairwise distance between any two partitions in `PAG(Partition Adjacency Graph)`. Two partitions are adjacent if the intersection set between the vertex lists of the two partitions is not empty. After building PAG, you can run shortest path algorithm such as `Floyed` or `Dijkstra` on PAG to calculate the pairwise distance.

The **ground truth** file from partition prediction stage should be provided when running PDD+. It stores 0/1 value, where the $i_{th}$ one represents whether the given query graph is in $i_{th}$ partition graph or not (0-not exist, 1-exist).

### Options & Configurations

```bash
-q: the path of query graph file
-g: the path of data graph file
-d: the path of folder which stores all partition graphs
-r: the path of folder which stores replicated nodes in each partition
-m: the path of distance matrix file
-t: the path of ground truth file
-k: the top-k value
-n: number of thread
-p: number of partition
```

You can switch between PDD and PDD+ and choose to write embeddings to `embedding.bin` in `configuration/Config.h` . File `show_embedding.cpp` helps to display all embeddings.

### NOTES:

Most of the offline processing tools can be found in `tools/`.

Vertex partitioning and hybrid partitioning are not supported due to the design of inter-partition matching algorithm.

The algorithm used in intra-partition matching phase can be replaced by any subgraph matching algorithm, in our experiemt we use a hybird way, which is GQL filter + GQL order + LFTJ  enumeration. 

