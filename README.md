# DiversifiedTopkSubgraphMatching
**NOTICE:** The APIs of the current version is **BADLY DEFINED** due to extensive need of debugging and testing. A well refined version will be submitted later on main branch. The experiment data will be uploaded to Google Drive later. 

## Graph Partitioning

We just use [Distributed NE](https://github.com/masatoshihanai/DistributedNE) for this step.

## Ordering-Embedding

### Environment Setup

Run `pip -r requirement.txt` to install dependencies automatically.

### Training

Train the encoder: `python3 -m subgraph_matching.train` .

The code uses balanced synthetic dataset by default, which generate same number of positive pairs and negative pairs. You can also use custom datasets, which should be organized in following format:

```
<src node> <dst node> <src node label> <dst node label>
```

### Testing

Test the encoder: `python3 -m subgraph_matching.test` . The APIs for query and dataset have not been defined, so the input informations are hard-coded.

For input, user should provide one query, $n$ partitions and $n$ anchor files (which are used to identify sub-structures, defined by user). For one execution, the program will output a binary vector, $i_{th}$ number represents whether this query is likely appears in $i_{th}$ partition or not. Query graph and partition graph should be organized in the format above.

## PDD & PDD+

### Environment Setup

The program should be executed on a 64-bit machine equipped with `GCC 11.4.0` or above, along with `CMake`. The directiry of a dataset should be organized as follows:

```
dataset_test/
dataset_test/query_graph
dataset_test/data_graph
dataset_test/map_info
dataset_test/partition_point
dataset_test/WG.txt
dataset_test/all_labels.txt
dataset_test/Distance_Matrix.txt
dataset_test/ground_truth.txt
```

The query graph folder stores query graphs in `*.graph`, while the data graph folder stores partition graphs in `*.graph` . `.graph` format is showed as follows:

``` 
t <vertex count> <edge count>
v <vertex ID> <vertex label> <vertex degree>
...
e <src vertex> <dst vertex>
...
```

The map_info folder stores the remap info of partition graphs, which will be used to output results mapping in data graph.

The partition_point folder stores the cut-off vertices of each partition in several text files, each file stores cuf-off vertex ID spliting by space.

WG.txt stores the edge list of data graph, while all_labels.txt stores label information of data graph, in a format of `<vertex ID> <label>` .

Distance_Matrix.txt stores a  $n\times n$ matrix, where $D_{i,j}$ represents the shortest distance between $P_i$ and $P_j$ in Partition Adjacency Graph.

ground_truth.txt stores the prediction results of whether each query appears in each partition, representing by several binary vectors.

### execute

For build, just simply run:

```bash
mkdir build && cd build
cmake ..
make
```

For testing, you can run :`./build/matching/SubgraphMatching.out -q $QUERY_GRAPH_PATH -d $DATA_GRAPH_PATH -k $K -p $NUMBER_OF_PARTITION` .

By default, the program will output Top-k matches along with running time.