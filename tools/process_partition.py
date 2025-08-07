import argparse
import os
import numpy as np
import heapq
from multiprocessing import Pool, shared_memory
import multiprocessing as mp
from typing import List, Set, Tuple
import warnings

def CSR_generator(edges, labels, path):
    deg = dict()
    n = 0
    m = len(edges)
    for u, v in edges:
        if u not in deg:
            n += 1
            deg[u] = 0
        deg[u] += 1
        if v not in deg:
            n += 1
            deg[v] = 0
        deg[v] += 1

    with open(path, 'w') as file:
        file.write(f't {n} {m}\n')
        for k, v in deg.items():
            file.write(f'v {k} {labels[k]} {v}\n')
        for u, v in edges:
            file.write(f'e {u} {v}\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', type=str, required=True, help="边文件目录")
    parser.add_argument('-o', type=str, required=True, help="输出目录")
    parser.add_argument('-l', type=str, required=True, help="标签文件")
    args = parser.parse_args()

    n = 0
    mxn = 0
    edge_list = []
    node_list = []

    print("Loading files...")
    for filename in sorted(os.listdir(args.d)):
        file_path = os.path.join(args.d, filename)
        edges = set()
        nodes = set()
        with open(file_path, 'r') as file:
            for line in file:
                parts = line.strip().split()
                if len(parts) != 2:
                    continue
                u, v = int(parts[0]), int(parts[1])
                edges.add((u, v))
                nodes.add(u)
                nodes.add(v)
                mxn = max(mxn, u, v)
        edge_list.append(edges)
        node_list.append(nodes)
        n += 1

    print(f"Loaded {n} partitions, max node ID: {mxn}")
    print("finish load file")

    labels = []
    with open(args.l, 'r') as file:
        for line in file:
            parts = line.strip().split()
            if len(parts) >= 2:
                labels.append(int(parts[1]))

    os.makedirs(f'{args.o}/partitions', exist_ok=True)
    for i in range(n):
        CSR_generator(edge_list[i], labels, f'{args.o}/partitions/{i}.graph')
    print("finish generate CSR")

    rep = [[] for _ in range(mxn + 1)]
    for i in range(n):
        for o in node_list[i]:
            rep[o].append(i)

    os.makedirs(f'{args.o}/node_replica', exist_ok=True)
    for i in range(n):
        with open(f'{args.o}/node_replica/{i}.txt', 'w') as file:
            for o in node_list[i]:
                if len(rep[o]) > 1:
                    file.write(f'{o} ')
            file.write('\n')
    print("finish replica")