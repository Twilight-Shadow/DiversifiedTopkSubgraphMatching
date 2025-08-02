nodes = [set() for _ in range(8000)]
avg_size = 0
for i in range(98):
    with open("./dataset/" + str(i) + ".txt", 'r') as file:
        for line in file:
            parts = line.strip().split()
            u = int(parts[0])
            v = int(parts[1])
            nodes[i].add(u)
            nodes[i].add(v)
    avg_size += len(nodes[i])
# print(avg_size / 98)
matrix = [[] for _ in range(8000)]
for i in range(98):
    for j in range(98):
        if  nodes[i].intersection(nodes[j]):
            matrix[i].append(1)
        else:
            matrix[i].append(0)
with open("./Partition_Matrix.txt", 'w') as file:
    for i in range(98):
        for j in range(98):
            file.write(str(matrix[i][j]) + ' ')
        file.write('\n')
