import numpy as np
import networkx as nx
import community as community_louvain  # python-louvain package

# Your adjacency matrix (similarity)
A = np.array([
    [1.0, 0.8, 0.1, 0.2],
    [0.8, 1.0, 0.0, 0.3],
    [0.1, 0.0, 1.0, 0.9],
    [0.2, 0.3, 0.9, 1.0]
])
B = np.array([
    [0.0, 0.85, 1e-10, 0.002],
    [0.85, 0.0, 0.05, 0.25],
    [1e-10, 0.05, 0.0, 0.0009],
    [0.002, 0.25, 0.0009, 0.0]
])
C=1-B
print(C)

# Small epsilon to avoid log(0)
epsilon = 1e-10

# Replace zeros or values below epsilon with epsilon
A_safe = np.clip(B, epsilon, 1.0)

# Compute distance matrix as -log(probabilities)
D = -np.log(A_safe)
D_min = D.min()
D_max = D.max()

A = (D - D_min) / (D_max - D_min)

print(A)

# Step 1: Convert adjacency matrix to NetworkX graph
# Threshold similarity to build edges (e.g., > 0.1 to avoid too many edges)
threshold = 0.1
G = nx.Graph()
n = A.shape[0]

for i in range(n):
    for j in range(i + 1, n):
        if A[i, j] > threshold:
            G.add_edge(i, j, weight=A[i, j])

# Step 2: Apply Louvain method for community detection
partition = community_louvain.best_partition(G, weight='weight')

# Organize nodes by cluster
clusters = {}
for node, comm_id in partition.items():
    clusters.setdefault(comm_id, []).append(node)

# Step 3: Find central element for each cluster
for comm_id, members in clusters.items():
    print(f"Cluster {comm_id}: members = {members}")

    min_sum_dist = float('inf')
    central_element = None

    for i in members:
        # sum distances to others = sum(1 - similarity)
        sum_dist = 0
        for j in members:
            if i != j:
                # If edge doesn't exist, similarity = 0
                sim = A[i, j] if A[i, j] > 0 else 0
                sum_dist += 1 - sim

        if sum_dist < min_sum_dist:
            min_sum_dist = sum_dist
            central_element = i

    print(f"  Central element: {central_element} (min total distance = {min_sum_dist:.3f})")
