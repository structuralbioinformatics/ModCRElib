import numpy as np
import networkx as nx
from networkx.algorithms.community import girvan_newman

# Example adjacency matrix of probabilities (similarities)
B = np.array([
    [1.0, 0.85, 0.1, 0.2],
    [0.85, 1.0, 0.05, 0.25],
    [0.1, 0.05, 1.0, 0.9],
    [0.2, 0.25, 0.9, 1.0]
])
C = np.array([
    [0.0, 0.85, 1e-10, 0.002],
    [0.85, 0.0, 0.005, 0.25],
    [1e-10, 0.005, 0.0, 0.0009],
    [0.002, 0.25, 0.0009, 0.0]
])

#A=1-C

# Small epsilon to avoid log(0)
epsilon = 1e-10

# Replace zeros or values below epsilon with epsilon
A_safe = np.clip(C, epsilon, 1.0)

# Compute distance matrix as -log(probabilities)
D = -np.log(A_safe)
D_min = D.min()
D_max = D.max()

A = (D - D_min) / (D_max - D_min)

print(A)

# Build graph from adjacency matrix
G = nx.Graph()
n = A.shape[0]

# Add edges for similarity above threshold (to reduce noise)
threshold = 0.1
for i in range(n):
    for j in range(i + 1, n):
        if A[i, j] > threshold:
            G.add_edge(i, j, weight=(A[i, j]))

# Girvan-Newman community detection generator
comp_generator = girvan_newman(G)

# Choose number of communities to extract from hierarchy (e.g., 2 or 3)
num_communities = 2
limited = None

for communities in comp_generator:
    limited = communities
    if len(communities) == num_communities:
        break

clusters = [list(c) for c in limited]
print(f"Detected {len(clusters)} clusters")

# Calculate central element in each cluster by minimal sum distance
# Distance = 1 - similarity (0 if i == j)
for idx, cluster in enumerate(clusters):
    print(f"Cluster {idx}: members = {cluster}")

    min_sum_dist = float('inf')
    central_element = None

    for i in cluster:
        sum_dist = 0
        for j in cluster:
            if i != j:
                # Use similarity from A if edge exists, else 0
                sim = A[i, j] if A[i, j] > 0 else 0
                sum_dist += 1 - sim

        if sum_dist < min_sum_dist:
            min_sum_dist = sum_dist
            central_element = i

    print(f"  Central element: {central_element} (min total distance = {min_sum_dist:.3f})")
