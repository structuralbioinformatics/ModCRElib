import numpy as np
from sklearn.cluster import AgglomerativeClustering

# Example adjacency matrix with probabilities (similarities)
B = np.array([
    [1.0, 0.85, 0.1, 0.2],
    [0.85, 1.0, 0.05, 0.25],
    [0.1, 0.05, 1.0, 0.9],
    [0.2, 0.25, 0.9, 1.0]
])
A = np.array([
    [0.0, 0.85, 1e-10, 0.002],
    [0.85, 0.0, 0.005, 0.25],
    [1e-10, 0.005, 0.0, 0.0009],
    [0.002, 0.25, 0.0009, 0.0]
])

# Convert similarity (probability) matrix A to distance matrix D
# Distance = 1 - similarity probability
# D = 1 - A
D=A

# Perform agglomerative clustering with distance threshold
clustering = AgglomerativeClustering(
    n_clusters=None,  # no fixed number of clusters
    metric='precomputed',  # we provide distance matrix directly
    linkage='average',
    distance_threshold=0.005  # tune this threshold to get clusters
)

labels = clustering.fit_predict(D)
n_clusters = labels.max() + 1
print(f"Number of clusters found: {n_clusters}")

# For each cluster, find the central element minimizing sum of distances
for cluster_id in range(n_clusters):
    members = np.where(labels == cluster_id)[0]
    print(f"Cluster {cluster_id}: members = {members}")

    min_sum_dist = float('inf')
    central_element = None

    for i in members:
        # Sum distance to other cluster members
        sum_dist = sum(D[i, j] for j in members if i != j)

        if sum_dist < min_sum_dist:
            min_sum_dist = sum_dist
            central_element = i

    print(f"  Central element: {central_element} (min total distance = {min_sum_dist:.3f})")
