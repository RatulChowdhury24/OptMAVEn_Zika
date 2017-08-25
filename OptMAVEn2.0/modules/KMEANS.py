import os
from collections import OrderedDict

import numpy as np

def kmeans(data, k, c=None, verbose=False):
    centroids = []

    if not isinstance(data, dict):
        data = OrderedDict(enumerate(data))

    centroids = randomize_centroids(data, centroids, k)  

    old_centroids = [[] for i in range(k)]

    iterations = 0
    while not (has_converged(centroids, old_centroids, iterations)):
        iterations += 1

        # assign data points to clusters
        clusters = euclidean_dist(data, centroids)

        # recalculate centroids
        for index, cluster in enumerate(clusters):
            old_centroids[index] = centroids[index]
            centroids[index] = np.mean(cluster.values(), axis=0).tolist()
    
    if verbose:
        for cluster in clusters:
            print(str(len(cluster)) + " members")
            print(np.array(cluster.values()).tolist())
            print("----------\n")

    return clusters
   
def euclidean_dist(data, centroids):
    #mfa#
    k = len(centroids)
    clusters = [OrderedDict() for i in range(k)]
    for label, coords in data.iteritems():
        mu_index = min([(i, np.linalg.norm(coords-centroid)) \
                            for i, centroid in enumerate(centroids)], key=lambda t:t[1])[0]
        clusters[mu_index][label] = coords

    # If any cluster is empty then assign one point
    # from data set randomly so as to not have empty
    # clusters and 0 means.
     
    for cluster in clusters:
        if len(cluster) == 0:
            #mfa#
            label, coords = data.items()[np.random.randint(0, len(data), size=1)]
            cluster[label] = coords
            
    return clusters

def randomize_centroids(data, centroids, k):
    for cluster in range(0, k):
        centroids.append(data.values()[np.random.randint(0, len(data), size=1)].flatten().tolist())
    return centroids

 
def has_converged(centroids, old_centroids, iterations):
    MAX_ITERATIONS = 1000
    if iterations > MAX_ITERATIONS:
        return True
    return old_centroids == centroids
