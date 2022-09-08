from sklearn.cluster import MiniBatchKMeans
from sklearn.metrics import silhouette_score, silhouette_samples
import pandas as pd

def kmeans_clustering_2d(xy, n_clusters=2, xlabel='x', ylabel='y'):
    '''
    Do KMeans clustering on 2D data and calculate Silhouette Score as quality metric.
    '''

    # Initialize the K-Means model
    kmeans = MiniBatchKMeans(n_clusters = n_clusters)

    # Fitting the model to training set
    kmeans.fit(xy)
    
    # calculate silhouette scores
    sil_avg = silhouette_score(xy, kmeans.labels_) # average score
    sil_samples = silhouette_samples(xy, kmeans.labels_) # score per sample

    # collect results
    results = pd.DataFrame()
    results[xlabel] = xy[:, 0]
    results[ylabel] = xy[:, 1]
    results['labels'] = kmeans.labels_
    results['inertia'] = kmeans.inertia_
    results['silhouette_avg'] = sil_avg
    results['silhouette_samples'] = sil_samples

    return results
