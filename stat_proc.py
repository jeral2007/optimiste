def per_clusters(arr, eps):
    processed = []
    clusters = []
    while(len(processed) < len(arr)):
        for ii in range(len(arr)):
            if ii not in processed:
                cur_cluster = [ii]
                processed += [ii]
                break
        transfers = 1
        while(transfers > 0):
            to_add, acc = [], []
            for ii in cur_cluster:
                to_add = [i for i, c in enumerate(arr) if abs(arr[ii]-c) < eps
                          and i not in to_add and i not in processed]
                acc += to_add
                processed += to_add

            cur_cluster += acc
            transfers = len(to_add)
        clusters += [cur_cluster]

    return [(sum(arr[i] for i in c)/len(c), c) for c in clusters]


def stat_cluster(cluster, arr):
    N = len(cluster)
    if N == 0:
        return 0e0, 0e0
    m = sum(arr[i] for i in cluster)/N
    d = sum(arr[i]**2 for i in cluster)/N
    return m, d-m**2
