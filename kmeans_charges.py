import sys
from turboutils import read_coord

def intervals(arr):
    if len(arr) == 0:
        return []
    cur_int = [arr[0]]
    res_aux = [cur_int]
    for x in arr[1:]:
        if x-cur_int[-1] == 1:
            cur_int += [x]
        else:
            cur_int = [x]
            res_aux += [cur_int]

    def lambda1(c):
        if len(c) == 1:
            return "{}".format(c[0])
        elif len(c) == 2:
            return "{},{}".format(*c)
        else:
            return "{}-{}".format(c[0], c[-1])

    return ",".join(lambda1(c) for c in res_aux)


def read_charges(charges_filename):
    err = lambda s: ValueError("Incorrect charge file {}, \
                               invalid string: {}".format(charges_filename, s))
    f = open(charges_filename)
    header = f.next()
    if '=======' not in header:
        raise err(header)
    header = f.next()
    if 'Valence' not in header:
        raise err(header)
    header = f.next()
    if '=======' not in header:
        raise err(header)
    res = []
    for line in f:
        if '============' in line:
            return res
        aux = line.split()
        if len(aux) != 5:
            raise err(line)
        try:
            res += [[aux[0]] + [float(x) for x in aux[1:]]]
        except Exception:
            raise err(line)


def kmeans(arr, N, maxiter=300):
    import random
    clusters = [[] for i in range(N)]

    def mean(cluster):
        if len(cluster) == 0:
            return 0
        return sum(arr[i] for i in cluster)/len(cluster)

    cc = 0
    rnd = range(len(arr))
    random.shuffle(rnd)
    for i in rnd:
        if cc == N:
            cc = 0
        clusters[cc] += [i]
        cc += 1
    means = [mean(cluster) for cluster in clusters]

    def nearest_cluster(x):
        return min(zip(means, clusters), key=lambda t: abs(t[0] - x))
    # until not changed
    iter = 1
    while(True):
        n_trans = 0
        if iter >= maxiter:
            break
        print "Iteration {}".format(iter)
        iter += 1
        for cluster, avg in zip(clusters, means):
            ii = 0
            while(ii < len(cluster)):
                x = arr[cluster[ii]]
                d, cluster2 = nearest_cluster(x)
                if abs(x-d) < abs(x - avg):
                    cluster2 += [cluster.pop(ii)]
                    n_trans += 1
                else:
                    ii += 1
        means = [mean(cluster) for cluster in clusters]
        if n_trans == 0:
            break
    return [(m, c) for m, c in zip(means, clusters) if len(c) > 0]


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


charges = read_charges(sys.argv[1])
qs = [q[2] for q in charges]
sigma = 0e0
# for a, c in kmeans(qs, int(sys.argv[2]), maxiter=100):
res = per_clusters(qs, 1e0/float(sys.argv[2]))
print "=="*40
print "{}  clusters".format(len(res))
print "of following sizes:"
print ", ".join(str(len(ac[1])) for ac in res)
print "=="*40
avs = []
for a, c in res:
    av, s_c = stat_cluster(c, qs)
    avs += [av]
    sigma += s_c
    print '='*10
    print "AVERAGE:{} SIZE: {} TOTAL: {}".format(av, len(c), av*len(c))
    print '='*10
    for i in sorted(c):
        print "{} {} {}".format(i, charges[i][0], charges[i][2])
print "SIGMA: {}".format(sigma)

ats = read_coord(sys.argv[3])
print "===="*10
print "$coord"

not_first_occ = []
for ind, at in enumerate(ats):
    if charges[ind][0].lower() != at[-1].lower():
        raise ValueError("Atom no. {} in charge is {} \
while in coord this atom is {}".format(ind,
                                       charges[ind][0].lower(), at[-1].lower()))
    clust_no = [ind in c for avv, c in res].index(True)
    if clust_no not in not_first_occ:
        print "{} {} {} {}{c_no:} {av:.5f}".format(*at, c_no=clust_no,
                                                     av=res[clust_no][0])
        not_first_occ += [clust_no]
    else:
        print "{} {} {} {}{c_no}".format(*at, c_no=clust_no)
print "$end"
print "=="*10
print "CLUSTER NUMBERS for control file"
for n, (av, c) in enumerate(res):
    print "{}{} {}".format(charges[c[0]][0], n, intervals([i+1 for i in c]))
