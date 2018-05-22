import sys
from turboutils import read_coord


def rules_xenotime(parent_ind, child_ind, coords):
    if parent_ind == child_ind:
        return True
    pname, cname = coords[parent_ind][3], coords[child_ind][3]
    if pname == 'Y':
        return cname == 'O'
    elif pname == 'P':
        return cname == 'O'
    elif pname == 'O':
        return cname == 'Y' or cname == 'P'
    else:
        raise ValueError("atom name is {}".format(coords[parent_ind][3]))


def dec_xenotime(parent_ind, coords):
    if coords[parent_ind][3] == 'P':
        return 0
    else:
        return -1


def fst_sphere(ii, coords, eps, rules):
    def dist2(i):
        return sum((x-x0)**2 for x, x0 in zip(coords[i][:-1], coords[ii][:-1]))
    sorted_inds = sorted(range(len(coords)), key=dist2)
    min_d2 = dist2(sorted_inds[2])
    return set(i for i in range(len(coords)) if dist2(i) <= min_d2*(eps+1e0)
               and rules(ii, i, coords))


def nth_sphere(n, ii, coords, eps=0.3, rules=rules_xenotime, dec=dec_xenotime):
    queue_inds = [ii]
    queue_levs = [n]
    res = set()
    while len(queue_inds) > 0:
        cind, clev = queue_inds.pop(0), queue_levs.pop(0)
        if clev+dec(cind, coords) < 0:
            res.add(cind)
            continue
        fsphere = fst_sphere(cind, coords, eps, rules)
        queue_levs += [clev+dec(cind, coords) for i1 in fsphere
                       if i1 not in queue_inds
                       and i1 not in res and i1 != cind]
        queue_inds += [i1 for i1 in fsphere if i1 not in queue_inds
                       and i1 not in res and i1 != cind]
        res = res.union(fsphere)
    return res


def print_atoms(inds, coords):
    for i in inds:
        print "{:.5f} {:.5f} {:.5f} {}".format(*coords[i])


def center_coords(ii, coords):
    at0 = coords[ii][:]
    return [[x-x0 for x, x0 in zip(at[:-1], at0[:-1])] + [at[-1]]
            for at in coords]


coords_name, atom_num, sphere_n = (sys.argv[1], int(sys.argv[2]),
                                   int(sys.argv[3]))
coords = center_coords(atom_num, read_coord(coords_name))

print_atoms(nth_sphere(sphere_n, atom_num, coords,
                       rules=lambda *args: True, eps=1e-1), coords)
