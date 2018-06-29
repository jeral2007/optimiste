from stat_proc import *
import standard_config as sconf

#config
coord_pat = 'coord_pat'
coord = 'coord'
control = 'control'
excls = ['zz3', 'zz4']
bonds_on = ['o10', 'o2']
tol = 0.05
neighbours_count=8
def read_coord_and_pat(coord=coord, coord_pat=coord_pat):
    f_coord = open(coord)
    f_pat = open(coord_pat)
    header_pat = f_pat.next()
    header_coord = f_coord.next()
    if '$coord' not in header_pat or '$coord' not in header_coord:
       raise  ValueError('invalid coord or coord_pat file')
    coords = []
    kinds = []
    coords_pat = []
    for line_coord, line_pat in zip(f_coord, f_pat):
        end_of_coord = ['$end' in line for line in [line_coord, line_pat]]
        if any(end_of_coord):
           if not all(end_of_coord):
              raise ValueError('Inconsistent coord and coord_pat files (different number of atoms)')
           else:
              return coords, coords_pat, kinds
        coords += [map(float,line_coord.split()[:3])]
        kinds += [line_pat.split()[3]]
        coords_pat += [map(float,line_pat.split()[:3])]
    raise ValueError("no $end in coord or coord_pat")


def calc_bonds(to_calc_kind, kinds_coords, neighbours_count=10):
     kinds, coords = kinds_coords
     inds = [i for (i,at_kind)
               in enumerate(kinds) 
               if at_kind == to_calc_kind]
     bonds = []
     bond_kinds = []
     def dist(ii, jj):
         return sum((x1-x2)**2 for x1,x2 in zip(coords[ii], coords[jj]))**0.5
     for ii in inds:
         aux = [dist(ii, jj) for jj in xrange(len(coords))]
         aux_inds = sorted(range(len(aux)), key=lambda ii:aux[ii])[1:neighbours_count+1]
         bonds +=[aux[jj] for jj in aux_inds]
         bond_kinds +=["{} - {}".format(to_calc_kind, kinds[jj]) for jj in aux_inds]
     return bonds, bond_kinds
     

coords, coords_pat, kinds = read_coord_and_pat()
print '-'*20
print "offsets"
print '-'*20
for c1, c2, k in zip(coords, coords_pat, kinds):
    if k in excls:
       continue
    print "{} {} {}".format(*[x2-x1 for (x2,x1) in zip(c2, c1)])
print '-'*20
print 'bonds'

print '-'*20
bonds = []
b_kinds = []
bonds_pat= []
b_kinds_pat=[]
for at in bonds_on:
    tmp_b, tmp_k = calc_bonds(at, (kinds, coords), neighbours_count)
    tmp_b_pat, tmp_k_pat = calc_bonds(at, (kinds, coords_pat), neighbours_count)
    bonds += tmp_b
    bonds_pat += tmp_b_pat
    b_kinds += tmp_k
    b_kinds_pat +=tmp_k_pat
mclusters = sorted(per_clusters(bonds, tol), key = lambda a:a[0])
mclusters_pat = sorted(per_clusters(bonds_pat, tol), key = lambda a:a[0])
print "Bond \t <l> \t sqrt(<sigma>)"
for m, c in mclusters_pat:
    length, sigma = stat_cluster(c, bonds_pat)
    bk = set(b_kinds_pat[i] for i in c)
    for b in bk:
        print "{} \t {:.5f} \t {:.5f}".format(b, length*sconf.au2ang, sconf.au2ang*(sigma+1e-12)**0.5)

print "Bond \t <l> \t sqrt(<sigma>)"
for m, c in mclusters:
    length, sigma = stat_cluster(c, bonds)
    bk = set(b_kinds[i] for i in c)
    for b in bk:
        print "{} \t {:.5f} \t {:.5f}".format(b, length*sconf.au2ang, sconf.au2ang*sigma**0.5)

print '-'*20
print "detailed"
print '-'*20

for m, c in mclusters:
    print "cluster"
    print "-"*10
    for i in c:
        print "{} {} {}".format(i,bonds[i]*sconf.au2ang, b_kinds[i])
