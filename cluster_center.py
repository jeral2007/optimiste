from  turboutils import read_coord
import sys
atoms = read_coord(sys.argv[1])

bx = -1e6, 1e6
by = -1e6, 1e6
bz = -1e6, 1e6

for x, y, z, at in atoms:
    if x<bx[0]:
        bx[0] = x
    if x>bx[1]:
        bx[1] = x
    if y<by[0]:
        by[0] = y
    if y>by[1]:
        by[1] = y
    if z<bz[0]:
        bz[0] = z
    if z>bz[1]:
        bz[1] = z

center = [(a+b)/2 for a,b in [bx, by, bz]]

def dist2(ii, atoms=atoms, center=center):
    return sum((x-x0)**2 for x, x0 in zip(atoms[ii][:-1], center))

sorted_inds = sorted(range(len(atoms)), key=dist2)
for i in sorted_inds:
    print "{:3d}  {:.5f} {:.5f} {:.5f} {}".format(i, *atoms[i])

