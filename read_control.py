from parse_coord import read_coord_pat, delta_grad, delta_energy, del_kinds
import sys
coord_pat_filename = sys.argv[1]
control_filename = sys.argv[2]
excls = []
if len(sys.argv)>3:
    excls = [el.lower() for el in sys.argv[3:]]
def read_control(filename, atoms):
    f = open(filename)
    for line in f:
        if "$grad          cartesian gradients" in line:
            break
    en = 0
    grads = []
    while(True):
        line = f.next()
        if '$end' in line:
            f.close()
            return grads, en
        en = float(line.split('=')[2].split()[0])
        for at, line in zip(atoms, f):
            if at[3] != line.split()[3]:
                print at[3], line.split()[3]
        grads = []
        for at, line in zip(atoms, f):
            aux = line.replace('D', 'e').replace('d', 'e').split()[0:3]
            grads += [[float(g) for g in aux]]

ats, qs, dqs, dq_kinds = read_coord_pat(coord_pat_filename)
grads, en = read_control(control_filename, del_kinds(ats))
dgs = delta_grad(ats, qs, dqs)
de = delta_energy(ats, qs, dqs)

print 'Total Energy: {} a.u.'.format(en + de)
print '-'*10
g2_tot =  0e0
for at, g, dg in zip(ats, grads, dgs):
    if at[3].lower() in excls:
        continue
    aux = [gx + dgx for gx,dgx in zip(g, dg)]
    g2_tot += sum( x**2 for x in aux)
    print "{0[0]} {0[1]} {0[2]} {1}".format(aux, at[3])
print '-'*10

print '|grad| = {}'.format(g2_tot**0.5)

