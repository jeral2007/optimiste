from parse_coord import read_coord_pat, delta_grad, delta_energy, del_kinds
import sys
from scipy import float64
import scipy as sc


def read_control(filename, atoms):
    f = open(filename)
    for line in f:
        if "$grad          cartesian gradients" in line:
            break
    en = 0
    grads = []
    first_step = True
    while(True):
        line = f.next()
        if not first_step and '$last' in line or '$actual' in line:
            f.close()
            return grads, en
        if '$end' in line:
            f.close()
            return grads, en
        if 'maximum norm' in line:
            f.close()
            return grads, en
        en = float64(line.split('=')[2].split()[0])
        for at, line in zip(atoms, f):
            pass
    #        if at[3] != line.split()[3]:
    #               print at[3], line.split()[3]
        grads = []
        for at, line in zip(atoms, f):
            aux = line.replace('D', 'e').replace('d', 'e').split()[0:3]
            grads += [[float64(g) for g in aux]]
        first_step = False


def print_grads(aux, at):
    print "{0[0]} {0[1]} {0[2]} {1} {2}".format(aux,
                                                sum(c**2
                                                    for c in aux[:3])**0.5,
                                                at[3])


def take_grads(coord_pat_filename, control_filename,
               excls, callback=print_grads):
    ats, qs, dqs, dq_kinds = read_coord_pat(coord_pat_filename)
    grads, en = read_control(control_filename, del_kinds(ats))
    dgs = delta_grad(ats, qs, dqs)
    de = delta_energy(ats, qs, dqs)

    print 'Total Energy: {} a.u.'.format(en + de)
    print '-'*10
    g2_tot = 0e0
    for at, g, dg in zip(ats, grads, dgs):
        if at[3].lower() in excls:
            continue
        aux = [gx + dgx for gx, dgx in zip(g, dg)]
        callback(aux, at)
        g2_tot += sum(x**2 for x in aux)
    print '-'*10

    print '|grad| = {}'.format(g2_tot**0.5)
    return g2_tot**0.5

class Grads(object):

    def read_coord_pat(self):
        self.labels = []
        self.atom_coords = []
        self.excl_inds = []
        self.atom_names = []
        f = open(self.coord_pat)
        header = f.next()
        if '$coord' not in header:
            raise ValueError('Invalid coord_pat file')
        for n, line in enumerate(f):
            if '$end' in line:
                return
            aux = line.split()
            self.atom_names += [aux[3]]
            if aux[3] not in self.excls:
                self.labels += [aux[3]]*3
                self.atom_coords += [(n, i) for i in range(3)]
            else:
                self.excl_inds += [n]

    def flatten_grads(self):
        res = [g for (n, g) in enumerate(self.grads)
               if n not in self.excl_inds]
        res = sc.array(res, dtype=sc.float64) # important
        return sc.reshape(res, 3*len(res)) # FIXME recalculation!

    def __init__(self, excls, coord_pat, control):
        self.excls = excls
        self.coord_pat = coord_pat
        self.control = control
        self.read_coord_pat()
        self.grads, self.en = read_control(control, self.atom_names)


if __name__=='__main__':
    coord_pat_filename = sys.argv[1]
    control_filename = sys.argv[2]
    excls = []
    if len(sys.argv)>3:
       excls = [el.lower() for el in sys.argv[3:]]
    take_grads(coord_pat_filename, control_filename, excls)

