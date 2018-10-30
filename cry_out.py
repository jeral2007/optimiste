import scipy as sc
import scipy.linalg as lng


def read_header(f):
    state = 0
    for line in f:
        if 'EEEEEEEEEE STARTING  DATE' not in line:
            continue
        state = 1
        break
    if state == 0:
        raise ValueError('Incorrect output stream')
    return f.next()


def read_type(f):
    for line in f:
        if 'CALCULATION' in line:
            return line.split()[0]

    raise ValueError('Incorrect output stream')


def read_space_group(f):
    for line in f:
        if 'SPACE GROUP' in line:
            return line.split(':')[1][2:]
    raise ValueError('Incorrect output stream')


def read_lat_pars(f):
    state = 0
    for line in f:
        if 'LATTICE PARAMETERS  (ANGSTROMS AND DEGREES)' in line:
            state = 1
            break
    if state == 0:
        raise ValueError('Incorrect output stream')
    f.next()
    par_names = 'a b c alpha beta gamma'.split()
    res = {k: float(v) for k, v in
           zip(par_names, f.next().split())}
    return res


def read_atoms_coord(f):
    state = 0
    for line in f:
        if 'NUMBER OF IRREDUCIBLE ATOMS' in line:
            state = 1
            n_atoms = int(line.split(':')[1])
            break
    if state == 0:
        raise ValueError('Incorrect output stream')
    for line in f:
        if 'ATOM AT. N.' in line:
            state = 2
            break
    if state == 1:
        raise ValueError('Incorrect output stream')
    coords, types = [], []
    for line, ia in zip(f, xrange(n_atoms)):
        aux = line.split()
        types += [int(aux[1])]
        coords += [sc.array([float(x) for x in aux[2:]])]
    if '*****************' not in f.next():
        raise ValueError('Incorrect output stream')
    return coords, types


def read_grads(f):
    state = 0
    for line in f:
        if 'CARTESIAN FORCES' in line:
            state = 1
            break
    if state == 0:
        raise ValueError('Incorrect output stream')
    f.next()
    res = []
    for line in f:
        aux = line.split()
        if len(aux) == 0:
            return sc.array(res).transpose()
        res += [map(float, aux[-3:])]


def read_symops(f):
    state = 0
    for line in f:
        if 'TRANSLATORS IN FRACTIONAL UNITS' in line:
            state = 1
            n_syms = int(line.split()[1])
            break
    if state == 0:
        raise ValueError('Incorrect output stream')
    f.next()
    f.next()
    res = []
    for isym, line in zip(xrange(n_syms), f):
        aux = map(float, line.split())
        rotm = sc.array(aux[2:11])
        tranv = sc.array(aux[11:])
        res += [(rotm.reshape((3, 3)), tranv)]
    return res


def read_lat_mat(f):
    state = 0
    for line in f:
        if 'DIRECT LATTICE VECTORS' in line:
            state = 1
            break
    if state == 0:
        raise ValueError('Incorrect output stream')
    f.next()
    res = []
    res += [map(float, f.next().split())]
    res += [map(float, f.next().split())]
    res += [map(float, f.next().split())]
    return sc.array(res)


class CrystallOut(object):
    def assert_primitive(self):
        if self.primitive_cell_atoms is None:
            raise ValueError('primitive cell params is undefined\
 call frac_coords first!')

    def __init__(self, output):
        f = open(output)
        self.header = read_header(f)
        self.type_of_calc = read_type(f)
        self.space_group = read_space_group(f)
        if self.type_of_calc == 'CRYSTAL':
            self.lat_pars = read_lat_pars(f)
        else:
            self.lat_pars = None
        self.frac_coords, self.at_types = read_atoms_coord(f)
        self.sym_ops = read_symops(f)
        if self.type_of_calc == 'CRYSTAL':
            self.lat_mat = read_lat_mat(f)
        else:
            self.lat_pars = None
        self.f = f
        self.grads = None
        self.primitive_cells_atoms = None

    def primitive_cell(self):
        res1, res2 = [], []
        for cs, at in zip(self.frac_coords, self.at_types):
            for rotm, v in self.sym_ops:
                res1 += [sc.dot(rotm, cs)+v]
                res2 += [at]
        self.primitive_cell_atoms = res2
        self.primitive_cell_coords = res1
        return res1, res2

    def cart_grads(self):
        self.assert_primitive()
        if self.grads is None:
            self.grads = read_grads(self.f)
        return self.grads

    def frac_grads(self):
        self.assert_primitive()
        if self.grads is None:
            self.grads = read_grads(self.f)
        return sc.dot(lng.inv(self.lat_mat), self.grads)


if __name__ == '__main__':
    import sys
    outfile = sys.argv[1]
    dr = float(sys.argv[2])
    cry_out = CrystallOut(outfile)
    print cry_out.header
    print cry_out.type_of_calc
    print cry_out.lat_mat
    print cry_out.frac_coords, cry_out.at_types

    print '-'*20
    for cs, at in zip(*cry_out.primitive_cell()):
        print "{0}   {1[0]:.5f} {1[1]:.5f} {1[2]:.5f}".format(at, cs)
    print '-'*20
    frac_grads = cry_out.frac_grads()
    for ii in xrange(len(cry_out.primitive_cell_atoms)):  #  FIXME
        print "{0}   {1[0]:.5f} {1[1]:.5f} {1[2]:.5f}".format(
            cry_out.primitive_cell_atoms[ii], frac_grads[:, ii])
        print "{0}   {1[0]:.5f} {1[1]:.5f} {1[2]:.5f}".format(
            cry_out.primitive_cell_atoms[ii], cry_out.grads[:, ii])
        print "recommended charge position and value"
        pos = cry_out.primitive_cell_coords[ii]
        eg = cry_out.grads[:, ii]
        eg_norm = sc.dot(eg, eg)**0.5
        eg = eg/sc.dot(eg, eg)**0.5
        pos = pos + sc.dot(lng.inv(cry_out.lat_mat), eg*dr)
        print pos, eg_norm**0.5*dr**2
