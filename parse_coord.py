import sys

def del_kinds(atoms):
    return [ at[:3]+[''.join(c for c in at[3] if not c.isdigit())
                     for at in atoms]

def write_coord(atoms):
    f = open('coord', 'w')
    f.write('$coord\n')
    for at in del_kinds(atoms):
        f.write("{0[0]:.6f}   {0[1]:.6f}  {0[2]:.6f}  {0[3]}\n".format(at))
    f.write('$end\n')


def read_coord_pat(filename):
    f = open(filename, 'r')
    first_line = f.next()
    assert('$coord' in first_line)
    atoms = []
    charges = []
    delta_charges = {}
    q_vs_name = {}
    dq_vs_name = {}
    i = 0
    for line in f:
        if line[0] == '#':
            continue
        if '$end' in line:
            return atoms, charges, delta_charges, dq_vs_name
        aux = line.split()
        assert(len(aux) > 3)
        aux[:3] = map(float, aux[:3])
        atoms += [aux[:4]]
        if aux[3] not in q_vs_name.keys():
            if len(aux) < 5:
                raise ValueError("""absent charge definition
                                 at line {} of {}""".format(i+2, filename))
            q = int(aux[4])
            q_vs_name[aux[3]] = q
        else:
            if len(aux) > 4:
                raise ValueError("""multiple charge definition
                                 at line {} of {}""".format(i+2, filename))
            q = q_vs_name[aux[3]]
        charges += [q]

        if aux[3] not in dq_vs_name.keys():
            if len(aux) > 5:
                dq = float(aux[5])
                dq_vs_name[aux[3]] = dq
                delta_charges[i] = dq
        else:
            if len(aux) > 5:
                raise ValueError("""multiple delta charge definition
                                 at line {} of {}""".format(i+2, filename))
            delta_charges[i] = dq_vs_name[aux[3]]
        i += 1

def delta_energy(atoms, qs, dqs):
    from math import sqrt
    res = 0e0
    for i in dqs:
        for j in xrange(i):
            dist = sqrt(sum((atoms[i][k]-atoms[j][k])**2 for k in xrange(3)))
            res += dqs[i]*qs[j]/dist
            if j in dqs:
                res += dqs[i]*dqs[j]/dist
    return res



def delta_grad(atoms, qs, dqs):
    from math import sqrt
    res = []
    for j in xrange(len(atoms)):
        aux = [0e0, 0e0, 0e0]
        q = qs[j]
        if j in dqs:
            q +=dqs[j]
            for i in xrange(len(qs)):
                if i == j:
                    continue
                dist3=sqrt(sum((atoms[i][k]-atoms[j][k])**2
                               for k in xrange(3)))**3
                aux = [g0 + qs[i]*q/dist3*(x-y)
                    for g0, x, y in zip(aux, atoms[i][0:3], atoms[j][0:3])]
        for i in dqs:
            if i == j:
                continue
            dist3=sqrt(sum((atoms[i][k]-atoms[j][k])**2 for k in xrange(3)))**3
            aux = [g0 + dqs[i]*q/dist3*(x-y)
                   for g0, x, y in zip(aux, atoms[i][0:3], atoms[j][0:3])]

        res += [aux]
    return res

def write_basis(basis_pat, dq_vs_name, eps=1e-8):
    f = open(basis_pat,'r')
    out = open('basis', 'w')
    for line in f:
        if '!-dq' not in line:
            out.write(line)
        else:
            atom_name = line.split()[1]
            out.write("          {}         1         {:.8f}\n".format(
                -dq_vs_name[atom_name], eps))


def main():
    res = read_coord_pat(sys.argv[1])
    atoms=res[0]
    print res[-1]

    print "delta_e = {}".format(delta_energy(*res[:-1]))

    print '-'*10
    for at, g in zip(atoms, delta_grad(*res[:-1])):
        print " ".join(str(ga) for ga in g)


    write_basis(sys.argv[2], res[-1])
    write_coord(atoms)

if __name__=='__main__':
    main()

