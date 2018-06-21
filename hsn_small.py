import sys
import scipy as sc
from time import sleep
from read_control import read_control
import subprocess as sp
coord = 'coord'
control = 'control'
realmos = 'realmos'
imagmos = 'imagmos'
mos = 'mos'
basis = 'basis'
embedding = 'embedding'
coord_pat = 'coord_pat'
excls = ' y o2 o10 zz3 zz4'.split()
node_dir_prefix = 'node'
srun_command = ''
#srun_command = 'srun -N1 -n24'
dscf_command = 'dscf>out'
grad_command = 'grad>grad.out'
hessian_file = 'HESS.DAT'
eps = 1e-1
delay = 10
strict = True
nrow = 5
fmt_float = "{:.10f}"
fmt_sci = "{:25.17e}"
#
number_of_nodes = int(sys.argv[1])
calc_files = [coord, control, realmos, imagmos, basis, embedding]


class Node(object):
    IDLE, RUN = range(2)

    def __init__(self, calc_files, node_dir):
        self.calc_files = calc_files[:]
        self.node_dir = node_dir
        for x in calc_files:
            sp.call("cp {0} {1}/{0}".format(x, node_dir), shell=True)
        self.process = None
        self.occupied = None
        self.status = Node.IDLE

    def run(self, command):
        self.process = sp.Popen(command, shell=True, cwd=self.node_dir)
        self.status = Node.RUN
        self.command = command

    def run_block(self, command):
        sp.call(command, shell=True, cwd=self.node_dir)

    def save_files(self):
        if self.status == Node.RUN:
            print "node {} is busy".format(self.node_dir)
            return -1
        print "save files in {}".format(self.node_dir)
        for x in calc_files:
            self.run('cp {0} {0}0'.format(x))
        return 0

    def restore_files(self):
        if self.status == Node.RUN:
            print "node {} is busy".format(self.node_dir)
            return -1
        print "restore files in {}".format(self.node_dir)
        for x in calc_files:
            self.run('mv {0}0 {0}'.format(x))
        return 0

    def update_status(self):
        if self.process is None:
            self.status = Node.IDLE
            return self.status
        if self.process.poll() is None:
            self.status = Node.RUN
            return self.status
        if self.process.poll() == 0:
            self.status = Node.IDLE
            return self.status
        else:
            sys.stderr.write("Something wrong\
on node {}\n".format(self.node_dir))
            sys.stderr.write("last command: {}".format(self.command))
            self.status = Node.IDLE
            return self.status

calc_files = [coord, control, realmos, imagmos, basis, embedding]


def make_indexes(coord_pat=coord_pat, excls=excls):
    labels = []
    atom_coords = []
    excl_inds = []
    f = open(coord_pat)
    header = f.next()
    if '$coord' not in header:
        raise ValueError('Invalid coord_pat file')
    for n, line in enumerate(f):
        if '$end' in line:
            return labels, atom_coords, excl_inds

        aux = line.split()
        if aux[3] not in excls:
            labels += [aux[3]]*3
            atom_coords += [(n, i) for i in range(3)]
        else:
            excl_inds += [n]


def exclude_grads_and_flatten(grads, excl_inds):
    res = [g for (n, g) in enumerate(grads) if n not in excl_inds]
    res = sc.array(res)
    return sc.reshape(res, 3*len(res))


def update_coord(i, eps, atom_coords, coord):
    nat, ix = atom_coords[i]
    changed = []
    print "updating {}-th coordinate by {}".format(i, eps)
    f = open(coord)
    header = f.next()
    if '$coord' not in header:
        raise ValueError('Invalid coord file {}'.format(coord))
    new_coord = header
    for n, line in enumerate(f):
        if n != nat:
            new_coord += line
            changed += [False]
            continue
        aux = line.split()
        aux[:3] = map(float, aux[:3])
        aux[ix] += eps
        aux[:3] = map(lambda x: "{:.8f}".format(x), aux[:3])
        new_coord += " ".join(aux)+'\n'
        changed += [True]
    f.close()
    f = open(coord, 'w')
    f.write(new_coord)
    f.close()
    return changed


def dgdxi(i, node, atom_coords, excl_inds,
          control=control, coord=coord,
          eps=eps, srun=srun_command, dscf=dscf_command,
          grad=grad_command):
    if node.occupied is None:
        node.occupied = i
    # wait
    while(node.occupied != i):
        yield None
    # wait
    while(node.update_status() != Node.IDLE):
        yield None
    my_coord   = node.node_dir+'/'+coord
    my_control = node.node_dir+'/'+control
    node.save_files()
    # wait
    while(node.update_status() != Node.IDLE):
        yield None
    changed = update_coord(i, eps, atom_coords, my_coord)
    print "node {} coord {} dscf+".format(node.node_dir, i)
    node.run(srun + dscf)
    while(node.update_status() != Node.IDLE):
        yield None
    print "node {} coord {} grad+".format(node.node_dir, i)
    node.run(srun + grad)
    # wait
    while(node.update_status() != Node.IDLE):
        yield None
    res = exclude_grads_and_flatten(read_control(my_control, changed)[0],
                                    excl_inds)

    update_coord(i, -2*eps, atom_coords)
    print "node {} coord {} dscf-".format(node.node_dir, i)
    node.run(srun + dscf)
    # wait
    while(node.update_status() != Node.IDLE):
        yield None
    print "node {} coord {} grad-".format(node.node_dir, i)
    node.run(srun + grad)
    # wait
    while(node.update_status() != Node.IDLE):
        yield None
    res = res - exclude_grads_and_flatten(read_control(my_control, changed)[0],
                                          excl_inds)
    res = res/(2*eps)
    node.occupied = None
    node.restore_files()
    yield res


labels, atom_coords, excl_inds = make_indexes()
N = len(labels)
print "There are {} coordinates".format(N)
nodes = [Node(node_dir_prefix+str(i), calc_files, strict=strict)
         for i in range(number_of_nodes)]
finished = [False] * N
tasks = [dgdxi(i, nodes[i % number_of_nodes], atom_coords, excl_inds) for i in
         range(N)]
hessian = sc.zeros((N, N))
while not all(finished):
    for ii, task in enumerate(tasks):
        res = task.next()
        if res is not None:
            hessian[ii, :] = res[:]
            finished[ii] = True
    sleep(delay)

print "HESSIAN stored in {}".format(hessian_file)
print "-"*20
hess_f = open(hessian_file, 'w')
delta = 0e0
fmt = (fmt_float+" \t ") * N
for ii in xrange(N):
    print fmt.format(*hessian[ii, :])
    for jj in xrange(ii):
        delta = delta + (hessian[ii, jj] - hessian[jj, ii])**2

print "delta = {}".format(delta)
hessian = (hessian + sc.transpose(hessian))*0.5
fmt = "{:2d}{:3d} " + fmt_sci*nrow+'\n'
for ii in xrange(N):
    for coln in xrange(N/nrow-1):
        hess_f.write(fmt.format(ii+1, coln+1,
                                *hessian[ii, nrow*coln:nrow*(coln+1)]))
    if coln is None:
        coln = 0
    else:
        coln += 1
    nnrow = N - (N/nrow)*nrow
    fmt = "{:2d}{:3d} " + fmt_sci*nnrow+'\n'
    hess_f.write(fmt.format(ii+1, coln+1, *hessian[ii, N-nnrow:]))
