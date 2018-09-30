import sys
import scipy as sc
<<<<<<< HEAD
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
=======
import scipy.linalg as scl
from scipy import float64, sqrt # sort of stupid
from time import sleep
from read_control import read_control
from standard_config import *
import subprocess as sp
excls = 'zz6 zz7'.split()
node_dir_prefix = 'node'
srun_command = 'srun -N1 -n1'

hessian_file = 'HESS.DAT'
omega_file = 'OMEGA.DAT.scipy'


eps = 1e-1
ncol = 5

srun_command +=' '#
number_of_nodes = int(sys.argv[1])
eps = float64(sys.argv[2])
>>>>>>> 109f3ccc1a8e0b48e3c0512753b53bb04da77d21
calc_files = [coord, control, realmos, imagmos, basis, embedding]


class Node(object):
    IDLE, RUN = range(2)

<<<<<<< HEAD
    def __init__(self, calc_files, node_dir):
        self.calc_files = calc_files[:]
        self.node_dir = node_dir
=======
    def __init__(self,node_dir, calc_files):
        self.calc_files = calc_files[:]
        self.node_dir = node_dir
        sp.call("mkdir {}".format(node_dir), shell=True)
>>>>>>> 109f3ccc1a8e0b48e3c0512753b53bb04da77d21
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
<<<<<<< HEAD
=======
            sleep(0.5)
>>>>>>> 109f3ccc1a8e0b48e3c0512753b53bb04da77d21
        return 0

    def restore_files(self):
        if self.status == Node.RUN:
            print "node {} is busy".format(self.node_dir)
            return -1
        print "restore files in {}".format(self.node_dir)
        for x in calc_files:
<<<<<<< HEAD
            self.run('mv {0}0 {0}'.format(x))
        return 0

=======
            self.run('cp {0}0 {0}'.format(x))
            sleep(0.5)
        return 0

    def make_log(self, out_files, prefix):
        logfile = '{}/{}{}.log'.format(self.node_dir, 
                                       prefix, self.occupied)
        my_out_files = ['{}/{}'.format(self.node_dir, fname) 
                        for fname in out_files]
        f_log = open(logfile,'w')
        f_outs = [open(fname) for fname in my_out_files]
        for f_out, fout_name in zip(f_outs, out_files):
            f_log.write('-'*20+'\n')
            f_log.write('file {} contents\n'.format(fout_name))
            f_log.write('-'*20+'\n')
            for line in f_out:
                f_log.write(line)
        f_log.close()

>>>>>>> 109f3ccc1a8e0b48e3c0512753b53bb04da77d21
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
<<<<<<< HEAD
    res = sc.array(res)
=======
    res = sc.array(res, dtype=sc.float64) # important
>>>>>>> 109f3ccc1a8e0b48e3c0512753b53bb04da77d21
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
<<<<<<< HEAD
    for n, line in enumerate(f):
        if n != nat:
            new_coord += line
            changed += [False]
            continue
        aux = line.split()
        aux[:3] = map(float, aux[:3])
        aux[ix] += eps
        aux[:3] = map(lambda x: "{:.8f}".format(x), aux[:3])
=======
    st = 0
    for n, line in enumerate(f):
        if n != nat:
            new_coord += line
            if not '$end' in line and st<1:
	        changed += [False]
            else:
                st =1    
            continue
        aux = line.split()
        aux[:3] = map(float64, aux[:3])
        aux[ix] += eps
        aux[:3] = map(lambda x: fmt_float.format(x), aux[:3])
>>>>>>> 109f3ccc1a8e0b48e3c0512753b53bb04da77d21
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
<<<<<<< HEAD
          grad=grad_command):
=======
          grad=grad_command): 
>>>>>>> 109f3ccc1a8e0b48e3c0512753b53bb04da77d21
    if node.occupied is None:
        node.occupied = i
    # wait
    while(node.occupied != i):
<<<<<<< HEAD
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
=======
        if node.occupied is None:
	   node.occupied = i
        yield (False, None)
    # wait
    while(node.update_status() != Node.IDLE):
        yield (False, None)
    my_coord   = node.node_dir+'/'+coord
    my_control = node.node_dir+'/'+control
    node.restore_files()
    # wait
    while(node.update_status() != Node.IDLE):
        yield (False, None)
>>>>>>> 109f3ccc1a8e0b48e3c0512753b53bb04da77d21
    changed = update_coord(i, eps, atom_coords, my_coord)
    print "node {} coord {} dscf+".format(node.node_dir, i)
    node.run(srun + dscf)
    while(node.update_status() != Node.IDLE):
<<<<<<< HEAD
        yield None
=======
        yield (False, None)
>>>>>>> 109f3ccc1a8e0b48e3c0512753b53bb04da77d21
    print "node {} coord {} grad+".format(node.node_dir, i)
    node.run(srun + grad)
    # wait
    while(node.update_status() != Node.IDLE):
<<<<<<< HEAD
        yield None
    res = exclude_grads_and_flatten(read_control(my_control, changed)[0],
                                    excl_inds)

    update_coord(i, -2*eps, atom_coords)
=======
        yield (False, None)
    grads, en = read_control(my_control, changed)
    print "{}-th + Energy {}".format(i,en)
    res = exclude_grads_and_flatten(grads, excl_inds)
    node.make_log(out_files=(coord, control), prefix='task_p')
    node.restore_files()
    update_coord(i, -eps, atom_coords, my_coord)
>>>>>>> 109f3ccc1a8e0b48e3c0512753b53bb04da77d21
    print "node {} coord {} dscf-".format(node.node_dir, i)
    node.run(srun + dscf)
    # wait
    while(node.update_status() != Node.IDLE):
<<<<<<< HEAD
        yield None
=======
        yield (False, None)
>>>>>>> 109f3ccc1a8e0b48e3c0512753b53bb04da77d21
    print "node {} coord {} grad-".format(node.node_dir, i)
    node.run(srun + grad)
    # wait
    while(node.update_status() != Node.IDLE):
<<<<<<< HEAD
        yield None
    res = res - exclude_grads_and_flatten(read_control(my_control, changed)[0],
                                          excl_inds)
    res = res/(2*eps)
    node.occupied = None
    node.restore_files()
    yield res
=======
        yield (False, None)
    grads, en = read_control(my_control, changed)
    print "{}-th - Energy {}".format(i, en)
    node.make_log(out_files=(coord, control), prefix='task_m')
    res = res - exclude_grads_and_flatten(grads,
                                          excl_inds)
    res = res/(2*eps)
    node.occupied = None
    yield (True, res)
>>>>>>> 109f3ccc1a8e0b48e3c0512753b53bb04da77d21


labels, atom_coords, excl_inds = make_indexes()
N = len(labels)
print "There are {} coordinates".format(N)
<<<<<<< HEAD
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
=======
nodes = [Node(node_dir_prefix+str(i), calc_files)
         for i in range(number_of_nodes)]
for node in nodes:
    node.save_files()
finished = [False] * N
tasks = [dgdxi(i, nodes[i % number_of_nodes], atom_coords, excl_inds) for i in
         range(N)]
hessian = sc.zeros((N, N),dtype=sc.float64)
while not all(finished):
    for ii, task in enumerate(tasks):
        if finished[ii]:
           continue
        task_finished, res = task.next()
        if task_finished:               # 
            hessian[ii, :] = res[:]
            finished[ii] = True
    sleep(delay)
#    for node in nodes:
#        print "{} status is {}".format(node.node_dir, node.update_status())

print "HESSIAN stored in {}, {}.scipy".format(hessian_file, hessian_file)
>>>>>>> 109f3ccc1a8e0b48e3c0512753b53bb04da77d21
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
<<<<<<< HEAD
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
=======
# write hessian
if N % ncol == 0:
    maxrow = N/ncol -1
    last_col = ncol
else:
    maxrow = N/ncol
    last_col = N % ncol
fmt2 = "{:2d}{:3d} " + fmt_sci*last_col+'\n'
fmt = "{:2d}{:3d} " + fmt_sci*ncol+'\n'
for ii in xrange(N):
    for rown in xrange(maxrow):
        hess_f.write(fmt.format(ii+1, rown+1,
                                *hessian[ii, ncol*rown:ncol*(rown+1)]))
   
    hess_f.write(fmt2.format(ii+1, maxrow+1, *hessian[ii, N-last_col:]))

sc.savetxt(hessian_file+'.scipy', hessian)

mass_vect = sc.array([1/sqrt(at_masses[at]) for at in labels])
mass_mat = sc.outer(mass_vect, mass_vect)/sau2au
omega_mat = hessian*mass_mat
sc.savetxt(omega_file, omega_mat)
energies,vectors = scl.eig(omega_mat)
with open(hessian_file+'.scipy', 'a') as hess_file:
     hess_file.write('N {}\n'.format(N))
     hess_file.write('eps {} delta {}\n'.format(eps, delta))
     hess_file.write('atom labels\n')
     fmt_labels = "{} "*N+'\n'
     hess_file.write("labels xyz ")
     hess_file.write(fmt_labels.format(*labels))
    
print [sqrt(abs(e))*au2cm for e in energies]


>>>>>>> 109f3ccc1a8e0b48e3c0512753b53bb04da77d21
