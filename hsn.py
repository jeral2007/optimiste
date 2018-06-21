#!/usr/bin/env python
import sys
import scipy as sc
from time import sleep
from read_control import read_control
from para_run import *
# config
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




def run_dscf(node, result, srun=srun_command, dscf=dscf_command):
    print "dscf run in {}".format(node.node_dir)
    node.run("{} {}".format(srun, dscf))


def run_grad(node, result, srun=srun_command, grad=grad_command):
    print "grad run in {}".format(node.node_dir)
    node.run("{} {}".format(srun, grad))


def make_update_coord(i, eps, atom_coords, coord=coord):
    nat, ix = atom_coords[i]

    def update_coord(node, result):
        node.save_files()
        print "updating {}-th coordinate by {}".format(i, eps)
        my_coord = "{}/{}".format(node.node_dir, coord)
        f = open(my_coord)
        header = f.next()
        if '$coord' not in header:
            raise ValueError('Invalid coord file {}'.format(my_coord))
        new_coord = header
        for n, line in enumerate(f):
            if n != nat:
                new_coord += line
                continue
            aux = line.split()
            aux[:3] = map(float, aux[:3])
            aux[ix] += eps
            aux[:3] = map(lambda x: "{:.6f}".format(x), aux[:3])
            new_coord += " ".join(aux)+'\n'
        f.close()
        f = open(my_coord, 'w')
        f.write(new_coord)
        f.close()

    return (update_coord, True)



def get_grad(node, result, control=control, coord=coord):
    my_control = "{}/{}".format(node.node_dir, control)
    f = open(coord)
    f.next()
    atoms = []
    for line in f:
        if '$end' in line:
            break
        atoms += [line.split()[3]]
    f.close()
    grads, en = read_control(my_control, atoms)
    result[0] = grads
    node.restore_files()


def exclude_grads_and_flatten(grads, excl_inds):
    res = [g for (n, g) in enumerate(grads) if n not in excl_inds]
    res = sc.array(res)
    return sc.reshape(res, 3*len(res))



labels, atom_coords, excl_inds = make_indexes()
N = len(labels)
print "There are {} coordinates".format(N)
print atom_coords
# make nodes and tasks
nodes = [Node(node_dir_prefix+str(i), calc_files, strict=strict)
         for i in range(number_of_nodes)]

tasks = []
for i in xrange(len(labels)):
    nd_ind = 2*i % number_of_nodes
    tasks += [Task(nodes[nd_ind], (make_update_coord(i, eps, atom_coords),
                                  (run_dscf, False),
                                  (run_grad, False),
                                  (get_grad, True)))]

    nd_ind = (2*i+1) % number_of_nodes
    tasks += [Task(nodes[nd_ind], (make_update_coord(i, -eps, atom_coords),
                                  (run_dscf, False),
                                  (run_grad, False),
                                  (get_grad, True)))]


# run tasks
finished = False
while not finished:
    for task in tasks:
        task.run()
    finished = all(t.finished for t in tasks)
    error = any(node.status() == 'ERROR' for node in nodes)
    for iin, node in enumerate(nodes):
        print "node {} status is {}".format(iin, node.status())
    if error and strict:
        raise ValueError("some node return error status")
    sleep(delay)


hessian = sc.zeros((N,N))
delta = 0e0
grads = sc.zeros(N)
for jj, task in enumerate(tasks):
    ii = jj / 2
    tp = task
    tm = tasks[ii+1]
    dg = (exclude_grads_and_flatten(tp.result[0], excl_inds) -
          exclude_grads_and_flatten(tm.result[0], excl_inds))/(2*eps)

    grads[:] = (exclude_grads_and_flatten(tp.result[0], excl_inds) +
                exclude_grads_and_flatten(tm.result[0], excl_inds))/2e0
    hessian[ii, :] = dg[:]

print "-"*20
print "grads"
print "-"*20
grads = sc.reshape(grads, (N/3, 3))
fmt = "{}. \t {} \t "+" \t ".join([fmt_float]*3)
for ii in xrange(N/3):
    print fmt.format(ii, labels[ii],*grads[ii,:])
print "-"*20
print "HESSIAN stored in {}".format(hessian_file)
print "-"*20
hess_f = open(hessian_file, 'w')

fmt = (fmt_float+" \t ") * N
for ii in xrange(N):
    print fmt.format(*hessian[ii,:])
    for jj in xrange(ii):
        delta = delta + (hessian[ii, jj] - hessian[jj, ii])**2

print "delta = {}".format(delta)
hessian = (hessian +sc.transpose(hessian))*0.5
fmt = "{:2d}{:3d} "+fmt_sci*nrow+'\n'
for ii in xrange(N):
    for coln in xrange(N/nrow-1):
        hess_f.write(fmt.format(ii+1, coln+1,*hessian[ii,nrow*coln:nrow*(coln+1)]))
    if coln is None:
       coln = 0
    else:
       coln +=1
    hess_f.write(fmt.format(ii+1,coln+1,*hessian[ii,N-nrow:]))

