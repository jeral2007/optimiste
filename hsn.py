#!/usr/bin/env python
import sys
import scipy as sc
import subprocess as sp
from read_control import read_control
# config
coord = 'coord'
control = 'control'
realmos = 'realmos'
imagmos = 'imagmos'
mos = 'mos'
coord_pat = 'coord_pat'
excls = 'zz3 zz4'.split()
node_dir_prefix = 'node'
srun_command = 'srun -N1 -n24'
dscf_command = 'dscf>out'
grad_command = 'grad>grad.out'
eps = 5e-3
#
number_of_nodes = int(sys.argv[1])

class Node(object):
    def __init__(self, node_dir, coord=coord,
                 control=control, realmos=realmos,
                 imagmos=imagmos, mos=mos):
        sp.call('mkdir {}'.format(node_dir), shell=True)
        self.calc_files = [node_dir+'/'+x for x in
                           [coord, control, realmos, imagmos, mos]]
        for nf, f in zip(self.calc_files,[coord, control, realmos,
                                          imagmos, mos]):
            sp.call('cp {} {}'.format(f, nf), shell=True)

        self.node_dir = node_dir
        self.task = None

    def save_files(self):
        for filename in self.calc_files:
            sp.call('cp {} {}0'.format(filename, filename))

    def restore_files(self):
        for filename in self.calc_files:
            sp.call('cp {}0 {}'.format(filename, filename))

    def free(self):
        sp.call('rm -rf {}'.format(self.node_dir))

    def run(self, command):
        if self.task is None or self.task.poll() == 0:
            self.task = sp.Popen(command, shell=True, cwd=self.node_dir)

    def status(self):
        if self.task is None:
            return 'IDLE'
        elif self.task.poll() is None:
            return 'rUNNING'
        elif self.task.poll() == 0:
            return 'FINISHED'
        else:
            return 'ERROR'


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


labels, atom_coords, excl_inds = make_indexes()


def update_coord(i, eps, coord, atom_coords=atom_coords):
    f = open(coord)
    nat, ix = atom_coords[i]
    header = f.next()
    if '$coord' not in header:
        raise ValueError('Invalid coord file')
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
    f = open(coord, 'w')
    f.write(new_coord)
    f.close()

def take_grads(control, coord=coord):
    f = open(coord)
    f.next()
    atoms = []
    for line in f:
        if '$end' in line:
            break
        atoms += [line.split()[3]]
    f.close()
    grads, en = read_control(control, atoms)
    return sc.array(grads)
def make_process_index(i, result, node, srun=srun_command,
                       dscf=dscf_command,
                       grad=grad_command, control=control, coord=coord,
                       eps=eps):
    state = 0
    result = []
    def process_index():
        if node.status() != 'IDLE' and node.status() != 'FINISHED':
            return

        if state == 0:
            print "first dscf run on index {}".format(i)
            node.save_files()
            update_coord(i, eps, coord="{}/{}".format(node.node_dir, coord))
            node.run('{} {}'.format(srun, dscf))
            state = 1
            return
        elif state == 1:
            print "first grad run on index {}".format(i)
            node.run('{} {}'.format(srun, grad))
            state +=1
            return
        elif state == 2:
            result = take_grads(node.node_dir+'/'+control)
            node.restore_files()
            update_coord(i, -eps, coord="{}/{}".format(node.node_dir, coord))
            print "second dscf run on index {}".format(i)
            node.run('{} {}'.format(srun, dscf))
            state += 1
            return
        elif state == 3:
            print "second grad run on index {}".format(i)
            node.run('{} {}'.format(srun, grad))
            state += 1
            return
        elif state == 4:
            result = result - take_grads(node.node_dir+'/'+control)
            result = result*0.5e0
            node.restore_files()
            state+=1
            return result
        else:
            return result
    return process_index


nodes = [Node(node_dir_prefix+str(i)) for i in range(number_of_nodes)]

