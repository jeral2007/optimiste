import sys
import scipy as sc
import scipy.linalg as scl
from scipy import float64, sqrt # sort of stupid
from time import sleep
from read_control import Grads
from standard_config import *
import subprocess as sp
import sys
dist = 1.5e2
excls = 'zz3 zz4 '.split()
embedding_pat ='embedding_pat'
calc_files = [coord, realmos, imagmos, basis, embedding]
dscf_command = 'part_dscf > out'
grad_command = 'part_grad >grad.out'
dipole_txt = 'dipole.dat'
def make_embed(params, embedding=embedding,embedding_pat=embedding_pat):
    fin  = open(embedding_pat)
    fout = open(embedding,'w')
    for line in fin:
        fout.write(line.format(*params))

def eval_dipole(charge, dist=dist):
    print "charge = {:.4f} a. u.".format(charge)
    field = -2*charge/dist**2
    Pm, dPm = [], []
    d2Pm = 0
    print "|E| = {}".format(field)
    backup()
    charges = [0e0]*6
    make_embed(charges)
    sp.call(dscf_command, shell=True)
    sp.call(grad_command, shell=True)
    grd_obj = Grads(excls, coord_pat, control)
    e0 = grd_obj.en
    d2Pm = -6*grd_obj.flatten_grads()
    print "E0 = {} a. u.".format(e0)
    restore()
    for ii in xrange(0, len(charges), 2):
        print "{}-th cartesian coordinate".format(ii/2)
        charges = [0e0]*6
        charges[ii] = charge
        charges[ii+1] = -charge
        make_embed(charges)
        sp.call(dscf_command, shell=True)
        sp.call(grad_command, shell=True)
        grd_obj = Grads(excls, coord_pat, control)
        ep = grd_obj.en
        gp = grd_obj.flatten_grads()
        print "E+ = {} a.u.".format(ep)
        backup(backup_suf='{}pbackup'.format(ii))
        restore()
        charges = [0e0]*6
        charges[ii] = -charge
        charges[ii+1] = charge
        make_embed(charges)
        sp.call(dscf_command, shell=True)
        sp.call(grad_command, shell=True)
        grd_obj = Grads(excls, coord_pat, control)
        em = grd_obj.en
        gm = grd_obj.flatten_grads()
        print "E- = {} a. u.".format(em)
        print "delta_E = {} a. u.".format(ep-em)
        Pm += [(ep-em)/(2*field)]
        dPm += [(gp-gm)/(2*field)]
        d2Pm +=(gp+gm)
        print "P_{} = {} a. u.".format(ii/2, Pm[-1])
        backup(backup_suf='{}mbackup'.format(ii))
        restore()
    return Pm,dPm,d2Pm/(field**2)

def backup(calc_files=calc_files, backup_suf='.backup0'):
    for filename in calc_files:
        sp.call("cp -r {0} {0}{1}".format(filename, backup_suf),
                shell=True)

def restore(calc_files=calc_files, backup_suf='.backup0'):
    for filename in calc_files:
        sp.call("cp -r {0}{1} {0}".format(filename, backup_suf),
                shell=True)

charge = float64(sys.argv[1])
atom_type = sys.argv[2]
grd_obj = Grads(excls, coord_pat, control)

de =  grd_obj.flatten_grads()
print "E0 = {}".format(grd_obj.en)
px,dpx,d2px =  eval_dipole(charge)
f = open(dipole_txt, 'w')
f.write('#p\n')
sc.savetxt(f, px)
for ii,dpi in enumerate(dpx):
    f.write('#dp {}\n'.format(ii))
    sc.savetxt(f, dpi)
f.write('#d2p\n')
sc.savetxt(f, d2px)

