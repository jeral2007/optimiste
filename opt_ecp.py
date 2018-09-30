#!/usr/bin/env python
import sys
import scipy as sc
from scipy.optimize import minimize
from turboutils import table, read_coord
from read_control import read_control
#config
filename = "basis"
coord_filename= 'coord'
bound = (-0.5e0, 2.5e0)
alpha_step = 2e-3
calc_command = "dscf >{out_filename} && dscf >>{out_filename} && grad >{out_filename}.grad"
#en config
pat_filename = filename +'_pat'
basis_filename = filename
control_filename='control'
out_filename = 'out'
task_num = 144
initial_args =sc.array(map(float, sys.argv[1:]))
calc_command = calc_command.format(task_num=task_num,
                                          out_filename=out_filename)
bounds = [bound]*len(initial_args)
energies = []
maxcol = len(initial_args) + 1
atoms = read_coord(coord_filename)

def read_energy(control_filename=control_filename, atoms=atoms):
    return read_control(control_filename, atoms)[1]

def make_basis_filename(args, pat_filename=pat_filename, basis_filename=basis_filename):
    f_in = open(pat_filename)
    f_out = open(basis_filename, 'w')
    state = 0
    r_ls = [0e0]*len(args)
    for line in f_in:
       if '!begin_ecpopt' in line:
          if state == 0:
             n_arg = int(line.split()[1]) - 1
             state = 1
             continue
          else:
             raise ValueError("{}\n start begin_ecpopt group before closing previous one".format(line))
       if "!end" in line:
          if state == 0:
             raise ValueError("{}\n closing begin_basopt group before opening".format(line))
          else:
             state = 0
             continue
       if state==1:
          aux = line.split()
          c, nr_ecp, alpha  = float(aux[0]), float(aux[2])
          k = args[n_arg]
          f_out.write("   {:.16f}  {}   {:.16f}  \n".format(c, nr_ecp, alpha/(1+k)**2))
          r_ls[n_arg] = 1./alpha**0.5*(1+k)
          continue

       f_out.write(line)
    return r_ls

def target(args, calc_command=calc_command, energies=energies):
    import subprocess as sp
    print args
    print "r_ls (a.u.)"
    try:
        r_ls = make_basis_filename(args)
        print r_ls
        sp.call(calc_command, shell=True)
        grads += [read_grad()]
    except Exception:
        save_res_and_report()
        raise
    return energies[-1]

def save_res_and_report(grads=grads, maxcol=maxcol):
    import subprocess as sp
    from datetime import datetime
    suf = datetime.now().strftime('%d-%m-%y-%H-%M')
    sp.call('cp basis basis_FINAL'+suf, shell=True)
    if len(grads) ==0:
       return
    print "-"*20
    print "finished at "+suf
    print "there are {} iterations".format(len(energies))
    print "total squared grads  on each iteration:"
    print table(energies, maxcol, field_fmt="{:.4f} ")
    print "-"*20

minimize(target, initial_args, bounds=bounds, method='L-BFGS-B', options={'eps':alpha_step,
        'gtol':1e-8, 'maxiter':20})
save_res_and_report()
