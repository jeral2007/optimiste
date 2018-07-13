#!/usr/bin/env python
import sys
import scipy as sc
from scipy.optimize import minimize
from turboutils import table
#config
filename = "fersmite"
bound = (-0.5e0, 2.5e0)
alpha_step = 3e-4
crystall_command = "mpirun -np {task_num} ./Pssc 2>{out_filename}"
#en config
pat_filename = filename +'.d12_pat'
d12_filename = 'INPUT'
f9_filename = filename +'.f9'
out_filename = filename+'.out'
task_num = 144
initial_args =sc.array(map(float, sys.argv[1:]))
crystall_command = crystall_command.format(task_num=task_num,
                                          out_filename=out_filename)
bounds = [bound]*len(initial_args)
energies = []
maxcol = len(initial_args) + 1
    
def read_energy(out_filename=out_filename):
    f_in = open(out_filename)
    state=0
    n_iter = 0
    for line in f_in:
        if 'CRYSTAL' in line:
          state = 1
          continue
        elif 'ETOT' in line and state==1:
          print line
          en = float(line.split()[3])
          n_iter+=1

    if state == 0:
       raise ValueError('{} is not crystal out file'.format(out_filename))
    print "exit after {} iterations".format(n_iter)
    print "ETOT = {}".format(en)
    return en

def read_grad(out_filename=out_filename):
    f_in = open(out_filename)
    state=0
    for line in f_in:
        if 'CRYSTAL09' in line:
          state = 1
          continue
        elif 'RMS GRADIENT' in line and state==1:
          print line
          return float(line.split()[2]) 

    if state == 0:
       raise ValueError('{} is not crystall09 out file'.format(out_filename))
    else:
       raise ValueError('{} is invalid crystall09 out file'.format(out_filename))
      
def make_d12_filename(args, pat_filename=pat_filename, d12_filename=d12_filename):
    f_in = open(pat_filename)
    f_out = open(d12_filename, 'w')
    state = 0
    r_ls = [0e0]*len(args)
    for line in f_in:
       if '!begin_basopt' in line:
          if state == 0:
             n_arg = int(line.split()[1]) - 1
             state = 1
             continue
          else:
             raise ValueError("{}\n start begin_basopt group before closing previous one".format(line))
       if "!end" in line:
          if state == 0:
             raise ValueError("{}\n closing begin_basopt group before opening".format(line))
          else:
             state = 0
             continue
       if state==1:
          aux = line.split()
          alpha, c  = float(aux[0]), float(aux[1]) 
          k = args[n_arg]
          f_out.write("   {:.16f}     {:.16f}  \n".format(alpha/(1+k)**2, c))
          r_ls[n_arg] = 1./alpha**0.5*(1+k)
          continue

       f_out.write(line)
    return r_ls
    
def target(args, crystall_command=crystall_command, energies=energies):
    import subprocess as sp 
    print args
    print "r_ls (a.u.)"
    try:
        r_ls = make_d12_filename(args)
        print r_ls
        sp.call(crystall_command, shell=True)
        sp.call('cp fort.9 fort.20', shell=True)
        energies += [read_energy()]
    except Exception:
        save_res_and_report()
        raise
    return energies[-1]

def save_res_and_report(energies=energies, maxcol=maxcol):
    import subprocess as sp
    from datetime import datetime
    suf = datetime.now().strftime('%d-%m-%y-%H-%M')
    sp.call('cp fort.9 wf_final'+suf, shell=True)
    sp.call('rm *.pe*', shell=True)
    sp.call('cp INPUT INPUT_FINAL'+suf, shell=True)
    if len(energies) ==0:
       return
    print "-"*20
    print "finished at "+suf
    print "there are {} iterations".format(len(energies))
    print "total energies on each iteration:"
    print table(energies, maxcol, field_fmt="{:.4f} ")
    print "-"*20

minimize(target, initial_args, bounds=bounds, method='L-BFGS-B', options={'eps':alpha_step,
        'gtol':1e-4, 'maxiter':20})
save_res_and_report()
