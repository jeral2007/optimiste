#!/usr/bin/env python
import sys
from read_control import take_grads
from parse_coord import make_coord_bas
import scipy as sc
from scipy.optimize import minimize
#config
c_args = {}
turbo_command = '(dscf>/dev/null) && (grad >/dev/null)'
c_args = 'q1 q2'.split()
initial_vals = sc.array(map(float,sys.argv[1:]))
bounds = [(0.2, 2.9)]*len(initial_vals)
coord_pat_pat = 'coord_pat_pat'
coord_pat = 'coord_pat'
basis_pat = 'basis_pat'
control = 'control'
excls = ['ne1', 'ne2']
#end config
def make_coord_pat(env, coord_pat_pat=coord_pat_pat, coord_pat=coord_pat):
    f = open(coord_pat_pat)
    f_out = open(coord_pat,'w')
    header = f.next()
    assert('$coord' in header)
    f_out.write(header)
    state = 1
    for line in f:
      if '!' in line:
         aux = line.split('!')
         line2 =  "{} {}". format(aux[0].strip(), eval(aux[1], env))
      else:
         line2 =  line.strip()
      if '$end' in line2:
         state = 2
      if state == 2:
         f_out.write(line2+'\n')
         continue

      if 'sc' == line2[:2]:
          aux = line2[2:].split()
          aux[:3] = [env['sc']*float(x) for x in aux[:3]]
          f_out.write(" ".join(str(x) for x  in aux)+'\n')
      else:
          f_out.write(line2+'\n')

def target(args, coord_pat=coord_pat, basis_pat=basis_pat, c_args=c_args,
           control=control,excls=excls):
    import subprocess as sp
    print args
    c_env = {k:v for k,v in zip(c_args, args)}
    make_coord_pat(c_env)
    make_coord_bas(coord_pat, basis_pat)
    sp.call(turbo_command, shell=True)
    res = take_grads(coord_pat, control, excls)
    return res


minimize(target, initial_vals, bounds=bounds, method='L-BFGS-B', options={'eps':5e-2, 
        'gtol':1e-4, 'maxfun':400, 'maxiter':20})
