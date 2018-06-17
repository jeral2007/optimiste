import sys
from read_control import take_grads
from parse_coord import make_coord_bas
import scipy as sc
from scipy.optimize import minimize

#config
turbo_command = '(/home/demidov/turbo/bin/amd64/part_dscf>out) && (/home/demidov/turbo/bin/amd64/part_grad >grad.out)'
initial_args = sc.array(map(float,sys.argv[1:]))
constr = [(-5e-1, 2.5e0)] * len(initial_args)
excls = ['ne3', 'ne4', 'ne5', 'q6', 'q7', 'q8','q9','o2']
cpat_filename = 'coord_pat'
contr_filename = 'control'
basis_pat_pat = 'basis_pat_pat'
bpat = 'basis'

def make_basis_pat(args, bpatpat=basis_pat_pat, bpat=bpat):
    f = open(bpatpat)
    outf = open(bpat,'w')
    state = 0
    for line in f:
       if '!ecp_begin' in line:
          assert(state==0)
          state=1
          n = int(line.split()[1])-1
          continue
       elif '!ecp_end' in line:
          assert(state==1)
          state=0
          continue
       elif state==0:
          outf.write(line)
          continue
       elif state==1:
          if any(c in line for c in 'spdfgh#'):
             outf.write(line)
             continue
          aux = line.split()
          x = args[n]
          aux[2] =str((1+x)**(-2)*float(aux[2]))
          p = int(aux[1])-2
#         aux[0] =str(float(aux[0])/(1+x)**p)
          print " ".join(aux+ [str(p)])
          outf.write("    ".join(aux)+'\n')

def target(args, turbo_command=turbo_command, excls=excls, 
            cpat_filename=cpat_filename, contr_filename=contr_filename,
            bpat_filename=bpat):
      import subprocess as sp
      print args
      make_basis_pat(args)
#      make_coord_bas(cpat_filename, bpat_filename)
      sp.call(turbo_command, shell=True)
      res = take_grads(cpat_filename, contr_filename, excls)
      return res


minimize(target, initial_args, bounds =constr, method='L-BFGS-B', options={'eps':5e-3, 
        'gtol':1e-3, 'maxfun':400, 'maxiter':20})


