import sys
from read_control import take_grads
import scipy as sc
from scipy.optimize import minimize

#config

#turbo_command = '(dscf>/dev/null) && (grad >/dev/null)'
c_args = 'q12 q3 q45'.split()

bounds_y = [(2.1, 3.2)] 
bounds_o = [(-1.5, -0.1)]

coord_pat = 'coord_pat'
basis_pat = 'basis_pat'
embedding_pat = 'embedding_pat'
embedding = 'embedding'
control = 'control'
turbo_command = '(/home/demidov/turbo/bin/amd64/part_dscf>out) && (/home/demidov/turbo/bin/amd64/part_grad >grad.out)'

excls = ['ne3', 'ne4', 'ne5', 'q6', 'q7', 'q8','q9', 'o2']
initial_vals = sc.array(map(float,sys.argv[1:]))
bounds = bounds_y*2 + bounds_o

def make_embedding(env, embedding_pat=embedding_pat, embedding=embedding):
    f_in = open(embedding_pat)
    f_out = open(embedding, 'w')
    for line in f_in:
        if '!' in line:
           aux = line.split('!')
           f_out.write("{} {}\n".format(aux[0], eval(aux[1], env)))
        else:
           f_out.write(line)

def target(args, coord_pat=coord_pat, basis_pat=basis_pat, c_args=c_args,
           control=control,excls=excls):
    import subprocess as sp
    print args
    c_env = {k:v for k,v in zip(c_args, args)}
    make_embedding(c_env)
    sp.call(turbo_command, shell=True)
    res = take_grads(coord_pat, control, excls)
    return res


minimize(target, initial_vals, bounds=bounds, method='L-BFGS-B', options={'eps':5e-2, 
        'gtol':1e-4, 'maxfun':400, 'maxiter':20})

