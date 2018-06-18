import sys
from read_control import take_grads
import scipy as sc
coord_pat = 'coord_pat'
control = 'control'
turbo_command = '(/home/demidov/turbo/bin/amd64/part_dscf>out) && (/home/demidov/turbo/bin/amd64/part_grad >grad.out)'
relax_command = '/home/demidov/turbo/bin/amd64/part_relax>relax.out'
coord = 'coord'
coord_old='coord_old'
excls = ['y', 'zz3', 'zz4']
maxiter=40
damping = 0.5
def del_num(s):
    return ''.join(c for c in s if c not in '0123456789')

def save_coords(coord=coord, coord_old=coord_old):
    from subprocess import call
    call ("cp {} {}".format(coord, coord_old), shell=True)
    
def write_control(control=control, coord_pat=coord_pat, excls=excls):

    f_control = open(control)
    f_pat = open(coord_pat)
    control_cont = ""
    header = f_pat.next()
    if '$coord' not in header:
       raise ValueError('Invalid coord_pat file')
    exclude_q = []
    
    for line in f_pat:
        if '$end' in line:
           break
        aux = line.split()
        exclude_q += [aux[3] in excls]
    f_pat.close()
    print len(exclude_q)
    for line in f_control:
        control_cont+=line
        if "$grad          cartesian gradients" in line:
            break
    first_step = True
    while(True):
        line = f_control.next()
        control_cont+=line
        if not first_step and '$last' in line or '$actual' in line:
           break 
        if '$end' in line:
           break 
        if 'maximum norm' in line:
           break 
        #coords
        for exc, line in zip(exclude_q, f_control):
           aux = line.split()
           if len(aux)<3: 
              raise ValueError("Invalid control file -- {} instead of atom coords".format(line))
           control_cont += line
        #gradients
        for exc, line in zip(exclude_q, f_control):
           aux = line.split()
           if len(aux) !=3: 
              raise ValueError("Invalid control file -- {} instead of atom gradient".format(line))
           if exc:
              control_cont += "     0.000000000  0.000000000 0.000000000 \n"
           else:
              control_cont+=line
        first_step=False
    for line in f_control:
        control_cont+=line
    f_control.close()
    f_control = open('control','w')
    f_control.write(control_cont)
    f_control.close()

def one_iter(turbo_command=turbo_command, relax_command=relax_command):
    import subprocess as sp
    print "dscf + grad run"
    sp.call(turbo_command, shell=True)
    print "relax run"
    save_coords()
    write_control()
    sp.call(relax_command, shell=True)

def restore_fixed(coord_pat=coord_pat, coord=coord, coord_old=coord_old,  excls = excls, damping=damping):
    f_pat = open(coord_pat, 'r')
    f_coord = open(coord, 'r')
    f_old = open(coord_old, 'r')
 
    header_pat = f_pat.next()
    header_coord = f_coord.next()
    header_old = f_old.next()
    
    if '$coord' not in header_pat or '$coord' not in header_coord:
       raise  ValueError('invalid coord or coord_pat file')
    if '$coord' not in header_old:
       raise ValueError('invalid coord_old file')

    res_coord = []
    delta = []
    for line_pat, line_coord, line_old in zip(f_pat, f_coord, f_old):
        end_of_coord = ['$end' in line for line in [line_coord, line_pat]]
        if any(end_of_coord):
           if not all(end_of_coord):
              raise ValueError('Inconsistent coord and coord_pat files (different number of atoms)')
           else:
              break
        aux_pat = line_pat.split()
        aux_coord = line_coord.split()
        if len(aux_coord)<4:
           aux_coord +=[del_num(aux_pat[3])]

        if aux_pat[3] in excls:
           aux_coord[:3] =aux_pat[:3]
           aux_coord[3] = del_num(aux_pat[3])
        else:
           aux_old = line_old.split()
           new_coords = [float(x0) + (float(x)-float(x0))*damping for x0, x in zip(aux_old[:3], aux_coord[:3])]
           aux_coord[:3] = map(str, new_coords)
           delta += [(sum((float(x)-float(x0))**2 for x, x0 in zip(aux_pat[:3], aux_coord[:3]))**0.5, 
                     aux_pat[3])]
        res_coord += [' '.join(aux_coord)+'\n']
    
    print " writing to coord"
    f_coord.close()
    f_pat.close()
    f_coord = open(coord, 'w')
    f_coord.write('$coord\n')
    for line in res_coord:
       f_coord.write(line)
    f_coord.write('$end\n')
    f_coord.close()
    return delta


save_coords()
restore_fixed()
for i in xrange(maxiter):
    one_iter()
    print restore_fixed()
    take_grads(coord_pat, control, excls)
