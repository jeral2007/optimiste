import scipy as sc
from scipy import float64, sqrt
from turboutils import coord_iterator, table
import scipy.linalg as scl
# config
from standard_config import *
hessian_file = 'HESS.DAT'+'.scipy'
omega_file = 'OMEGA.DAT.scipy'
excls = ['zz3', 'zz4']
gaus_file = "GAUSS"
firefly_fname = 'gran.in'

import stat_proc as stp
# HESS.DAT format description:
# N - number of coordinates (number of atoms * 3)
# eps - step value for hessian calculation
# delta - estimated error (Trace of Err**2, Err[i,j] = H[i,j] - H[j, i])
#
#with open(hessian_file+'.scipy', 'a') as hess_file:
#     hess_file.write('N {}\n'.format(N))
#     hess_file.write('eps {} delta {}\n'.format(eps, delta))
#     hess_file.write('atom labels\n')
#     fmt_labels = "{} "*N+'\n'
#     hess_file.write("labels xyz ")
#     hess_file.write(fmt_labels.format(*labels))


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

def read_hess_mat(hessian_file=hessian_file):
    for line in open(hessian_file):
        if 'N' in line:
           break
        yield line

def read_hess_params(hessian_file=hessian_file):
    f_in = open(hessian_file)
    for line_tmp in f_in:
        if 'N' in line_tmp:
           line = line_tmp
           break
    N = int(line_tmp.split()[1])
    line = f_in.next()
    eps, delta = (float64(x) for x in [line.split()[1, 3]])
    line = f_in.next()
    line = f_in.next()
    labels = [lab for lab in line.split() if lab!='labels' and lab!='xyz']

def reduced_mass(i, vectors, labels, at_masses=at_masses):
    av =  sum(1e0/at_masses[lab]*vectors[j, i]**2 for (j,lab) in enumerate(labels))
    res = 1e0/sum(1e0/at_masses[lab]**2*vectors[j, i]**2/av for (j,lab) in enumerate(labels))
    return res   

def make_firefly_input(hess, ncol=5, firefly_fname=firefly_fname, 
                             coord=coord, coord_pat=coord_pat, excls = excls):
    import fortranformat as ff
    gess_fmt_fortran="(I2,I3,1P,{}E15.8)".format(ncol)
    gess_row = ff.FortranRecordWriter(gess_fmt_fortran)
    f_c = coord_iterator(coord)
    f_pat = coord_iterator(coord_pat)
    f_firefly = open(firefly_fname,'w')
    f_firefly.write("""\
 $contrl scftyp=UHF runtyp=HESSIAN mplevl=0 units=BOHR icharg=0 mult=2 $end
 $system timlim=30 memory=500000 $end
 $basis gbasis=DZV  $end
 $guess guess=huckel $end
 $FORCE RDHESS = .TRUE. PURIFY=.true. temp(1)=0.001, 298.15, $end
""")
    #masses and coords
    coord_section="""\
 $data
oh pbe0
C1
"""
    h_ind=5
    at_names = []
    for x, xpat in zip(f_c, f_pat):
        if xpat[3] in excls:
            continue
        at_names += [xpat[3]]
        coord_section +="h {:1d}    {:.7f}   {:.7f}    {:.7f}\n".format(h_ind, *xpat[:3])
        h_ind=6
    coord_section += '$end\n'
    print at_names
    fmt_mass =  " {:.5f},"*len(at_names)
    mass_section =("""\
 $MASS
   AMASS(1)= FIXME  """+fmt_mass+"""
 $END
""").format(*[at_masses[at] for at in at_names])
    f_firefly.write(mass_section)
    f_firefly.write(coord_section)
    f_firefly.write("""\
$HESS
 ! Total energy      =      -339.0052743001
""")
    N = hess.shape[0]
    if N % ncol == 0:
        maxrow = N/ncol -1
        last_col = ncol
    else:
        maxrow = N/ncol
        last_col = N % ncol
    for trow in range(N):
        for srow in range(maxrow):
            line = gess_row.write([trow+1,srow+1]+list(hess[trow,srow*ncol:(srow+1)*ncol]))
            f_firefly.write(line+'\n')
        line = gess_row.write([trow+1,srow+1]+ list(hess[trow, N-last_col:]))

        f_firefly.write(line+'\n')
    f_firefly.write('$end\n')
    f_firefly.close()


labels, atc_i, excl_inds = make_indexes()   
hessian = sc.loadtxt(read_hess_mat())
mass_vect = sc.array([1/sqrt(at_masses[at]) for at in labels])
mass_mat = sc.outer(mass_vect, mass_vect)/sau2au
omega_mat = hessian*mass_mat
sc.savetxt(omega_file, omega_mat)
eigs, vectors = scl.eig(omega_mat)
print hessian
print '---------------'
freqs_and_masses = [(sqrt(sc.real(e))*au2cm, reduced_mass(i, vectors, labels)) for (i,e) in enumerate(eigs)]
freqs_and_masses.sort(key= lambda t: t[0])
for n,(f,w) in enumerate(freqs_and_masses):
    print 
print freqs_and_masses

freqs = sc.array([f for (f, m) in freqs_and_masses])
masses = sc.array([m for (f, m) in freqs_and_masses])
ordering = [ w  for (w,m) in freqs_and_masses]
print 'analysis'
print '---------------'
f_m_clusters = stp.per_clusters(ordering, 35)
gaussian_freqs = []
f_m_clusters.sort(key=lambda t: t[0])
clust_an = []
for n, c_m in enumerate(f_m_clusters):
    c = c_m[1]
    ws = freqs[c] 
    ms = masses[c]
    gaussian_freqs +=[sc.mean(ws)]
    avg_w =sc.mean(ws)
    err_w = 0.5*(max(ws)-min(ws))
    avg_m = sc.mean(ms)
    err_m = 0.5*(max(ms)-min(ms))
    clust_an += ["{}. {} modes from {}-th mode:w = {:.1f} += {:.1f} cm^-1,\
m = {:.1f} +- {:.1f} s.a.u.".format(n+1, len(c),  c[0]+1, avg_w, err_w, avg_m, err_m)]   

print table(clust_an, 2)

print "firefly file: {}".format(firefly_fname)
make_firefly_input(hessian)

print "gaussian_an"
f_gaus = open(gaus_file,'w')   

if len(f_m_clusters)*2 % 3 == 0:
   n_blocks = len(f_m_clusters)*2 / 3
   last_block_len = 3
else:
   last_block_len = len(f_m_clusters)*2 % 3


to_print = []
for n, c_m in enumerate(f_m_clusters):
    fi, li = c_m[1][0], c_m[1][-1]
    to_print +=[(fi, li, gaussian_freqs[n])]

kk=1
for ii1, ii2, omega in to_print:
   f_gaus.write("                   {}           \n".format(kk))
   f_gaus.write("                   A            \n")
   f_gaus.write("Frequencies -- {}\n ".format(omega))
   f_gaus.write("Atom  AN      X      Y      Z \n")
   for jj in range(len(vectors[:,ii1])/3):
       f_gaus.write(("{:3d} {:3d} "+"{:.2f}  "*3+'\n').format(jj+1, at_numbers[labels[jj]],*vectors[jj:jj+3,ii1]))
   kk+=1
   f_gaus.write("                   {}           \n".format(kk))
   f_gaus.write("                   A          \n")
   f_gaus.write("Frequencies -- {}\n ".format(omega))
   f_gaus.write("Atom  AN      X      Y      Z \n")
   for jj in range(len(vectors[:,ii2])/3):
       f_gaus.write(("{:3d} {:3d} "+"{:.2f}  "*3+'\n').format(jj+1, at_numbers[labels[jj]],*vectors[jj:jj+3,ii2]))
   kk+=1
#                   1                      2                      3
"                   A                      A                      A"
#Frequencies --    57.7280                75.0007                75.0254
# Atom  AN      X      Y      Z        X      Y      Z        X      Y      Z
#    1  43     0.00   0.00   0.00     0.01  -0.02   0.00    -0.02  -0.01   0.00
#    2   6     0.01  -0.01   0.20     0.05   0.00  -0.04     0.14   0.15   0.01
#    3   6    -0.01   0.01   0.20     0.05   0.00   0.04     0.14   0.15  -0.01
#    4   6     0.01   0.01  -0.20    -0.15   0.14  -0.01     0.00  -0.05  -0.04
#    5   6    -0.01  -0.01  -0.20    -0.15   0.14   0.01     0.00  -0.05   0.04
#    6   6     0.00   0.00   0.00     0.07  -0.10   0.00    -0.09  -0.07   0.00
#    7   8    -0.02   0.02   0.46     0.09   0.04   0.11     0.39   0.40  -0.02
#    8   8     0.00   0.00   0.00     0.16  -0.23   0.00    -0.23  -0.16   0.00
#    9   8     0.02   0.03  -0.46    -0.40   0.39  -0.01     0.05  -0.09  -0.11
#   10   8     0.02  -0.03   0.45     0.09   0.05  -0.11     0.39   0.40   0.02
#   11   8    -0.03  -0.02  -0.46    -0.40   0.39   0.02     0.05  -0.09   0.11
#   12  17     0.00   0.00   0.00     0.22  -0.31   0.00    -0.31  -0.22   0.00


