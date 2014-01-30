#!/usr/bin/env python
import sys
import dirac_bindings as db
import argparse
import quadratic as qu
import os.path
   
parser  = argparse.ArgumentParser()
parser.add_argument('inp',type = str, help = "dirac input  filename (Atom.inp for example")
parser.add_argument('mol',type = str, help = "dirac xyz or mol filename(that to be constructed)")
parser.add_argument('opt_mol',type = str, help = "optimiste mol filename(from that mol constructed)")
args = parser.parse_args()
dirac_output="{0}_{1}.out".format(os.path.splitext(args.inp)[0],os.path.splitext(args.mol)[0])
energy = db.MakeGetEnergyFunc(db.dirac_en)
rp = db.MakeRunProgram('pam-dirac --inp={}'.format(args.inp),db.transform_input_eval)

def optimum(dist):
   rp(args.opt_mol,args.mol,dist) #run program 
   en = energy(dirac_output)
   sys.stderr.write("dist={0}, en = {1}\n".format(dist,en))
   return en #get_energy

x = qu.quadr_optimize(optimum,[3.6,3.9,4.1],eps=1e-2,maxiter = 5)
print ("eq. distance is {}".format(x))

