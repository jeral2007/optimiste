#!/usr/bin/env python
import sys
sys.path+=['../'] 
import dirac_bindings as db
energy = db.MakeGetEnergyFunc(db.dirac_en)
rp = db.MakeRunProgram('pam-dirac --inp=dir.inp',db.transform_input_diatomic)
print "transform_input_diatomic_test"
rp('PbO.xyz','in.mol',2.0)
print "en test"
print energy('dirac_out')

print "transform_input_eval_test"
rp = db.MakeRunProgram('pam-dirac --inp=dir.inp',db.transform_input_eval)
rp('PbF4.opt','in.mol',2.0)
