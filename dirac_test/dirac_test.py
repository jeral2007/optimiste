#!/usr/bin/env python
import sys
sys.path+=['../'] 
import dirac_bindings as db
energy = db.MakeGetEnergyFunc(db.dirac_en)
rp = db.MakeRunProgram('pam-dirac --inp=dir.inp',db.transform_input_diatomic)
print "transform_input_diatomic_test"
with open('PbO.xyz','r') as inp:
    rp(inp,2.0)
print "en test"
with open('dirac_out','r') as out:
    print energy(out)
