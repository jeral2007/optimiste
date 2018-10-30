from cry_out import CrystallOut
import scipy as sc
import sys
outfile = sys.argv[1]
dp1, dp2 = map(int, sys.argv[2:])
cry_out = CrystallOut(outfile)
print "sum gradient"
print '-'*20
cs, ats = cry_out.primitive_cell()
N = len(ats)
grads = cry_out.cart_grads()
sum_grads = []
sum_ats = []
dp1_inds = []
for ii in xrange(N):
    if ats[ii] == dp1:
        sum_grads += [grads[:, ii]]
        dp1_inds += [len(sum_grads)-1]
        sum_ats += [ats[ii]]
    elif ats[ii] == dp2:
        ind = dp1_inds.pop(0)
        sum_grads[ind] += grads[:, ii]
    else:
        sum_grads += [grads[:, ii]]
        sum_ats += [ats[ii]]
rms = sum(sc.dot(v, v) for v in sum_grads)/len(sum_grads)
for v, at in zip(sum_grads, sum_ats):
    print "{0} {1[0]} {1[1]} {1[2]}".format(at, v)
print "RMS GRAD: {}".format(rms**0.5)

