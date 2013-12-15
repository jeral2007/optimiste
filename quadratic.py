#!/usr/bin/env python 
#This is quadratic optimisation one dimensional alghoritm realization
#Lomachuk Y.V, first version on 9 Dec 2013 
import scipy,scipy.linalg
#Linear dependency error
class LinDep(Exception):
    def __init__(self,mat):
            self.mat = mat
#Not Minimum error
class NotMin(Exception):
        def __init__(self,coefs):
            self.coefs = coefs
class PointsTooNear(Exception):
    pass

def quadr_optimize(func,x0vec, eps=1e-5, maxiter = 40):
    """ quadr_optimize takes callable func(x), that takes one float value and return float value,
and set of 3 initial points x0vec and returns value of the minimum, 

optional parameters:
eps -- absolute accuracy to achieve
maxiter -- maximum iterations count"""
    # check and setup values
    if (abs(x0vec[0]-x0vec[1])<eps or abs(x0vec[0]-x0vec[2])<eps):raise PointsTooNear
    assert(eps>0)
    xvec = x0vec
    xmin = xvec[0]
    yval = scipy.array([func(x) for x in xvec])
    curiter = 1
    dx = eps*2
    while(curiter<maxiter and dx>eps):
        xmat = scipy.array([[x**i for i in xrange(3)] for x in xvec])
        #find coefs array coefs_i is the coefficent value for the x**i term
        try:
            coefs =scipy.linalg.solve(xmat,yval)
        except scipy.linalg.LinAlgError:
            print (xmat,oxmin,xmin)
            raise LinDep(xmat)
        
        if coefs[2]<=0: # if f(x) is nearly a_2 x**2 + a_1 x + a_0 and a[2]<0, there is no minimum. 
            raise NotMin(coefs)

        oxmin,xmin = xmin,-coefs[1]/(2*coefs[2]) #update xmin value and store old of this in oxmin.
        #update xvec and yvec
        ind1 = yval.argmax()
        xvec[ind1] = xmin
        yval[ind1] = func(xmin)
        # stop if minimum position remain unchanged or distance between points lesser than eps
        dx = min(abs(oxmin-xmin),abs(xvec[0]-xvec[1]),abs(xvec[0]-xvec[2]),
                        abs(xvec[1]-xvec[2]))
        curiter+=1

    return xmin # FIXME with maxiter account 

