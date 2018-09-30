"""\
Bunch of functions for evaluating minimum energy path
"""
import scipy as sc
import scipy.linalg as lin


def force(f0s, xs, k):
    """\
evaluates forces on the images on path (see DOI: 10.1063/1.2841941) for
details, INPUT: f0 -- potential forces, xs -- positions of the images, k -
string stifnes OUTPUT: res -- resulting force.

f0, xi are scipy float arrays with shape (Npoints, Ndim), where Npoints is
number of images, Ndim is the dimensionality of configurational space , k -
scalar float value.  resulting force res is the sum of component of potential
force f0 transversal component and linear force that is tangential to the
string of images."""
    assert(f0s.shape == xs.shape)
    drs = sc.zeros_like(xs)
    drs[:-1] = xs[1:] - xs[:-1]
    fs = sc.zeros_like(f0s)
    for i in xrange(1, len(f0s) - 1):
        tau = drs[i]/lin.norm(drs[i])
        ftan = tau*k*(lin.norm(drs[i]-drs[i-1]) - lin.norm(drs[i+1]-drs[i]))
        fs[i] = f0s[i] - sc.dot(tau, f0s[i])*tau - ftan
    return fs


def lin_interp(xs, xg):
    return sc.array([sc.argmin((x-xg)**2) for x in xs])
