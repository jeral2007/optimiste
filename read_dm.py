from standard_config import *
import scipy as sc


def make_ls_js_ms(header):
    import re
    ket_re = re.compile(r'\|([SPDFGH]),([0-9]+.5),(-?[0-9]+.5)\>')
    kets = header.split()
    ls, js, ms = [], [], []
    for ket in kets:
        aux = ket_re.match(ket).groups()
        ls += ['SPDFGH'.index(aux[0])]
        js += [float(aux[1])]
        ms += [float(aux[2])]

    ns = [0]
    for ii in xrange(1, len(ls)):
        if ls[ii] != ls[ii-1] or (abs(js[ii]-js[ii-1]) > 1e-3) or (ms[ii] > ms[ii-1]):
            ns += [ns[-1] + 1]
        else:
            ns += [ns[-1]]
    return ls, js, ms, ns


def read_dm(jred_dm_re, jred_dm_im):

    def my_chop(stream):
        for line in stream:
            yield line.split('>')[1]  # FIXME

    fre = open(jred_dm_re)
    fim = open(jred_dm_im)
    N = int(fre.next())
    assert(int(fim.next()) == N)
    header = fre.next()
    assert(fim.next() == header)
    dm_re = sc.loadtxt(my_chop(fre))
    dm_re.shape = (N, N)
    im_re = sc.loadtxt(my_chop(fim))
    im_re.shape = (N, N)
    return (dm_re + 1j*im_re,)+make_ls_js_ms(header)


if __name__ == '__main__':
    dm, ls, js, ms, ns = read_dm(jred_dm_re.format(Rc='0.7'), jred_dm_im.format(Rc='0.7'))
    dm = sc.real(dm)
    print '-'*20
    print dm
    print '-'*20
    print ls
    print ns
    print '-'*20
    for n in xrange(ns[-1]+1):
        inds = [ii for ii, n1 in enumerate(ns) if n1 == n]
        print dm[inds, inds]
        print "Tr |{0}>P<{0}| = {1:.5f}".format(n, sc.sum(dm[inds, inds]))
