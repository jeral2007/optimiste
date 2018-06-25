from scipy import float64
def read_coord(filename):
    f = open(filename)
    res = []
    header = f.next()
    if '$coord' not in header:
        raise ValueError
    for line in f:
        if '$end' in line:
            return res
        res += [map(float64, line.split()[:-1]) + [line.split()[-1]]]

def coord_iterator(filename):
    f = open(filename)
    res = []
    header = f.next()
    if '$coord' not in header:
        raise ValueError('Invalid coord file {}'.format(filename))
    for line in f:
        if '$end' in line:
             break
        res= map(float64, line.split()[:3]) + [line.split()[3]]
        yield res
    f.close()


def table(arr, maxcol,field_fmt="{}", sep="\t"):
    N = len(arr)
    if N % maxcol == 0:
        nrows = N/maxcol - 1
        lastcol = maxcol
    else:
        nrows = N / maxcol
        lastcol = N % maxcol
    fmt_row = sep.join([field_fmt]*maxcol)
    fmt_last_row = sep.join([field_fmt]*lastcol)
    res=""
    for srow in range(nrows):
        res += fmt_row.format(*arr[srow*maxcol:(srow+1)*maxcol])+"\n"
    res += fmt_last_row.format(*arr[N-lastcol:])
    return res

