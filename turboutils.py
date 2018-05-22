def read_coord(filename):
    f = open(filename)
    res = []
    header = f.next()
    if '$coord' not in header:
        raise ValueError
    for line in f:
        if '$end' in line:
            return res
        res += [map(float, line.split()[:-1]) + [line.split()[-1]]]

