class table_pp(object):
    def __init__(self, cell_specs, trows_per_page=100000, cols_per_page=1, header=""):
        self.cell_specs = cell_specs
        self.trows_per_page = trows_per_page
        self.cols_per_page = cols_per_page
        self.cur_row = 0
        self.cur_col = 0
        self.table = []
        self.changed = True
    def append(self, arr):
        res = []
        for el, cell_spec in zip(arr, self.cell_specs):
            aux =  = cell_spec['fmt'].format(el)
            maxlength = cell_spec['length']
            align = cell_spec['align']
            res += [""]
            for rcell in aux.split('\n'):
                raw_cell = rcell.strip()
                if len(raw_cell) % maxlength == 0:
                    nrows = len(raw_cell) / maxlength - 1
                    lastcol = maxlength
                    lrskip = (0, 0)
                else:  
                    nrows = len(raw_cell) / maxlength
                    lastcol = len(raw_cell) % maxlength
                    if align == 'l':
                        lrskip = (0, maxlength - lastcol)
                    elif align == 'r':
                        lrskip = (maxlength - lastcol, 0)
                    else:
                        t = (maxlength - lastcol) / 2
                        lrskip = (t, maxlength - t)
                for row in range(nrows):
                    res[-1] += [raw_cell[row*maxlength:(row+1)*maxlength]]
                    res[-1] += [" "*lrskip[0]+raw_cell[-lastcol:]+" "*lrskip[1]]
            self.table +=[res]
       def __process_x__(self):
           self.changed == False
           if 
       def __process_y__(self):
            self.changed = False
            for row in self.table:
                max_disp_rows = max(len(c) for c in row)
                if all(len(c) == max_disp_rows for c in row):
                    continue
                self.changed = True
                for n in range(len(row)):
                    cell_spec = self.cell_specs[n]
                    if len(row[n]) 
                    
             
