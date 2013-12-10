#!/usr/bin/python -W ignore
# -*- coding: utf-8 -*-
from optparse import OptionParser 
from scipy import optimize
import re, os
from numpy import array, array_str
import copy, math

parser = OptionParser()
parser.add_option("-i", "--in", dest="infile",
                  help="file with paremetres", metavar="FILE")
(options, args) = parser.parse_args()
infile_fh = open(options.infile)
params = {}
for line in infile_fh:
	line = line.strip()
	if re.search("^#",line) or line == '':
		continue
	(key, val) = re.split("\s", line, maxsplit=1)
	params[key] = val
infile_fh.close()
#смотрю, сколько параметров
params['init'] = params['init'].strip()
init = []
init = re.split("\s", params['init'])
n = init.pop(0)
params['mode'] = params['mode'].strip()
print params['mode']
if params['mode'] == 'lr':
	for i in range(len(init)):
		init[i] = float(init[i])
		init[i] = math.log(init[i])
x0 =array(init, dtype=float)
params['step'] = float(params['step'].strip())
step_arr = [] 
step_arr.append(params['step'])
step_arr.append(params['mode'])
step =array(step_arr)
def func(x,*args):
	mode = args[1]
	if mode == 'lr':
		e = copy.deepcopy(x)
		e[len(e)-1] = float(e[len(e) - 1])
		for i in reversed(range(len(e) - 1)):
			e[i] = float(e[i])
			#e[i]=e[i+1]+e[i]
		for i in range(len(e)):
			e[i] =math.exp(e[i])
		x = copy.deepcopy(e)
	print x
	xx_xx = open('xx.xx', 'w')
	par = array_str(x)
	g = re.match(r"\[(.*?)\]", par)
 	new_g = g.group(1)
	new_g = new_g.strip()
	xx_xx.write(new_g)
	xx_xx.close()
	os.system('nrj.pl {}'.format(options.infile))
	e_file = open('xx.yy')
	str = e_file.read()
	e_file.close()
	E = re.split("\n", str)[1]
	E = float(E)
	print 'E = ', E
	return E
def grad_wul(x,*args):
	grad_file = open('xx.grad')
	str = grad_file.read()
	str = str.strip()
	grad_str = re.split("\s", str)
	myg = array(grad_str, dtype=float)
        return myg
def grad_zai(x, *args):
	step = args[0]
	step = float(step)
	grad = []
	for i in range(len(x)):
		y = copy.deepcopy(x)
		print y[i] + 0.5*step
		y[i] = y[i] + 0.5*step
		point1 = func(y, *args)
		y[i] = y[i] - step
		point2 = func(y, *args)
		grad.append((point1 -point2)/step)
	myg = array(grad, dtype=float)
	print "grad: ", myg
        return myg
# как вычислять градиенты
grad_mode = params['grad']
grad_mode = grad_mode.strip()
if grad_mode == 'y':
	grad = grad_wul
if grad_mode == 'n':
	grad = grad_zai

#Broyden-Fletcher-Goldfarb-Shanno optimization routine
#retval = optimize.fmin_cg(func, x0, grad, maxiter=10, full_output=True, disp=True)
# Powell (direction set) optimization
retval = optimize.fmin_powell(func, x0, maxiter=10, full_output=True, disp=True)
# Nelder-Mead simplex algorithm
#retval = optimize.fmin(func, x0, maxiter=10, full_output=True, disp=True)
#line-search Newton conjugate gradient optimization routine
##retval = optimize.fmin_ncg(func, x0, grad, maxiter=10, full_output=True, disp=True)
# limited-memory bound-constrained BFGS algorithm
#retval = optimize.fmin_l_bfgs_b(func, x0, grad, maxfun=100)
#retval = optimize.fmin_l_bfgs_b(func, x0, grad, args=(step), maxfun= 10)
(params, fopt, d) = retval
print '\nMy report:'
print 'X: ',retval[0]
print 'F: ',retval[1]
print 'N calls: ', d['funcalls']


