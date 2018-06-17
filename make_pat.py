#!/usr/bin/env python
import sys
q1 = float(sys.argv[1])
q2 = float(sys.argv[2])
sc = 1+ float(sys.argv[3])
f = open('coord_pat_pat')
header = f.next()
print header.strip()
assert('$coord' in header)
state = 1
for line in f:
  if '!' in line:
     aux = line.split('!')
     line2 =  "{} {}". format(aux[0].strip(), eval(aux[1]))
  else:
     line2 =  line.strip()
  if '$end' in line2:
     state = 2
  if state == 2:
     print line2
     continue
  if 'sc' == line2[:2]:
      aux = line2[2:].split()
      aux[:3] = [sc*float(x) for x in aux[:3]]
      print " ".join(str(x) for x  in aux)
  else:
      print line2 
   
