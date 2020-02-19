#!/usr/bin/python3


import sys
from statistics import *

iparents = sys.argv[1]
iseqs = sys.argv[2]
nof_roots = int(sys.argv[3])
ofile = sys.argv[4]


childs = dict()

parent_g =set()
all_g = set()

for line in open(iparents, 'r'):
	cc = line.strip().split(' ')
	print(cc)
	p = int(cc[2])
	c = int(cc[1])
	parent_g.add(p)
	all_g.add(c)
	all_g.add(p)
	if not ( p in childs ):
		childs[p] = list()
	childs[p].append(c)

print('-'*40)
print("nof all", len(all_g))
print("nof parents", len(parent_g))
print("nof leafs",  len( all_g - parent_g  ))

roots = [0]
for i in range(0, nof_roots):
	gi = roots[i]
	if gi in childs:
		for c in childs[gi]:
			roots.append(c)
roots = roots[:nof_roots]


with open(ofile, 'w') as off:
	i = 0
	take = False
	sdescr = ""
	for line in open(iseqs, 'r'):
		if i%2 == 0:
			g = int( (line.split('\t')[0]).split('_')[1] )
			if g in roots:
				take = True
				sdescr = line 
		else:
			if take:
				off.write(sdescr)
				off.write(line)
			take = False
		i += 1

