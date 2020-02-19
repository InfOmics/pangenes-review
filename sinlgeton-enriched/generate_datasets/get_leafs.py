#!/usr/bin/python3


import sys
from statistics import *

iparents = sys.argv[1]
iseqs = sys.argv[2]
nof_leafs = int(sys.argv[3])
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


def depth_distribution(n, depth, childs, distr, depths):
	if n in childs:
		for c in childs[n]:
			distr[depth + 1] = distr.get(depth + 1, 0) + 1
			depths[c] = depth + 1
			depth_distribution(c, depth +1, childs, distr, depths)

depth_distr = dict()
depths = {0:0}
depth_distribution(0,0,childs, depth_distr, depths)
print('-'*40)
print("nof genomes per depth")
for k,v in sorted(depth_distr.items()):
	print(k,v)


print('-'*40)
vv = [ len(x) for x in childs.values() ]
print("average/stdev degree", mean(vv), stdev(vv))

print('-'*40)
print('degree by depth')
print(0, len(childs[0]), 0)
for k in sorted(depth_distr.keys()):
	vv = [ len( childs[x]) for x in parent_g if depths[x] == k ]
	if len(vv) > 1:
		print(k, mean(vv), stdev(vv))
	elif len(vv) > 0:
		print(k, mean(vv), 0)

print('-'*40)
vv = sorted([ (depths[g] ,g) for g in ( all_g - parent_g  ) ])
#print(vv)
leafs = [v[1] for v in vv[- nof_leafs :]]


with open(ofile, 'w') as off:
	i = 0
	take = False
	sdescr = ""
	for line in open(iseqs, 'r'):
		if i%2 == 0:
			g = int( (line.split('\t')[0]).split('_')[1] )
			if g in leafs:
				take = True
				sdescr = line 
		else:
			if take:
				off.write(sdescr)
				off.write(line)
			take = False
		i += 1

