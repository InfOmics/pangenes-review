#!/usr/bin/python3

import sys

ifile = sys.argv[1]
ofile = sys.argv[2]

lineno = 0
with open(ofile,'w') as of:
	for line in open(ifile, 'r'):
		if lineno % 2 == 0:
			of.write(">"+line)
		else:
			of.write(line)
	
		lineno += 1
