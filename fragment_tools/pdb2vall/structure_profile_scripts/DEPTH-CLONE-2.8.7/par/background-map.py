# script to mix maps and generate background map

weights = {"ALA":14607, "ARG":7643, "ASN":6339, "ASP":9938, "CYS":1580, "GLN":6171, "GLU":11767, "GLY":12558, "HIS":5692, "ILE":9217, "LEU":15705, "LYS":12872, "MET":3199, "PHE":7058, "PRO":5895, "SER":9124, "THR":8998, "TRP":2474, "TYR":5543, "VAL":11291}
denominator = 0.0
for key in weights.keys():
	denominator = denominator + weights[key]
# end for

import sys
from clique_analysis import read_table, print_table
z = [['undef' for x in range(150)] for y in range(150)] # maximum ranges for both dimensions
survive_n = sys.argv[1]
for key in weights.keys():
	t = read_table(key+'_survive'+survive_n+'.odd.map', "", True)
	print key
	for i in range(len(t)):
		for j in range(len(t[i])):
			if t[i][j] == 'nan':
				pass
			else:
				if z[i][j] == 'undef':
					z[i][j] = t[i][j]*weights[key]
				else:
					z[i][j] = z[i][j] + t[i][j]*weights[key]
				# end if
			# end if
		# end for
	# end for
# end for

for i in range(len(z)):
	for j in range(len(z[i])):
		try:
			z[i][j] = z[i][j] / denominator
		except TypeError:
			pass
		# end try
	# end for
# end for

print_table(z, sys.argv[2])
