from clique_analysis import PDB, read_table
import sys

mdl = PDB(sys.argv[1])
p = dict(read_table(sys.argv[2]))

for i in range(len(mdl)):
	try:
		mdl.T()[i] = p[mdl.chainID()[i]+":"+mdl.resSeq()[i]]
	except KeyError:
		mdl.T()[i] = 0.0
	# end try
# end for

mdl.write(sys.argv[3])

