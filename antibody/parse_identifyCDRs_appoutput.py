#Script to for parsing pdb files for CDR data
#author: Laura Beth Fulton


import re
import argparse

# def str2dict(filename = "antibodycluster.raw"):
#   antibody_dict ={}
#   keys=[]
    #antibody_dict={'key': 'info'}
    
#reading in file into list of lines
with open("identify_cdrsapp.out", 'r') as my_file:
    lines = my_file.readlines()

# current key
c_key=""

# results
results_dict={}
errors_dict={}

# temporary string
remarks = []

#e_key=""
#temp string
eremarks= []
#error_dict= {}
#errors

    
#loop through lines to get keys
#for line in lines:
for i in range(len(lines)):

    line = lines[i]

    # first pdb ever!
    if line.startswith("protocols.jd2.PDBJobInputter: filling pose from PDB"):
        # new pdb detected
        if c_key != '': # not the first time I have found a pdb
            # do stuff for the old key
            results_dict[c_key] = remarks
            remarks = []

        c_key = line[-17:-13]   #attempt to get the '12e8 portion of pdb'
    # any other pdb!
    # lines in between
    elif line.startswith('REMARK'):
        # hey, maybe it makes sense to split thhis line up and parse out the clutster.
        line_components = line.rstrip().split()
        cluster_id = line_components[2]
        cluster_distance = line_components[3]
        remarks.append('{} {}'.format(cluster_id, cluster_distance))

    elif line.startswith("[ERROR]"):
        err_lines = [line.strip()]
        while not line.startswith("protocols"):
            i += 1
            line=lines[i]
            err_lines.append(line.strip())

            errors_dict[c_key] = err_lines[1:-1]


#print [x for x in results_dict['12e8'].split('\n') if x.startswith("REMARK CLUSTER")]

#for line in lines:
#    if line.startswith("[ERROR] Exception caught by JobsDistributor for job " ):
#            if e_key != '':
#                results_dict[e_key] = eremarks
#                eremarks = []
#            e_key= line

    #important for all lines
    #for line in lines:
        #if line.startswith("REMARK CLUSTER"):
            #results_dict.split('\n')

#print lines
print "PDB"
#print "PDB, H1, H2, H3, L1, L2, L3"
for pdb, remark_lines in results_dict.iteritems():
    #print '{}:\n{}\n'.format(pdb, '\n'.join(remark_lines))
   

#    print pdb
#    print remark_lines

    if len(remark_lines) == 0:
        print '{} | \"{}\"'.format(pdb, ' '.join([x for x in errors_dict[pdb]]))

    #print  '{},{}'.format(pdb, e_key)

    #     continue
    # print '{}, {}'.format(pdb, ",".join([x.split()[0] for x in remark_lines]))