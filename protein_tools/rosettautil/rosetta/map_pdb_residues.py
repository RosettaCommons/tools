import sys
import os
import tempfile
import commands
import traceback
import sqlite3
from optparse import OptionParser # deprecated since Python 2.7

# Python functions to map PDB residue IDs to Rosetta/pose IDs by using the features database.
# Direct complaints to shane.oconnor@ucsf.edu
# Warning: The inputs to the commands.getstatusoutput (one of these is an executable) are unsanitized. Only use these functions if you trust the caller.
#
# Sample command line:
#   python map_pdb_residues.py -d ~/rosetta/main/database -e ~/rosetta/main/source/bin/rosetta_scripts.static.linuxgccrelease -f 1QG8.pdb -c A

script = '''<ROSETTASCRIPTS>
  <MOVERS>
    <SavePoseMover name=init_struct reference_name=init_struct/>
    <ReportToDB name="features_reporter" database_name="%s">
        <feature name="ResidueFeatures"/>
        <feature name="PdbDataFeatures"/>
    </ReportToDB>
  </MOVERS>
  <PROTOCOLS>
    <Add mover_name=init_struct/>
    <Add mover_name=features_reporter/>
  </PROTOCOLS>
</ROSETTASCRIPTS>'''

def get_pdb_to_post_residue_map(pdb_path, rosetta_scripts_path, rosetta_database_path):
    '''Takes a path to a PDB file, the RosettaScripts executable, and the Rosetta database and then uses the features database to map PDB residue IDs to pose residue IDs.
       On success, (True, the residue mapping) is returned. On failure, (False, a list of errors) is returned.'''
    mapping = {}
    errors = []
    exit_code = 0
    F, script_path = tempfile.mkstemp(dir=".")
    script_handle = os.fdopen(F, "w")
    try:
        db_path = script_path + ".db3"
        script_handle.write(script % db_path)
        script_handle.close()
        command_line = '%s -database %s -constant_seed -ignore_unrecognized_res -in:file:s %s -parser:protocol %s -overwrite -out:nooutput' % (rosetta_scripts_path, rosetta_database_path, pdb_path, script_path)
        exit_code, stdout = commands.getstatusoutput(command_line)
        if exit_code != 0:
            errors.append("An error occured during execution. The exit code was %d. The output was:\n\n%s" % (exit_code, stdout))
        else:
            conn = sqlite3.connect(db_path)
            try:
                results = conn.cursor().execute('''
SELECT chain_id, pdb_residue_number, insertion_code, residues.struct_id, residues.resNum, residues.name3, residues.res_type
FROM residue_pdb_identification
INNER JOIN residues ON residue_pdb_identification.struct_id=residues.struct_id AND residue_pdb_identification.residue_number=residues.resNum
''')
                for r in results:
                    mapping["%s%s%s" % (r[0], str(r[1]).rjust(4), r[2])] = {'pose_residue_id' : r[4], 'name3' : r[5], 'res_type' : r[6]}
            except Exception, e:
                errors.append("The features database does not seem to have been correctly created. Check to see if the command '%s' is correct." % command_line)
    except Exception, e:
        errors.append(str(e))
        errors.append(traceback.format_exc())
        exit_code = 1

    if os.path.exists(script_path):
        os.remove(script_path)
    if os.path.exists(db_path):
        os.remove(db_path)
    if exit_code or errors:
        return False, errors
    return True, mapping

def strip_pdb(pdb_path, chains = [], strip_hetatms = False):
    '''Takes a PDB file and strips all lines except ATOM and HETATM records. If chains is specified, only those chains are kept. If strip_hetatms is True then HETATM lines are also stripped.
       Returns (True, a path to the stripped PDB file) on success and (False, a list of errors) on failure.'''
    chains = set(chains)
    contents = open(pdb_path).read().split("\n") # file handle should get garbage collected
    if strip_hetatms:
        atom_lines = [l for l in contents if l.startswith("ATOM  ") and l[21] in chains]
    else:
        atom_lines = [l for l in contents if (l.startswith("ATOM  ") or l.startswith("HETATM")) and l[21] in chains]
    existing_chains = set([l[21] for l in atom_lines])
    if chains.difference(existing_chains):
        return False, ["Error: The following chains do not exist in the PDB file - %s" % ", ".join(list(chains.difference(existing_chains)))]

    F, temp_pdb_path = tempfile.mkstemp(dir=".")
    temp_pdb_handle = os.fdopen(F, "w")
    temp_pdb_handle.write("\n".join(atom_lines))
    temp_pdb_handle.close()
    return True, temp_pdb_path

def get_stripped_pdb_to_post_residue_map(input_pdb_path, rosetta_scripts_path, rosetta_database_path, chains = [], strip_hetatms = False):
    '''Takes a path to an input PDB file, the path to the RosettaScripts executable and Rosetta database, an optional list of chains to strip the PDB down to, and an optional flag specifying whether HETATM lines should be stripped from the PDB.
       On success, a pair (True, mapping between PDB and pose residues) is returned. On failure, a pair (False, a list of errors) is returned.'''
    success, result = strip_pdb(input_pdb_path, chains = chains, strip_hetatms = strip_hetatms)
    if success:
        assert(os.path.exists(result))
        success, mapping = get_pdb_to_post_residue_map(result, rosetta_scripts_path, rosetta_database_path)
        os.remove(result)
        if success:
            return True, mapping
        else:
            return False, mapping
    return False, result

if __name__ == '__main__':
    chains = []
    
    parser = OptionParser()
    parser.add_option("-e", "--executable", dest="rosetta_scripts_path", help="The location of the RosettaScripts executable e.g. ~/bin/rosetta_scripts.linuxgccrelease", metavar="EXECUTABLE")
    parser.add_option("-d", "--database", dest="rosetta_database_path", help="The location of the Rosetta database", metavar="DATABASE")
    parser.add_option("-f", "--file", dest="filename", help="The input PDB", metavar="FILE")
    parser.add_option("-c", "--chains", dest="chains", default=[], help="A comma-separated list of chains to keep (all other chains will be discarded). The default behavior is to keep all chains.")
    parser.add_option("-s", "--strip_hetatms", dest="strip_hetatms", action="store_true", default=False, help="Use this option to strip HETATM lines from the input PDB file. The default behavior is to keep HETATM lines.")
    (options, args) = parser.parse_args()
    parser.set_usage(None)
 
    filename = options.filename
    rosetta_database_path = options.rosetta_database_path
    rosetta_scripts_path = options.rosetta_scripts_path
    chains = options.chains
    strip_hetatms = options.strip_hetatms
    if not filename:
        print("\nError: A filename must be specified.\n")
        parser.print_help()
        sys.exit(1)
    elif not(os.path.exists(filename)):
        print("\nError: File '%s' does not exist.\n" % filename)
        sys.exit(1)
    if not rosetta_database_path:
        print("\nError: The path to the Rosetta database corresponding with the RosettaScripts executable must be specified.\n")
        parser.print_help()
        sys.exit(1)
    elif not(os.path.exists(rosetta_database_path)):
        print("\nError: The path '%s' does not exist.\n" % rosetta_database_path)
        sys.exit(1)
    if not rosetta_scripts_path:
        print("\nError: The path to the RosettaScripts executable must be specified.\n")
        parser.print_help()
        sys.exit(1)
    elif not(os.path.exists(rosetta_scripts_path)):
        if os.path.exists(os.path.join(os.getcwd(), rosetta_scripts_path)):
            rosetta_scripts_path = "./%s" % os.path.join(os.getcwd(), rosetta_scripts_path)
	else:
	    print("\nError: The path '%s' does not exist.\n" % rosetta_scripts_path)
            sys.exit(1)
    rosetta_scripts_path = os.path.abspath(rosetta_scripts_path)
    if chains:
        chains = chains.split(",")
        for c in chains:
            if not len(c) == 1:
                print("\nError: Chain ID '%s' is invalid. PDB chain identifiers are one character in length.\n" % c)
                sys.exit(1)

    success, result = get_stripped_pdb_to_post_residue_map(filename, rosetta_scripts_path, rosetta_database_path, chains = chains, strip_hetatms = strip_hetatms)
    if success:
        print("{")
        for k, v in sorted(result.iteritems()):
            print("'%s': %s," % (k, v))
        print("}")
    else:
        print("\n".join(result))
        sys.exit(1)

