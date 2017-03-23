#!/usr/bin/env python

from beautify_compiled_files_w_fork import *
import subprocess, sys
try:
    import blargs
except ImportError:
    # if this script is in the Rosetta/tools/xsd_xrw/ directory
    # blargs is in the ../external/ directory. Add that to the path. and re-import
    blargs_path = os.path.join( os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'external')
    sys.path.append(blargs_path)
    import blargs


if __name__ == "__main__" :
    with blargs.Parser(locals()) as p :
        p.int( "num_cpu" ).shorthand("j").default(1)
        p.str( "pound_if_setting" ).default("take_if").described_as( "Can be take_if or take_else.")
        p.flag( "dry_run" ).described_as( "Only attempt to beautify without changing the files.")
        p.str( "ref_branch" ).default("origin/master").described_as( "Which branch to use for 'master' when calculating changes." )

    closest_master_command = [ "git", "merge-base", ref_branch, "HEAD" ]
    rev_to_diff_against = subprocess.Popen(closest_master_command, stdout=subprocess.PIPE).communicate()[0].strip()
    if ( rev_to_diff_against == '' ): # If, for some reason, the branch doesn't exist
        sys.stderr.write( "ERROR: Branch '" + ref_branch + "' doesn't seem to be a valid branch in this repository - not beautifying.\n" )
        sys.exit(-1)
    #print "rev to diff: " + rev_to_diff_against
    bash_command = [ "git", "diff", "--relative", "--name-only", rev_to_diff_against, "HEAD" ]
    file_list = subprocess.Popen(bash_command, stdout=subprocess.PIPE).communicate()[0]
    file_list = file_list.splitlines()
    file_list = [ x for x in file_list if x[-3:]==".hh" or x[-3:]==".cc" ]
    # print file_list
    # sys.exit(0)

    fbm = beautify_files_in_parallel( file_list, not dry_run, num_cpu, pound_if_setting )
    exit_following_beautification( fbm )

