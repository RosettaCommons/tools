#!/usr/bin/env python

# Run this script from somewhere beneath the source/ directory in the main repository;
# it will look for "source" in the CWD and then move into the root directory for the
# repository to do its business


from beautify_compiled_files_w_fork import *
import subprocess, sys
try:
    from python_cc_reader.external.blargs import blargs
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
        p.flag( "quiet" ).shorthand("q").described_as( "Silent all output produced by this script" )

    # ok, let's CD into the root of the git repository
    # which should house the "source" directory
    pwd = os.getcwd();
    pwd_parts = pwd.rpartition( "source" )
    assert( pwd_parts[1] == "source" ) # you must run this script from the main/source directory or one of its subdirectories
    if not quiet :
        print("cd'ing to " + pwd_parts[0])
    os.chdir( pwd_parts[0] )


    closest_master_command = [ "git", "merge-base", ref_branch, "HEAD" ]
    rev_to_diff_against = subprocess.Popen(closest_master_command, stdout=subprocess.PIPE).communicate()[0].strip()
    if ( rev_to_diff_against == '' ): # If, for some reason, the branch doesn't exist
        sys.stderr.write( "ERROR: Branch '" + ref_branch + "' doesn't seem to be a valid branch in this repository - not beautifying.\n" )
        sys.exit(-1)
    #print("rev to diff: " + rev_to_diff_against)
    # --name-status -- Simple "files changed" listing which annotates files which have been removed.
    # --no-renames -- turn the more complicated rename line into two add & delete lines
    # --relative -- make the name listing relative to the current subdirectory.
    bash_command = [ "git", "diff", "--relative", "--name-status", rev_to_diff_against, "HEAD" ]
    file_list = subprocess.Popen(bash_command, stdout=subprocess.PIPE).communicate()[0].decode('ascii')
    #print("Initial list\n", file_list)
    file_list = [str(x) for x in file_list.splitlines()]
    file_list = [ x.split(None,1) for x in file_list ]

    # pare down this list to the set of files that should be beautified at all
    # Files with status of D have been deleted, and don't need beautification.
    file_list = [ x for s, x in file_list if (x[-3:]==".hh" or x[-3:]==".cc") and s != 'D' ]
    file_list = [ x for x in file_list if x.find("source/src/ui/") == -1 and x.find( "source/code_templates/" ) == -1  ]
    if not quiet :
        print("Preparing to beautify: " + ", ".join( file_list ))
    # sys.exit(0)

    fbm = beautify_files_in_parallel( file_list, not dry_run, num_cpu, pound_if_setting, quiet )
    exit_following_beautification( fbm )
