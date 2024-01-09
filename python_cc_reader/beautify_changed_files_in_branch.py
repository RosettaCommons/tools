#!/usr/bin/env python

from __future__ import print_function

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


def parse_changedlist(rev_to_diff_against):
    '''Look at the result of the status diff, and return the list of files to process.'''
    # --name-status -- Simple "files changed" listing which annotates files which have been removed.
    # --no-renames -- turn the more complicated rename line into two add & delete lines
    # --relative -- make the name listing relative to the current subdirectory.
    bash_command = [ "git", "diff", "--relative", "--name-status", rev_to_diff_against, "HEAD" ]
    command_output = subprocess.Popen(bash_command, stdout=subprocess.PIPE).communicate()[0].decode('ascii')
    #print("Initial list\n", command_output)

    file_list = []
    for line in command_output.splitlines():
        entries = line.split('\t')
        mod = entries[0]
        if mod == 'D':
            # Files with status of D have been deleted, and don't need beautification.
            continue
        # For added, modified, renamed and copied files, the last filename in the line is what we want to use.
        filename = entries[-1]
        if not ( filename.endswith(".hh") or filename.endswith(".cc") ):
            # Don't bother to beautify not hh/cc files.
            continue
        if filename.find("source/src/ui/") != -1 or filename.find( "source/code_templates/" ) != -1:
            # Ignore these directories.
            continue
        file_list.append( filename )

    return file_list

if __name__ == "__main__" :
    with blargs.Parser(locals()) as p :
        p.int( "num_cpu" ).shorthand("j").default(1)
        p.str( "pound_if_setting" ).default("take_if").described_as( "Can be take_if or take_else.")
        p.flag( "dry_run" ).described_as( "Only attempt to beautify without changing the files.")
        p.str( "ref_branch" ).default("origin/main").described_as( "Which branch to use for 'main' when calculating changes." )
        p.flag( "quiet" ).shorthand("q").described_as( "Silent all output produced by this script" )

    # ok, let's CD into the root of the git repository
    # which should house the "source" directory
    pwd = os.getcwd();
    pwd_parts = pwd.rpartition( "source" )
    assert( pwd_parts[1] == "source" ) # you must run this script from the main/source directory or one of its subdirectories
    if not quiet :
        print("cd'ing to " + pwd_parts[0])
    os.chdir( pwd_parts[0] )


    closest_main_command = [ "git", "merge-base", ref_branch, "HEAD" ]
    rev_to_diff_against = subprocess.Popen(closest_main_command, stdout=subprocess.PIPE).communicate()[0].strip()
    if ( rev_to_diff_against == '' ): # If, for some reason, the branch doesn't exist
        sys.stderr.write( "ERROR: Branch '" + ref_branch + "' doesn't seem to be a valid branch in this repository - not beautifying.\n" )
        sys.exit(-1)
    #print("rev to diff: " + rev_to_diff_against)
    file_list = parse_changedlist(rev_to_diff_against)
    if not quiet :
        print("Preparing to beautify: " + ", ".join( file_list ))
    # sys.exit(0)

    fbm = beautify_files_in_parallel( file_list, not dry_run, num_cpu, pound_if_setting, quiet )
    exit_following_beautification( fbm )
