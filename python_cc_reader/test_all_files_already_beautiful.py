from beautify_compiled_files_w_fork import *
import os, sys, subprocess
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
        p.flag( "dry_run" )

    # 
    parts = os.getcwd().partition("/source");
    if parts[1] != "/source" or parts[2] != ""  :
        print("this script must be run from within the main/source directory")
        sys.exit(1)

    os.chdir( "src" )
    file_list = [ "src/" + x for x in files_in_src_to_beautify() ]
    os.chdir( "../test" )
    file_list.extend( [ "test/" + x for x in files_in_test_to_beautify() ] )
    os.chdir( ".." )
    fbm = beautify_files_in_parallel( file_list, not dry_run, num_cpu, "take_if", True )

    ok = True
    if not fbm.all_files_beautified :
        for fname in fbm.files_that_failed :
            print("File", fname, "could not be beautified")
        ok = False
        print("One of the most likely reasons your file doesn't compile is that you have used a macro (e.g. TS_ASSERT) but have not followed it with a semicolon")

    git_status_cmd = [ "git", "ls-files", "-m" ]
    modified_files = str(subprocess.Popen(git_status_cmd, stdout=subprocess.PIPE, encoding='utf8').communicate()[0])

    print("modified files:", modified_files) 
    
    if modified_files == ""  :
        sys.exit( 0 if ok else 1 )
    else :
        modified_files = modified_files.splitlines()
        print("The following files changed when run through the beautifier")
        for fname in modified_files :
            print(fname)
        sys.exit(1)
