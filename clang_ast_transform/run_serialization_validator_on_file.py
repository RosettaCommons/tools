import os
import sys
import subprocess
import json

sys.path.insert( 0, os.path.realpath(__file__).rpartition("/")[0]+"/../external" )

import blargs

# This script may be executed on its own by passing a particular file in with the --filename flag.  It should
# be run from somewhere within the Rosetta/main repository.  Alternatively, this script may be passed to the
# run_on_all_ccfiles_w_fork.py script to test all of the .cc files that are compiled in Rosetta.
#
# To run on this script on all cc files in the Rosetta source tree, use the following command from within
# the Rosetta/main/source/src directory:
#
# python /path/to/tools/clang_ast_transform/run_on_all_ccfiles_w_fork.py -e "python /path/to/tools/clang_ast_transform/run_serialization_validator_on_file.py --filename" -n 10
#
# or if you want to pass extra arguments, just make sure that --filename is the last argument provided:
#
# python /path/to/tools/clang_ast_transform/run_on_all_ccfiles_w_fork.py -e "python /path/to/tools/clang_ast_transform/run_serialization_validator_on_file.py
# --executable_path /path/to/clang/build/bin --json_output_path /path/to/directory/for/json_output --filename" -n 10

# Assumptions for how this script is run below:
# 1: assume that this is being run from somewhere within the rosetta "main" repository
# 2: assume that this script lives in the Rosetta/tools/clang_ast_transform directory.
#    The rosetta-refactor-tool must live in the
#    Rosetta/tools/clang_ast_transform/clang/build/bin
#    directory

if __name__ == "__main__" :
    with blargs.Parser(locals()) as p :
        p.str("filename").required()
        p.str("executable_path") #the path to the serialization_validator executable; if unspecified, looks in the clang/build/bin directory
        p.str("json_output_path").default(".") #if unspecified, writes to the Rosetta/main/source directory


    # get the directory where this script lives
    # CLANG_AST_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
    # CLANG_BIN=$CLANG_AST_DIR/clang/build/bin
    if not executable_path :
        script_dir = os.path.realpath( __file__ ).rpartition("/run_serialization_validator_on_file.py")[0]
        clang_exec_dir = script_dir + "/clang/build/bin"
        serialization_validator_executable = clang_exec_dir + "/serialization_validator"
    else :
        serialization_validator_executable = executable_path



    # get the directory where we're executing this script
    # SOURCE=$( pwd | sed 's:/src/: :' | awk '{print $1}' )
    rosetta_source_dir = os.getcwd().partition( "/src/" )[0]

    # the file to "compile" should be the one and only argument to this script.
    #MATCHERS=$1
    # FILE=$1
    fname = filename
    json_outfname = json_output_path + "/" + fname.replace("/","_")[ len('src/'): ] + ".json"

    #echo "matchers", $MATCHERS

    #cd $SOURCE
    os.chdir( rosetta_source_dir )
    command =  serialization_validator_executable + " " + fname + " -- " + \
              "clang++ -D__STDC_CONSTANT_MACROS -D__STDC_LIMIT_MACROS " + \
              "-std=c++11 " + \
              "-pipe " + \
              "-ffor-scope " + \
              "-DEXIT_THROWS_EXCEPTION " + \
              "-isystem " + rosetta_source_dir + "/external/boost_submod/ " + \
              "-isystem " + rosetta_source_dir + "/external/ " + \
              "-isystem " + rosetta_source_dir + "/external/include/ " + \
              "-isystem " + rosetta_source_dir + "/external/dbio/ " + \
              "-isystem " + rosetta_source_dir + "/external/rdkit/ " + \
              "-DUNUSUAL_ALLOCATOR_DECLARATION " + \
              "-stdlib=libstdc++ " + \
              "-DBOOST_ERROR_CODE_HEADER_ONLY " + \
              "-DBOOST_SYSTEM_NO_DEPRECATED " + \
              "-DNDEBUG " + \
              "-DPTR_STD -DCXX11 " + \
              "-DSERIALIZATION " + \
              "-I" + rosetta_source_dir + "/src " + \
              "-I" + rosetta_source_dir + "/external " + \
              "-I" + rosetta_source_dir + "/external/include " + \
              "-I" + rosetta_source_dir + "/src/platform/linux/64/clang/3.5-1ubuntu1 " + \
              "-I" + rosetta_source_dir + "/src/platform/linux/64/clang " + \
              "-I" + rosetta_source_dir + "/src/platform/linux/64 " + \
              "-I" + rosetta_source_dir + "/src/platform/linux " + \
              "-I" + rosetta_source_dir + "/external/dbio " + \
              "-I" + rosetta_source_dir + "/external/libxml2/include " + \
              "-I/usr/include " + \
              "-I/usr/local/include"
    # print("command\n", command)
    command_list = command.split()
    # print command_list
    p = subprocess.Popen( command_list, stdout = subprocess.PIPE, stderr = subprocess.PIPE, encoding='utf8')

    output, errors = (str(x) for x in p.communicate())
    exit_code = p.returncode

    outlines = []

    outdict = {}
    outdict[ fname ] = {}
    outdict[ fname ][ "log" ] = "Processing file " + fname + ':\n' + output + '\n' + errors
    outdict[ fname ][ "state" ] = "failed" if exit_code else "passed"

    if output :
        lines = output.split("\n")
        for line in lines :
            #print "line:\"", line,"\""
            cols = line.split()
            if len(cols) < 10 : continue
            classname = cols[4]
            field = cols[6]
            load_save = "not_" + cols[9]
            outdict[ fname ][ "state" ] = "failed"
            outdict[ fname ][ "results" ] = {}
            if classname not in outdict[ fname ][ "results" ] :
                outdict[ fname ][ "results" ][ classname ] = {}
            if load_save not in outdict[ fname ][ "results" ][ classname ] :
                outdict[ fname ][ "results" ][ classname ][ load_save ] = []
            outdict[ fname ][ "results" ][ classname ][ load_save ].append( field )

    # print("output", output)
    # print("errors", errors)
    # print("json_outfname", json_outfname)
    with open(json_outfname, 'w') as f:
        json.dump(outdict, f, sort_keys=True, indent=2)

    sys.exit( 1 if output else 0 )
