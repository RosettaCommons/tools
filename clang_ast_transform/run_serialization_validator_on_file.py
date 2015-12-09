import os
import sys
import subprocess
import json

# To run on this script on all cc files in the Rosetta source tree, use the following command from within
# the Rosetta/main/source/src directory:
#
# python /path/to/tools/clang_ast_transform/run_on_all_ccfiles_w_fork.py -e "python /path/to/tools/clang_ast_transform/run_serialization_validator_on_file.py" -n 10
#

# Assumptions for how this script is run below:
# 1: assume that this is being run from somewhere within the rosetta "main" repository
# 2: assume that this script lives in the Rosetta/tools/clang_ast_transform directory.
#    The rosetta-refactor-tool must live in the
#    Rosetta/tools/clang_ast_transform/clang/build/bin
#    directory

# $1 == the file to run the serialization validator on

# get the directory where this script lives
# CLANG_AST_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
script_dir = os.path.realpath( __file__ ).rpartition("/run_serialization_validator_on_file.py")[0]
# CLANG_BIN=$CLANG_AST_DIR/clang/build/bin
clang_exec_dir = script_dir + "/clang/build/bin"

# get the directory where we're executing this script
# SOURCE=$( pwd | sed 's:/src/: :' | awk '{print $1}' )
rosetta_source_dir = os.getcwd().partition( "/src/" )[0]

# the file to "compile" should be the one and only argument to this script.
#MATCHERS=$1
# FILE=$1
fname = sys.argv[1]
json_outfname = "serialization_validator_" + fname.replace("/","_") + ".json"

#echo "matchers", $MATCHERS

#cd $SOURCE
os.chdir( rosetta_source_dir )
command = clang_exec_dir + "/serialization_validator " + fname + " -- " + \
          "clang++ -D__STDC_CONSTANT_MACROS -D__STDC_LIMIT_MACROS " + \
          "-std=c++11 " + \
          "-pipe " + \
          "-ffor-scope " + \
          "-DEXIT_THROWS_EXCEPTION " + \
          "-isystem " + rosetta_source_dir + "/external/boost_1_55_0/ " + \
          "-isystem " + rosetta_source_dir + "/external/include/ " + \
          "-isystem " + rosetta_source_dir + "/external/dbio/ " + \
          "-DUNUSUAL_ALLOCATOR_DECLARATION " + \
          "-stdlib=libstdc++ " + \
          "-DBOOST_ERROR_CODE_HEADER_ONLY " + \
          "-DBOOST_SYSTEM_NO_DEPRECATED " + \
          "-DNDEBUG " + \
          "-DPTR_STD -DCXX11 " + \
          "-DSERIALIZATION " + \
          "-I" + rosetta_source_dir + "/src " + \
          "-I" + rosetta_source_dir + "/external/include " + \
          "-I" + rosetta_source_dir + "/src/platform/linux/64/clang/3.5-1ubuntu1 " + \
          "-I" + rosetta_source_dir + "/src/platform/linux/64/clang " + \
          "-I" + rosetta_source_dir + "/src/platform/linux/64 " + \
          "-I" + rosetta_source_dir + "/src/platform/linux " + \
          "-I" + rosetta_source_dir + "/external/boost_1_55_0 " + \
          "-I" + rosetta_source_dir + "/external/dbio " + \
          "-I/usr/include " + \
          "-I/usr/local/include"
command_list = command.split()
#print command_list
p = subprocess.Popen( command_list, stdout = subprocess.PIPE )
output = p.stdout.read()
outlines = []

outdict = {}
outdict[ fname ] = {}
outdict[ fname ][ "log" ] = "Human-readable log for " + fname
outdict[ fname ][ "state" ] = "passed"

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

outlines = json.dumps( outdict )
open( json_outfname, "w" ).writelines( outlines )

sys.exit( 1 if output else 0 )
