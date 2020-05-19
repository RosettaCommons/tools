#!/bin/bash

# To run on this script on all headers in the Rosetta source tree, use the following command from within
# the Rosetta/main/source/src directory:
#
# python /path/to/tools/clang_ast_transform/run_on_all_headers_w_fork.py -e /path/to/tools/clang_ast_transform/extract_serialization_data.sh -n 4
#
# Once this has completed, cat all .def files into one file, and pass it to
# make_serialize_templates.py to insert serialization stubs

# 1: assume that this is being run from somewhere within the rosetta source tree
# 2: assume that this script lives in the Rosetta/tools/clang_ast_transform directory.
#    The rosetta-refactor-tool must live in the
#    Rosetta/tools/clang_ast_transform/clang/build/bin
#    directory

# $1 == the comma separated list of the matchers to run on a particular file
# $2 == the file to run the matchers on

# get the directory where this script lives
CLANG_AST_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
CLANG_BIN=$CLANG_AST_DIR/clang/build_debug/bin

# get the directory where we're executing this script
SOURCE=$( pwd | sed 's:/src/: :' | awk '{print $1}' )

# the file to "compile" should be the one and only argument to this script.
MATCHERS=$1
FILE=$2

echo "matchers", $MATCHERS

cd $SOURCE

$CLANG_BIN/rosetta-refactor-tool -matchers=$MATCHERS $FILE -- \
	clang++ -D__STDC_CONSTANT_MACROS -D__STDC_LIMIT_MACROS \
	-std=c++11 \
	-isystem $SOURCE/external/boost_submod/ \
	-isystem $SOURCE/external \
	-isystem $SOURCE/external/include/ \
	-isystem $SOURCE/external/dbio/ \
	-isystem $SOURCE/external/libxml2/include/ \
	-isystem $SOURCE/external/rdkit/ \
	-DUNUSUAL_ALLOCATOR_DECLARATION \
	-stdlib=libstdc++ \
	-DBOOST_ERROR_CODE_HEADER_ONLY \
	-DBOOST_SYSTEM_NO_DEPRECATED \
	-DNDEBUG \
	-DPTR_MODERN -DPTR_STD -DCXX11 \
	-DSERIALIZATION \
	-I$SOURCE/src \
	-I$SOURCE/external \
	-I$SOURCE/external/include \
	-I$SOURCE/src/platform/linux/64/clang/3.5-1ubuntu1 \
	-I$SOURCE/src/platform/linux/64/clang \
	-I$SOURCE/src/platform/linux/64 \
	-I$SOURCE/src/platform/linux \
	-I$SOURCE/external/dbio \
	-I/usr/include \
	-I/usr/local/include
