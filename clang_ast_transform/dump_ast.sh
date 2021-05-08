
# 1: assume that this is being run from somewhere within the rosetta source tree
# 2: assume that this script lives in the Rosetta/tools/clang_ast_transform directory.
#    The rosetta-refactor-tool must live in the
#    Rosetta/tools/clang_ast_transform/clang/build/bin
#    directory

# get the directory where this script lives
CLANG_AST_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
CLANG_BUILD=$CLANG_AST_DIR/clang/build
CLANG_LLVM=$CLANG_AST_DIR/clang/llvm

# get the directory where we're executing this script
ROSETTA_SOURCE=$( pwd | sed 's:/src/: :' | awk '{print $1}' )


#BD=/home/andrew/GIT/Rosetta/tools/clang_ast_transform/clang/build
#SRC=/home/andrew/GIT/Rosetta/tools/clang_ast_transform/clang/llvm

#ROSETTA_SOURCE=/home/andrew/GIT/Rosetta/main/source
#cd /home/andrew/GIT/Rosetta/main/source

$CLANG_BUILD/bin/clang-check -ast-dump $1 -ast-dump-filter=$2 -- \
	clang++ -D__STDC_CONSTANT_MACROS -D__STDC_LIMIT_MACROS \
	-std=c++11 \
	-isystem $ROSETTA_SOURCE/external/boost_submod/ \
	-isystem $ROSETTA_SOURCE/external/ \
	-isystem $ROSETTA_SOURCE/external/include/ \
	-isystem $ROSETTA_SOURCE/external/dbio/ \
	-isystem $ROSETTA_SOURCE/external/libxml2/include \
	-isystem $ROSETTA_SOURCE/external/rdkit/ \
	-DUNUSUAL_ALLOCATOR_DECLARATION \
	-stdlib=libstdc++ \
	-DBOOST_ERROR_CODE_HEADER_ONLY \
	-DBOOST_SYSTEM_NO_DEPRECATED \
	-DNDEBUG \
	-DPTR_MODERN -DPTR_BOOST \
	-DSERIALIZATION \
	-I$ROSETTA_SOURCE/src \
	-I$ROSETTA_SOURCE/external \
	-I$ROSETTA_SOURCE/external/include \
	-I$ROSETTA_SOURCE/src/platform/linux/64/clang/3.5-1ubuntu1 \
	-I$ROSETTA_SOURCE/src/platform/linux/64/clang \
	-I$ROSETTA_SOURCE/src/platform/linux/64 \
	-I$ROSETTA_SOURCE/src/platform/linux \
	-I$ROSETTA_SOURCE/external/dbio \
	-I/usr/include \
	-I/usr/local/include
