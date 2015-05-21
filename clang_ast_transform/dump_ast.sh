BD=/home/andrew/GIT/Rosetta/tools/clang_ast_transform/clang/build
SRC=/home/andrew/GIT/Rosetta/tools/clang_ast_transform/clang/llvm

ROSETTA_SOURCE=/home/andrew/GIT/Rosetta/main/source
#cd /home/andrew/GIT/Rosetta/main/source

$BD/bin/clang-check -ast-dump $1 -ast-dump-filter=$2 -- \
	clang++ -D__STDC_CONSTANT_MACROS -D__STDC_LIMIT_MACROS \
	-std=c++11 \
	-isystem $ROSETTA_SOURCE/external/boost_1_55_0/ \
	-isystem $ROSETTA_SOURCE/external/include/ \
	-isystem $ROSETTA_SOURCE/external/dbio/ \
	-DUNUSUAL_ALLOCATOR_DECLARATION \
	-stdlib=libstdc++ \
	-DBOOST_ERROR_CODE_HEADER_ONLY \
	-DBOOST_SYSTEM_NO_DEPRECATED \
	-DNDEBUG \
	-DPTR_MODERN -DPTR_BOOST \
	-I$ROSETTA_SOURCE/src \
	-I$ROSETTA_SOURCE/external/include \
	-I$ROSETTA_SOURCE/src/platform/linux/64/clang/3.5-1ubuntu1 \
	-I$ROSETTA_SOURCE/src/platform/linux/64/clang \
	-I$ROSETTA_SOURCE/src/platform/linux/64 \
	-I$ROSETTA_SOURCE/src/platform/linux \
	-I$ROSETTA_SOURCE/external/boost_1_55_0 \
	-I$ROSETTA_SOURCE/external/dbio \
	-I/usr/include \
	-I/usr/local/include
