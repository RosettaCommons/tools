#!/bin/bash

CLANG_BIN=/local/luki/clang/build/bin
SOURCE=/local/luki/main/source
FILE=$1

$CLANG_BIN/rosetta-refactor-tool ${*:2} $FILE -- \
	clang++ -D__STDC_CONSTANT_MACROS -D__STDC_LIMIT_MACROS \
	-std=c++11 \
	-isystem $SOURCE/external/boost_submod/ \
	-isystem $SOURCE/external/ \
	-isystem $SOURCE/external/include/ \
	-isystem $SOURCE/external/dbio/ \
	-isystem $SOURCE/external/libxml2/include \
	-isystem $SOURCE/external/rdkit/ \
	-DUNUSUAL_ALLOCATOR_DECLARATION \
	-stdlib=libstdc++ \
	-DBOOST_ERROR_CODE_HEADER_ONLY \
	-DBOOST_SYSTEM_NO_DEPRECATED \
	-DNDEBUG \
	-DPTR_MODERN -DPTR_STD -DCXX11 \
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
