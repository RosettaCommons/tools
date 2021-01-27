#!/bin/bash

CLANG_BIN=/local/luki/clang/build/bin
SOURCE=/local/luki/main/source
OUT_DIR=/local/luki/main-copy/source
FILE=$1

cd $SOURCE

cp -a $FILE $OUT_DIR/$FILE

$CLANG_BIN/rosetta-refactor-tool -matchers=rewrite $OUT_DIR $FILE -- \
	clang++ -D__STDC_CONSTANT_MACROS -D__STDC_LIMIT_MACROS \
	-std=c++11 \
	-isystem external/boost_submod/ \
	-isystem external/ \
	-isystem external/include/ \
	-isystem external/dbio/ \
	-isystem external/libxml2/include \
	-isystem external/rdkit/ \
	-DUNUSUAL_ALLOCATOR_DECLARATION \
	-stdlib=libstdc++ \
	-DBOOST_ERROR_CODE_HEADER_ONLY \
	-DBOOST_SYSTEM_NO_DEPRECATED \
	-DNDEBUG \
	-DPTR_REFCOUNT \
	-Isrc \
	-Iexternal \
	-Iexternal/include \
	-Isrc/platform/linux/64/clang/3.5-1ubuntu1 \
	-Isrc/platform/linux/64/clang \
	-Isrc/platform/linux/64 \
	-Isrc/platform/linux \
	-Iexternal/dbio \
	-Iexternal/cxxtest/ \
	-Itest/ \
	-I./ \
	-I/usr/include \
	-I/usr/local/include
