#!/bin/bash

# Run on all headers: python run_on_all_headers_w_fork.py -e $PWD/extract_serialization_data.sh -n 4
# Then cat all .def files into one file, and pass it to make_serialize_templates.py to insert serialization stubs

CLANG_BIN=/local/luki/clang/build/bin
SOURCE=/local/luki/main/source
FILE=$1

cd $SOURCE
$CLANG_BIN/rosetta-refactor-tool -matchers=find_record_decl,find_constructor_decl,find_field_decl $FILE -- \
	clang++ -D__STDC_CONSTANT_MACROS -D__STDC_LIMIT_MACROS \
	-std=c++11 \
	-isystem $SOURCE/external/boost_1_55_0/ \
	-isystem $SOURCE/external/include/ \
	-isystem $SOURCE/external/dbio/ \
	-DUNUSUAL_ALLOCATOR_DECLARATION \
	-stdlib=libstdc++ \
	-DBOOST_ERROR_CODE_HEADER_ONLY \
	-DBOOST_SYSTEM_NO_DEPRECATED \
	-DNDEBUG \
	-DPTR_MODERN -DPTR_STD -DCXX11 \
	-I$SOURCE/src \
	-I$SOURCE/external/include \
	-I$SOURCE/src/platform/linux/64/clang/3.5-1ubuntu1 \
	-I$SOURCE/src/platform/linux/64/clang \
	-I$SOURCE/src/platform/linux/64 \
	-I$SOURCE/src/platform/linux \
	-I$SOURCE/external/boost_1_55_0 \
	-I$SOURCE/external/dbio \
	-I/usr/include \
	-I/usr/local/include > $FILE.def
