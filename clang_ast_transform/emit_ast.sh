BD=/local/luki/clang/build
SRC=/local/luki/clang/llvm

$BD/bin/clang++ $1 -emit-ast -D__STDC_CONSTANT_MACROS -D__STDC_LIMIT_MACROS \
	-std=c++11 \
	-isystem external/boost_1_55_0/ \
	-isystem external/ \
	-isystem external/include/ \
	-isystem external/dbio/ \
	-DUNUSUAL_ALLOCATOR_DECLARATION \
	-stdlib=libstdc++ \
	-DBOOST_ERROR_CODE_HEADER_ONLY \
	-DBOOST_SYSTEM_NO_DEPRECATED \
	-DNDEBUG \
	-Isrc \
	-Iexternal \
	-Iexternal/include \
	-Isrc/platform/linux/64/clang/3.5-1ubuntu1 \
	-Isrc/platform/linux/64/clang \
	-Isrc/platform/linux/64 \
	-Isrc/platform/linux \
	-Iexternal/boost_1_55_0 \
	-Iexternal/dbio \
	-I/usr/include \
	-I/usr/local/include
