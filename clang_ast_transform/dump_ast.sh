BD=/local/luki/clang/build
SRC=/local/luki/clang/llvm

cd /local/luki/main/source
$BD/bin/clang-check -ast-dump $1 -ast-dump-filter=$2 -- \
	clang++ -D__STDC_CONSTANT_MACROS -D__STDC_LIMIT_MACROS \
	-std=c++11 \
	-isystem external/boost_1_55_0/ \
	-isystem external/include/ \
	-isystem external/dbio/ \
	-DUNUSUAL_ALLOCATOR_DECLARATION \
	-stdlib=libstdc++ \
	-DBOOST_ERROR_CODE_HEADER_ONLY \
	-DBOOST_SYSTEM_NO_DEPRECATED \
	-DNDEBUG \
	-Isrc \
	-Iexternal/include \
	-Isrc/platform/linux/64/clang/3.5-1ubuntu1 \
	-Isrc/platform/linux/64/clang \
	-Isrc/platform/linux/64 \
	-Isrc/platform/linux \
	-Iexternal/boost_1_55_0 \
	-Iexternal/dbio \
	-I/usr/include \
	-I/usr/local/include
