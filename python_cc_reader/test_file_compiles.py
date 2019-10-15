# This file will test the build
# You will have to modify three variables here for your own
# computer:  os, nbits, and compiler.

import subprocess, re, time, sys
from python_cc_reader.cpp_parser.code_utilities import expand_includes_for_file, load_source_tree
from python_cc_reader.inclusion_removal.reinterpret_objdump import relabel_sections, compare_objdump_lines
from python_cc_reader.inclusion_removal.test_compile import *
from python_cc_reader.inclusion_removal.dont_remove_include import DontRemoveInclude


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python test_compile.py <filename>")
        sys.exit(1)
    print("First testing compilation directly from .cc file")

    target = sys.argv[1]
    
    dri = DontRemoveInclude()
    objdump = None
    if target in dri.surrogates:
        print("testing compilation with surrogates")
        compiles = test_compile_w_surrogates(target, dri.surrogates[target])
    else:
        compiles, objdump = generate_objdump_for_file(target)
    print("file compiles?", compiles)
    #print(objdump)

    if objdump:
        compiles2 = test_compile_for_file_extreme(sys.argv[1], objdump)
        print("file's objdump matches?", compiles2)
    
    sys.exit( not compiles)

    # compiled = test_compile(sys.argv[1], True)
    # 
    # print("Now testing compilation using python-expanded #includes")
    # compilable_files, all_includes, file_contents = load_source_tree()
    # print("...source tree loaded")
    # if sys.argv[1] not in file_contents:
    #     print("File", sys.argv[1], "not found in source tree")
    #     sys.exit(1)
    # compiled = test_compile_from_lines(
    #     expand_includes_for_file(sys.argv[1], file_contents), verbose=True
    # )
    # if compiled:
    #     print("Success")
    # else:
    #     test_compile(sys.argv[1], True)
    #     print("Failed")
