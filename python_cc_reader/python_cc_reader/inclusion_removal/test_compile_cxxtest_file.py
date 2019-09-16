from .test_compile import *
import sys

# $1 cxxtest file to test compilation of
cxxfname = sys.argv[1]
if cxxtest_test_compile(cxxfname, True):
    print("success")
else:
    print("failed")
