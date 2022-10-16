# This file will test the build
# You will have to modify three variables here for your own
# computer:  os, nbits, and compiler.

import subprocess, re, time, sys
from ..cpp_parser.code_utilities import expand_includes_for_file, load_source_tree
from .reinterpret_objdump import relabel_sections, compare_objdump_lines

EXTERNAL_INCLUDE_DIRS = [
    'external',
    'external/include',
    'external/boost_submod',
    'external/dbio',
    'external/libxml2/include',
    'external/rdkit',
]

def no_empty_args(command_list):
    clprime = []
    for arg in command_list:
        if arg != "":
            clprime.append(arg)
    return clprime


def central_compile_command():

    # The os string is the name of your operating system.  This is used for the src/platform includes.
    # check what scons uses if you're uncertain.

    # os = "macos"
    os = "linux"
    # os = "cygwin"
    # os = "windows"

    # The nbits string is the number of bits used in pointers in your system; in 99%
    # of cases, it will be either 64 or 32

    nbits = "64"
    # nbits = "32"

    # The compiler string is the executable that compiles your code.
    compiler = "g++"
    # compiler = "g++-mp-4.3"
    # compiler = "mpic++"

    # -S for g++ avoids the code-generation step in complation.
    # if your compiler doesn't support -S, you might see this error:
    # cc1plus: error: output filename specified twice

    include_directories = (
        ' '.join( "-isystem ../"+d for d in EXTERNAL_INCLUDE_DIRS if d != "external/libxml2/include" )
        # weird warning from the libxml2 -isystem line:
        # "g++: warning: ../external/libxml2/include: linker input
        # file unused because linking not done"
        + " -I./ -I../external -I../external/include -Iplatform/"
        + os
        + "/"
        + nbits
        + "/gcc -Iplatform/"
        + os
        + "/"
        + nbits
        + " -Iplatform/"
        + os
        + ' ' + ' '.join( "-I../"+d for d in EXTERNAL_INCLUDE_DIRS )
        + " -I/usr/local/include -I/usr/include/ "
    )

    generic_command = (
        " -c -std=c++11 -pipe -ffor-scope -pedantic -Wno-long-long -Werror -O0 -ffloat-store -DPTR_MODERN -DPTR_STD "
        + include_directories
    )
    return compiler, generic_command


# to be executed in the rosetta_source/src directory
# follow this command with 1) the name of the output (.cpp) file to be generated and 2) the name of the input .cxxtest.hh file
def cxxtest_testgen_command():
    return "../external/cxxtest/cxxtestgen.py --have-std --part -o "


def cxxtest_gcc_compile_command():
    return ( "g++ -c -std=c++0x -ffor-scope "
        + ' '.join( "-isystem ../"+d for d in EXTERNAL_INCLUDE_DIRS )
        + " -isystem ../external/cxxtest/ -pipe -Wall -Wextra -pedantic -Werror -Wno-long-long -Wno-strict-aliasing -march=core2 -mtune=generic -O0 -g -ggdb -ffloat-store -DBOOST_ERROR_CODE_HEADER_ONLY -DBOOST_SYSTEM_NO_DEPRECATED -DBOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS -DPTR_STD -D_GLIBCXX_DEBUG -Iexternal/cxxtest -I../. -I. "
        + ' '.join( "-I../"+d for d in EXTERNAL_INCLUDE_DIRS )
        + " -Iplatform/linux/64/gcc/7 -I../test -Iplatform/linux/64/gcc -Iplatform/linux/64 -Iplatform/linux -I/usr/include -I/usr/local/include -o "
    )


def cxxtest_test_compile(cxx_hh, verbose=False, id=""):
    out_log = "out.log"
    err_log = "err.log"
    temp_o = "temp.o"
    if id != "":
        out_log = out_log + "." + str(id)
        err_log = err_log + "." + str(id)

    errfile = open(out_log, "w")
    logfile = open(err_log, "w")
    first_compile_command = (
        cxxtest_testgen_command() + "cxx1_tmp" + str(id) + ".cpp " + cxx_hh
    )
    # print("first compile command", first_compile_command)
    job1 = subprocess.Popen(
        no_empty_args(first_compile_command.split(" ")), stderr=errfile, stdout=logfile
    )
    _, _2 = job1.communicate()
    if job1.returncode == 0:
        second_compile_command = (
            cxxtest_gcc_compile_command()
            + "cxx2_tmp"
            + str(id)
            + ".o cxx1_tmp"
            + str(id)
            + ".cpp"
        )
        # print("second compile command", second_compile_command)
        job2 = subprocess.Popen(
            no_empty_args(second_compile_command.split(" ")),
            stderr=errfile,
            stdout=logfile,
        )
        _, _2 = job2.communicate()
        return job2.returncode == 0
    return False


def test_compile(cc_file, verbose=False, id="", devnull=False, silent=False):

    compiler, generic_command = central_compile_command()

    out_log = "out.log"
    err_log = "err.log"
    temp_o = "temp.o"
    if devnull:
        temp_o = "/dev/null"
    if id != "":
        out_log = out_log + "." + str(id)
        err_log = err_log + "." + str(id)
        temp_o = temp_o + "." + str(id)

    command = compiler + " -o " + temp_o + generic_command + " " + cc_file
    command_list = command.split(" ")
    errfile = open(out_log, "w") if id else sys.stderr
    logfile = open(err_log, "w") if id else sys.stdout

    command_list = no_empty_args(command_list)

    if verbose:
        print(command)

    return_code = subprocess.call(command_list, stderr=errfile, stdout=logfile)

    if verbose:
        print("return code: ", return_code, type(return_code))

    if id:
        errfile.close()
        logfile.close()

    if return_code == 0:
        return True
    else:
        # print file(out_log).read(), file(err_log).read()
        if not silent:
            # print( command )
            print(
                "To compile this header locally run following command: " +
                "cd source/src && python ../../tools/python_cc_reader/" +
                "test_all_headers_compile_w_fork.py --headers",
                cc_file,
                "\n\n",
            )
        return False

    # return return_code == 0


def test_compile_from_lines(filelines, verbose=False):
    # print("file lines:")
    # print("".join(filelines))
    # print()
    # print("----")
    # print()
    compiler, generic_command = central_compile_command()

    command = compiler + " -o /dev/null " + generic_command + " -x c++ -"
    command_list = command.split(" ")
    command_list = no_empty_args(command_list)

    # outfile = open("test_compile.log", "w")
    # errfile = open("test_compile.err", "w")
    # with open("example_input.cc", "w") as fid:
    #     fid.writelines(filelines)

    job = subprocess.Popen(
        command_list,
        stderr=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stdin=subprocess.PIPE,
        encoding='utf8'
    )
    job.communicate("".join(filelines))
    job.wait()

    if job.returncode == 0:
        return True
    else:
        if verbose:
            print("test_compile_from_lines return code:", job.returncode)
            print(command)
            open("blah", "w").writelines(filelines)
            print(job.stderr)
        return False


# C_cc may be either a header or a cc file; either will compile
def test_compile_from_stdin(C_cc, file_contents):

    if C_cc not in file_contents:
        print(C_cc, "file not found in file_contents")
        return False

    return test_compile_from_lines(expand_includes_for_file(C_cc, file_contents))


def generate_objdump_for_file(fname, id=""):
    compiler, generic_command = central_compile_command()
    temp_o = "temp.o"
    temp_objdump = "temp.objdump"
    if id != "":
        temp_o = temp_o + "." + str(id)
    command = " ".join([compiler,"-o", temp_o, generic_command, fname])
    # print(command)
    command_list = no_empty_args(command.split(" "))

    job = subprocess.Popen(
        command_list,
        stderr=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stdin=subprocess.PIPE
    )
    out, err = job.communicate()

    if job.returncode == 0:
        if len(fname) > 3 and fname[-3:] == ".cc":
            command2 = " ".join(["objdump -d", temp_o])
            command_list2 = no_empty_args(command2.split(" "))
            job2 = subprocess.Popen(
                command_list2,
                stderr=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stdin=subprocess.PIPE,
                encoding='utf8'
            )
            out, err = job2.communicate()
            #print("out", out)
            #print("err", err)
            objdump = relabel_sections(out.splitlines(True))
            return True, objdump
        else:
            return True, []

    return False, None


def generate_objdump(filelines, id=""):
    compiler, generic_command = central_compile_command()
    temp_o = "temp.o"
    temp_objdump = "temp.objdump"
    if id != "":
        temp_o = temp_o + "." + str(id)

    command = compiler + " -o " + temp_o + " " + generic_command + " -x c++ -"
    # print command
    command_list = command.split(" ")
    command_list = no_empty_args(command_list)

    job = subprocess.Popen(
        command_list,
        stderr=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stdin=subprocess.PIPE,
    )
    job.communicate("".join(filelines))
    job.wait()

    if job.returncode == 0:
        # OK -- it compiles -- but does it compile the same way?
        command2 = "objdump -d " + temp_o
        # print command2
        command_list2 = command2.split(" ")
        command_list2 = no_empty_args(command_list2)
        job2 = subprocess.Popen(
            command_list2,
            stderr=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stdin=subprocess.PIPE,
        )
        out, err = job2.communicate()
        objdump = relabel_sections(out.splitlines(True))
        return True, objdump

    return False, None


def test_compile_extreme(filelines, gold_objdump, id=""):
    builds, test_objdump = generate_objdump(filelines, id)
    # if ( test_objdump ) : open( "test_objdump.objdump", "w" ).writelines( test_objdump ) #temp debug
    if not builds:
        # print "test compile extreme: build fails" #temp debug
        # open( "failed_filelines.txt","w").writelines( filelines );
        return False
    else:
        if compare_objdump_lines(gold_objdump, test_objdump):
            return True
        else:
            # print "test compile extreme: objdump comparison fails"
            return False

def test_compile_for_file_extreme(fname, gold_objdump, id=""):
    builds, test_objdump = generate_objdump_for_file(fname, id)
    # if ( test_objdump ) : open( "test_objdump.objdump", "w" ).writelines( test_objdump ) #temp debug
    if not builds:
        # print "test compile extreme: build fails" #temp debug
        # open( "failed_filelines.txt","w").writelines( filelines );
        return False
    else:
        if compare_objdump_lines(gold_objdump, test_objdump):
            return True
        else:
            # print "test compile extreme: objdump comparison fails"
            return False

def test_compile_w_surrogates(fname, surrogates, id=""):
    compiles = test_compile(fname, id=id, silent=True)
    if compiles:
        for surrogate in surrogates:
            # print("testing compilation of surrogate", surrogate)
            if not test_compile(surrogate, id=id, silent=True):
                return False
        return True
    else:
        return False


def tar_everything(tar_file_name):
    dirs_to_tar = ["core", "devel", "apps", "protocols"]
    if len(tar_file_name) < 8 or tar_file_name[len(tar_file_name) - 7 :] != ".tar.gz":
        tar_file_name = tar_file_name + ".tar.gz"
    command = "tar -czf " + tar_file_name
    command_list = command.split(" ")
    command_list = no_empty_args(command_list)
    command_list.extend(dirs_to_tar)
    subprocess.call(command_list)


def tar_together_files(tar_file_name, filelist):
    if len(tar_file_name) < 8 or tar_file_name[len(tar_file_name) - 7 :] != ".tar.gz":
        tar_file_name = tar_file_name + ".tar.gz"
    command = "tar -czf " + tar_file_name
    command_list = command.split(" ")
    command_list.extend(filelist)
    subprocess.call(command_list)
