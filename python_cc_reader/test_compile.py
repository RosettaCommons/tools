# This file will test the build
# You will have to modify three variables here for your own
# computer:  os, nbits, and compiler.

import subprocess, re, time, sys
from code_utilities import expand_includes_for_file, load_source_tree
from reinterpret_objdump import relabel_sections, compare_objdump_lines


def no_empty_args( command_list ) :
   clprime = []
   for arg in command_list :
      if arg != "" :
         clprime.append( arg )
   return clprime

def central_compile_command() :

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

   include_directories = " -I./ -I../external -I../external/include -Iplatform/" + os + "/" + nbits + "/gcc -Iplatform/" + os + "/" + nbits +  " -Iplatform/" + os + " -I../external/boost_1_55_0 -I/usr/local/include -I/usr/include/ -I../external/dbio -I../external/libxml2/include"

   generic_command = " -c -std=c++11 -pipe -ffor-scope -w -pedantic -Wno-long-long -O0 -ffloat-store -DPTR_MODERN -DPTR_STD" + include_directories
   return compiler, generic_command

# to be executed in the rosetta_source/src directory
# follow this command with 1) the name of the output (.cpp) file to be generated and 2) the name of the input .cxxtest.hh file
def cxxtest_testgen_command() :
   return "../external/cxxtest/cxxtestgen.py --have-std --part -o "


def cxxtest_gcc_compile_command() :
   return "g++ -c -isystem ../external/boost_1_55_0/boost/ -O0 -g -ggdb -ffloat-store -I../external/cxxtest -I../. -I../test -I../src -I../external -I../external/include -Iplatform/linux/64/gcc -Iplatform/linux/64 -Iplatform/linux -I../external/boost_1_55_0 -I../external/dbio -I/usr/local/include -I/usr/include -o "

def cxxtest_test_compile( cxx_hh, verbose=False, id="" ) :
   out_log = "out.log"
   err_log = "err.log"
   temp_o  = "temp.o"
   if id != "" :
     out_log = out_log + "." + str( id )
     err_log = err_log + "." + str( id )

   errfile = open( out_log, "w" )
   logfile = open( err_log, "w" )
   first_compile_command = cxxtest_testgen_command() + "cxx1_tmp" + str(id) + ".cpp " + cxx_hh
   return_code = subprocess.call( no_empty_args( first_compile_command.split(" ")), stderr=errfile, stdout=logfile )
   if return_code ==  0:
      second_compile_command = cxxtest_gcc_compile_command() + "cxx2_tmp" + str(id) + ".o cxx1_tmp" + str(id) + ".cpp"
      return_code = subprocess.call( no_empty_args( second_compile_command.split(" ")), stderr=errfile, stdout=logfile )
      return return_code == 0
   return False


def test_compile( cc_file, verbose=False, id="", devnull=False ) :

   compiler, generic_command = central_compile_command()

   out_log = "out.log"
   err_log = "err.log"
   temp_o = "temp.o"
   if devnull : temp_o  = "/dev/null"
   if id != "" :
     out_log = out_log + "." + str( id )
     err_log = err_log + "." + str( id )
     temp_o  = temp_o  + "." + str( id )

   command = compiler + " -o " + temp_o + generic_command + " " + cc_file
   command_list = command.split(" ")
   errfile = open( out_log, "w" ) if id else sys.stderr
   logfile = open( err_log, "w" ) if id else sys.stdout

   command_list = no_empty_args( command_list )

   if (verbose) :
      print command

   return_code = subprocess.call( command_list, stderr=errfile, stdout=logfile )

   if (verbose) :
      print "return code: ", return_code

   if id: errfile.close();  logfile.close()

   if return_code == 0:
      return True
   else:
      #print file(out_log).read(), file(err_log).read()
      print 'To compile this header locally run following command: cd source/src && python ./../../../tools/python_cc_reader/test_all_headers_compile_w_fork.py --headers', cc_file, '\n\n'
      return False

   #return return_code == 0

def test_compile_from_lines( filelines, verbose=False ) :
   compiler, generic_command = central_compile_command()

   command = compiler + " -o /dev/null " + generic_command + " -x c++ -"
   command_list = command.split(" ")
   command_list = no_empty_args( command_list )

   job = subprocess.Popen( command_list, stderr=subprocess.PIPE, stdout=subprocess.PIPE, stdin=subprocess.PIPE )
   job.communicate( "".join( filelines ))
   job.wait()

   if ( job.returncode == 0 ) :
      return True
   else :
      if verbose :
         print "test_compile_from_lines return code:", job.returncode
         print command
         open("blah","w").writelines(filelines)
         print job.stderr
      return False

# C_cc may be either a header or a cc file; either will compile
def test_compile_from_stdin( C_cc, file_contents ) :

   if not file_contents.has_key( C_cc ) :
      print C_cc, "file not found in file_contents"
      return False

   return test_compile_from_lines( expand_includes_for_file( C_cc, file_contents ) )

def generate_objdump( filelines, id = "" ) :
   compiler, generic_command = central_compile_command()
   temp_o  = "temp.o"
   temp_objdump = "temp.objdump"
   if id != "" :
     temp_o  = temp_o  + "." + str( id )

   command = compiler + " -o " + temp_o + " " + generic_command + " -x c++ -"
   #print command
   command_list = command.split(" ")
   command_list = no_empty_args( command_list )

   job = subprocess.Popen( command_list, stderr=subprocess.PIPE, stdout=subprocess.PIPE, stdin=subprocess.PIPE )
   job.communicate( "".join( filelines ))
   job.wait()

   if ( job.returncode == 0 ) :
      # OK -- it compiles -- but does it compile the same way?
      command2 = "objdump -d " + temp_o
      #print command2
      command_list2 = command2.split(" ")
      command_list2 = no_empty_args( command_list2 )
      job2 = subprocess.Popen( command_list2, stderr=subprocess.PIPE, stdout=subprocess.PIPE, stdin=subprocess.PIPE )
      out, err = job2.communicate()
      objdump = relabel_sections( out.splitlines(True) )
      return True, objdump

   return False, None


def test_compile_extreme( filelines, gold_objdump, id = "" ) :
   builds, test_objdump = generate_objdump( filelines, id )
   #if ( test_objdump ) : open( "test_objdump.objdump", "w" ).writelines( test_objdump ) #temp debug
   if not builds :
      #print "test compile extreme: build fails" #temp debug
      #open( "failed_filelines.txt","w").writelines( filelines );
      return False
   else :
      if compare_objdump_lines( gold_objdump, test_objdump ) :
         return True
      else :
         #print "test compile extreme: objdump comparison fails"
         return False

def tar_everything( tar_file_name ) :
   dirs_to_tar = ["core", "devel", "apps", "protocols" ]
   if len( tar_file_name ) < 8 or tar_file_name[ len(tar_file_name)-7:] != ".tar.gz" :
      tar_file_name = tar_file_name + ".tar.gz"
   command = "tar -czf " + tar_file_name
   command_list = command.split(" ")
   command_list = no_empty_args( command_list )
   command_list.extend( dirs_to_tar )
   subprocess.call( command_list )

def tar_together_files( tar_file_name, filelist ) :
   if len( tar_file_name ) < 8 or tar_file_name[ len(tar_file_name)-7:] != ".tar.gz" :
      tar_file_name = tar_file_name + ".tar.gz"
   command = "tar -czf " + tar_file_name
   command_list = command.split(" " )
   command_list.extend( filelist )
   subprocess.call( command_list )


if __name__ == "__main__" :
   if len(sys.argv) < 2 :
      print "Usage: python test_compile.py <filename>"
      sys.exit(1)
   print "First testing compilation directly from .cc file"
   compiled = test_compile( sys.argv[1], True )

   print "Now testing compilation using python-expanded #includes"
   compilable_files, all_includes, file_contents = load_source_tree()
   print "...source tree loaded"
   if sys.argv[1] not in file_contents :
      print "File", sys.argv[1], "not found in source tree"
      sys.exit(1)
   compiled = test_compile_from_lines( expand_includes_for_file( sys.argv[1], file_contents), verbose=True )
   if ( compiled ) :
      print "Success"
   else :
      test_compile( sys.argv[1], True )
      print "Failed"
