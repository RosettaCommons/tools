# This file will test the build 
# You will have to modify three variables here for your own
# computer:  os, nbits, and compiler.

import subprocess, re, time
from code_utilities import expand_includes_for_file
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

   include_directories = " -I./ -I../external/include -Iplatform/" + os + "/" + nbits + "/gcc -Iplatform/" + os + "/" + nbits +  " -Iplatform/" + os + " -I../external/boost_1_46_1 -I/usr/local/include -I/usr/include/ -I../external/dbio"

   generic_command = " -c -std=c++98 -pipe -ffor-scope -w -pedantic -Wno-long-long -O0 -ffloat-store" + include_directories
   return compiler, generic_command

def test_compile( cc_file, verbose=False, id="" ) :

   compiler, generic_command = central_compile_command()

   out_log = "out.log"
   err_log = "err.log"
   temp_o  = "temp.o"
   if id != "" :
     out_log = out_log + "." + str( id )
     err_log = err_log + "." + str( id )
     temp_o  = temp_o  + "." + str( id )
   command = compiler + " -o " + temp_o + generic_command + " " + cc_file
   command_list = command.split(" ")
   errfile = open( out_log, "w" )
   logfile = open( err_log, "w" )

   command_list = no_empty_args( command_list )

   if (verbose) :
      print command

   return_code = subprocess.call( command_list, stderr=errfile, stdout=logfile )

   if (verbose) :
      print "return code: ", return_code

   errfile.close()
   logfile.close()

   if ( return_code == 0 ) :
      return True
   else :
      return False

   #return return_code == 0

def test_compile_from_lines( filelines ) :
   compiler, generic_command = central_compile_command()

   command = compiler + " -o /dev/null " + generic_command + " -x c++ -"
   #print command
   command_list = command.split(" ")
   command_list = no_empty_args( command_list )

   job = subprocess.Popen( command_list, stderr=subprocess.PIPE, stdout=subprocess.PIPE, stdin=subprocess.PIPE )
   job.communicate( "".join( filelines ))
   job.wait()

   if ( job.returncode == 0 ) :
      return True
   else :
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
   if not builds :
      return False
   else :
      if compare_objdump_lines( gold_objdump, test_objdump ) :
         return True
      else :
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

