import sys, signal
import willclang
import willclang.test_look_for_inherritance as tlfi
from clang.cindex import Index,CursorKind,TypeKind;

import os, commands, re, subprocess, time
from os import path
from optparse import OptionParser

sys.path.insert( 0, os.path.realpath(__file__).rpartition("/")[0]+"/../python_cc_reader" )

from python_cc_reader.cpp_parser import code_utilities


class OI:
    def __init__(self, **entries): self.__dict__.update(entries)


class Runner:
    def __init__(self):
        self.jobs = []    # list of spawned process pid's
        self.output = ''  # output of current job
        self.exclude = None

    def set_excludes( self, excludes ) :
        self.exclude = excludes

    def log(self, message):
        self.output += message
        if not Options.quiet: print message

    def mfork(self):
        ''' Check if number of child process is below Options.jobs. And if it is - fork the new pocees and return its pid.
        '''
        while len(self.jobs) >= Options.jobs :
            for p in self.jobs[:] :
                try :
                    r = os.waitpid(p, os.WNOHANG)
                except os.OSError :
                    print "pid", p, "was not found, assume it exited with status 0?"
                    sys.stdout.flush()
                    r = (p, 0)
                if r == (p, 0):  # process have ended without error
                    self.jobs.remove(p)
                elif r[0] == p :  # process ended but with error, special case we will have to wait for all process to terminate and call system exit.
                    for p in self.jobs: os.waitpid(p, 0)
                    print 'Some of the unit test suite terminate abnormally!'
                    sys.exit(1)

            if len(self.jobs) >= Options.jobs: time.sleep(.5)
        pid = os.fork()
        if pid: self.jobs.append(pid) # We are parent!
        return pid


    def signal_handler(self, signal_, f):
        print 'Ctrl-C pressed... killing child jobs...'
        for pid in self.jobs:
            os.killpg(os.getpgid(pid), signal.SIGKILL)


    def runCommandLine(self, file_, line_, command_line):
        self.log("Running %s:%s %s" % (file_, line_, command_line) )

        f = subprocess.Popen(command_line, bufsize=0, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).stderr
        for line in f:
            self.log(line)
            sys.stdout.flush()
        f.close()


    def run(self):
        print "Loading source tree"
        compilable_files, all_includes, file_contents = code_utilities.load_source_tree()
        all_rosetta_files_compiled = file_contents.keys() # include .cc files and .hh files
        #all_rosetta_files_compiled = ["utility/sql_database/DatabaseSessionManager.hh"]
        #all_rosetta_files_compiled = [
        #    "protocols/init/init.MoverRegistrators.ihh",
        #    "protocols/backrub/BackrubSidechainMoverCreator.hh",
        #]
        if self.exclude :
            before_len = len( all_rosetta_files_compiled )
            all_rosetta_files_compiled = filter( lambda x : x not in self.exclude, all_rosetta_files_compiled )
            after_len = len( all_rosetta_files_compiled )
            print "Excluding", before_len - after_len, "files from the source tree"

        file_contents = None; compilable_files = None; all_includes = None; # trash all this
        signal.signal(signal.SIGINT, self.signal_handler)
        count = 0
        nfiles = len(all_rosetta_files_compiled)
        for srcfile in all_rosetta_files_compiled :
            count += 1
            if srcfile[-4:] == ".ihh" : continue
            if srcfile[-4:] == ".ipp" : continue
            print "Examining file:", srcfile, "#", count, "of", nfiles, "(", "%.1f" % ( 100 * float(count) / nfiles), "% )"
            sys.stdout.flush()
            pid = self.mfork()
            if not pid:  # we are child process
                fn = os.path.abspath(srcfile);
                srcdir=os.path.abspath(".")
                rosetta_source_basedir="/home/andrew/scr2/GIT/rosetta/rosetta_source/"
                clangargs=[
                    "-I%s"%(srcdir),
                    "-I%s"%(rosetta_source_basedir + "src/platform/macos"),
                    "-I%s"%(rosetta_source_basedir + "external"),
                    "-I%s"%(rosetta_source_basedir + "external/boost_1_46_1"),
                    "-I%s"%(rosetta_source_basedir + "external"),
                    "-I%s"%(rosetta_source_basedir + "external/include"),
                    "-I%s"%(rosetta_source_basedir + "external/dbio"),
                    "-I%s"%("/usr/lib/clang/2.9/include")
                ];
                src = willclang.SourceFile(fn,clangargs);
                try :
                    ast = src.get_ast()
                except willclang.ASTException as e :
                    print "problem with file", srcfile
                    print e
                    sys.exit(0)
                cr = tlfi.CodeReader(ast)
                ast.root.treemap( cr.find_class_declarations, allchild=True )
                cr.note_classes_deriving_from_reference_count()
                rfptrfinder = tlfi.VarDecToRefCountSubclassPointerFinder( cr )
                ast.root.treemap( rfptrfinder.save_lines_w_rawptrs_assigned_opdata, allchild=False )
                if rfptrfinder.lines_w_rawpts_assigned_opdata :
                    open( srcfile + ".ptrprobs", "w" ).writelines( rfptrfinder.lines_w_rawpts_assigned_opdata )
                sys.exit(0)

        for p in self.jobs:
            try :
                os.waitpid(p, 0)  # waiting for all child process to termintate...
            except os.OSError :
                print "problem waiting for process", p, "to exit -- apparently no such process"



def main(args):
    ''' Script to run Jobs in parallel.
    '''
    parser = OptionParser(usage="usage: %prog [OPTIONS] file_with_command_lines [file2] [file3] ...")
    parser.set_description(main.__doc__)


    parser.add_option("-j", "--jobs",
      default=1,
      type="int",
      help="Number of processors to use when running testss (default: 1)",
    )

    parser.add_option("-p", '--prefix',
      default='',
      action="store",
      help="Specify the prefix for files name where output is saved. Default is '' - which mean no output is saved.",
    )

    parser.add_option('-q', "--quiet", action="store_true", dest="quiet", default=False,
      help="Suppress (mute) output to std::out."
    )

    parser.add_option('-x', "--exclude", action="store", dest="exclude_file", default="",
      help="Give a file listing all source files that should be excluded from traversal"
    )

    (options, args) = parser.parse_args(args=args[1:])

    global Options;  Options = options

    R = Runner()
    if Options.exclude_file :
        exclude_set = set( [ x.strip() for x in open( Options.exclude_file ).readlines() ] )
        R.set_excludes( exclude_set )
    R.run()



if __name__ == "__main__":
    main(sys.argv)
