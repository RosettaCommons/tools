import sys, signal
import willclang
import test_look_for_inherritance as tlfi
from clang.cindex import Index,CursorKind,TypeKind;

import os, commands, re, subprocess, time
from os import path
from optparse import OptionParser

sys.path.insert( 0, os.path.realpath(__file__).rpartition("/")[0]+"/../python_cc_reader" )
from python_cc_reader.cpp_parser import code_utilities


class OI:
    def __init__(self, **entries): self.__dict__.update(entries)


class CodeQualityScanner:
    def __init__(self):
        self.jobs = []    # list of spawned process pid's
        self.output = ''  # output of current job
        self.exclude = None
        self.code_checkers = []

    def set_excludes( self, excludes ) :
        self.exclude = excludes

    def log(self, message):
        self.output += message
        if not Options.quiet: print message

    def mfork(self):
        '''
        Check if number of child process is below Options.jobs. And if it is - fork the new pocees and return its pid
        and a pipe that the child process can write to.
        '''
        while len(self.jobs) >= Options.jobs :
            for p,readpipe in self.jobs[:] :
                try :
                    r = os.waitpid(p, os.WNOHANG)
                except OSError :
                    print "pid", p, "was not found, assume it exited with status 0?"
                    sys.stdout.flush()
                    r = (p, 0)
                if r == (p, 0):  # process have ended without error
                    newtxt = readpipe.read()
                    if newtxt :
                        self.output += newtxt
                    self.jobs.remove( (p,readpipe) )
                elif r[0] == p :  # process ended but with error, special case we will have to wait for all process to terminate and call system exit.
                    self.jobs.remove( (p,readpipe) )
                    for p,readpipe in self.jobs:
                        try :
                            os.waitpid(p, 0)
                        except OSError :
                            print "pid", p, "was not found; assuming it exited with status 0."
                    print 'One of the jobs exited with non-zero status'
                    sys.exit(1)

            if len(self.jobs) >= Options.jobs: time.sleep(.5)
        readpipe, writepipe = os.pipe()
        pid = os.fork()
        if pid:
            # we are parent!
            os.close(writepipe)
            readpipe = os.fdopen(readpipe)
            #print "forked", pid
            sys.stdout.flush()
            self.jobs.append( (pid,readpipe) )
        else :
            # we are the child!
            os.close( readpipe )
            writepipe = os.fdopen( writepipe, "w" )
        return ( pid, writepipe )

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
        all_library_files = code_utilities.find_library_files()
        #all_library_files = [ "protocols/unfolded_state_energy_calculator/UnfoldedStateEnergyCalculatorMover.cc" ]
        #all_library_files = ["utility/sql_database/DatabaseSessionManager.hh"]
        #all_library_files = [
        #    "protocols/init/init.MoverRegistrators.ihh",
        #    "protocols/backrub/BackrubSidechainMoverCreator.hh",
        #]

        # discard injection headers
        all_library_files = filter( lambda x : len(x) < 4 or ( x[-4:] != ".ipp" and x[-4:] != ".ihh" ), all_library_files )

        if self.exclude :
            before_len = len( all_library_files )

            # discard any files that were requested to be skipped on the command line
            all_library_files = filter( lambda x : x not in self.exclude, all_library_files )

            after_len = len( all_library_files )
            print "Excluding", before_len - after_len, "files from the source tree"

        signal.signal(signal.SIGINT, self.signal_handler)
        count = 0
        nfiles = len(all_library_files)
        for srcfile in all_library_files :
            count += 1
            print "Examining file:", srcfile, "#", count, "of", nfiles, "(", "%.1f" % ( 100 * float(count) / nfiles), "% )"
            sys.stdout.flush()
            pid,writepipe = self.mfork()
            if not pid:  # we are child process
                self.examine_ast_for_file( srcfile, writepipe )
                sys.exit(0)

        for (p,readpipe) in self.jobs:
            try :
                os.waitpid(p, 0)  # waiting for all child process to termintate...
                newtxt = readpipe.read()
                if newtxt :
                    self.output += newtxt
            except OSError :
                print "problem waiting for process", p, "to exit -- apparently no such process"
        if self.output :
            print "Error messages:"
            print self.output
            sys.exit(1)

    def examine_ast_for_file( self, srcfile, writepipe ) :
        fn = os.path.abspath(srcfile);
        srcdir=os.path.abspath(".")
        rosetta_source_basedir=srcdir[:-3]
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

        for code_checker in self.code_checkers :
            code_checker.examine_ast( ast )
            if code_checker.has_problem() :
                writepipe.write( code_checker.problems() )

        writepipe.close()

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

    scanner = CodeQualityScanner()
    if Options.exclude_file :
        exclude_set = set( [ x.strip() for x in open( Options.exclude_file ).readlines() ] )
        scanner.set_excludes( exclude_set )
    scanner.code_checkers.append( tlfi.CodeQualityChecker_FindStackDeclaredObjectsReservedForTheHeap() )
    scanner.code_checkers.append( tlfi.CodeQualityChecker_FindRawPtrs_to_RefcountSubclasses() )
    scanner.run()



if __name__ == "__main__":
    main(sys.argv)
