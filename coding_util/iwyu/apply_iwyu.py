#!/usr/bin/env python

'''apply_iwyu.py - a script to read the output of run_iwyu and apply the listed changed.

Note - this assumes that no changes to the relevant file have happened between the run_iwyu.py
and the current entry. If that's not the case, you could get some odd results.

Run from the Rosetta/main/source directory and provide the script with files or directories,
and the .riwyuf therein will be applied.
'''

from __future__ import print_function

import sys, os
import subprocess
import codecs
import json

from optparse import OptionParser

########## Internal config #############################

#These are the clang commandline flags for debug mode, stripped of warning issues
commandline_flags_linux = '''-c -std=c++11 -isystem external/boost_submod/ -isystem external/ -isystem external/include/ -isystem external/dbio/ -isystem external/rdkit -isystem external/libxml2/include -isystem external/cxxtest/ -pipe -Qunused-arguments -DUNUSUAL_ALLOCATOR_DECLARATION -ftemplate-depth-256 -stdlib=libstdc++ -Wno-long-long -Wno-strict-aliasing -O0 -g -DBOOST_ERROR_CODE_HEADER_ONLY -DBOOST_SYSTEM_NO_DEPRECATED -DBOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS -DPTR_STD -Isrc -I./ -Itest/ -Isrc/platform/linux -ferror-limit=1 -DIWYU_SCAN'''.split()

commandline_flags = commandline_flags_linux

####################################

DEBUG = False

def find_insertion_position(contents):
    # Now we need to find the insertion line.
    # For now, we heuristically put it right after the last preprocessor directive,
    # but before any non-comment/non-preprocessor line.
    # We also preferentially put it with any previously added auto lines.
    insertion_pos = 0
    if_loop_level = 0
    last_auto_include = 0
    for ii, line in enumerate(contents):
        line = line.strip()
        if len(line) == 0:
            continue
        if line.startswith("//"):
            continue
        if line.startswith("/*") or line.startswith("*/") or line.startswith("*"):
            continue # C-style comment blocks
        if line.startswith("#"):
            if line.startswith('#if') and not line.startswith('#ifndef INCLUDED'):
                if_loop_level += 1
                continue
            if line.startswith('#endif'):
                if_loop_level -= 1
            if line.startswith("#include") and "AUTO IWYU" in line:
                last_auto_include = ii+1
            if line.startswith("#include") and if_loop_level == 0:
                insertion_pos = ii+1
        elif if_loop_level == 0 and insertion_pos != 0 :
            break
        elif line.startswith("namespace"):
            break

    if last_auto_include != 0:
        return last_auto_include
    else:
        return insertion_pos

LIB_LEVELS = {
    "apps":0,
    "devel":1,
    "protocols":2,
    "core":3,
    "basic":4,
    "numeric":5,
    "utility":6,
    "ObjexxFCL":7,
    "cifparse":8,
    "rdkit":9,
    "libxml2":10,
    "libzmq":11,
    "cppdb":12,
    "sqlite3":13,
}

def include_sort_key(entry):
    parts = entry.split('/')
    return ( LIB_LEVELS.get(parts[0], 99), entry )

class CompileTest:
    def __init__(self, filename, instructions, options):
        self.filename = filename
        with open(filename) as f:
            self.contents = f.readlines()
        self.modified = False
        self.failed = False

        self.deletions = instructions['deletions'] # { headername: [lineno] }
        self.possible_additions = [] # [ (sortorder, line) ]
        self.known_additions = [] # [ (sortorder, line) ]

        just_includes = [ l for l in self.contents if l.startswith("#include") ]
        for fn in instructions['additions']:
            whys = instructions['additions'][fn]
            newline = '#include <'+fn+'> // AUTO IWYU For ' + ' '.join(whys) + '\n'
            if newline not in just_includes:
                self.possible_additions.append( ( include_sort_key(fn), newline ) )
            else:
                if DEBUG: print("Ignoring addition of",fn,"as we already included it.")
        self.possible_additions.sort()


        if len(self.possible_additions) != 0:
            self.insert_pos = find_insertion_position( self.contents )
            if self.insert_pos == 0:
                print("ERROR: Can't find insertion position --", filename)
                self.insert_pos = None
        else:
            self.insert_pos = None


        if options.compiler is None:
            self.compiler = "clang++"
        else:
            self.compiler = options.compiler
        self.testfilename = filename[:-3] + ".test" + filename[-3:]
        if options.overwrite:
            self.outfilename = filename
        else:
            self.outfilename = filename[:-3] + ".trim" + filename[-3:]

    def test_compile(self, contents):
        with open(self.testfilename, 'w') as f:
            f.writelines(contents)

        command = [ self.compiler ] + commandline_flags + [ self.testfilename, '-o', self.testfilename +'.o' ]

        run = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = run.communicate()
        #stdout = codecs.decode( stdout, "UTF-8", "replace")
        #stderr = codecs.decode( stderr, "UTF-8", "replace")

        if os.path.exists( self.testfilename+'.o' ):
            os.remove( self.testfilename+'.o' )

        return run.returncode == 0

    def inserted_contents(self, contents, insertion_pos, additions=None):
        if insertion_pos is None or insertion_pos == 0:
            return contents

        new_contents = contents[:insertion_pos]
        if contents[insertion_pos-1] != '\n' and "AUTO IWYU" not in contents[insertion_pos-1]:
            new_contents.append('\n')
        new_contents += [l for (o,l) in self.known_additions]
        if additions is not None:
            new_contents += [l for (o,l) in additions]
        if contents[insertion_pos] != '\n' and "AUTO IWYU" not in contents[insertion_pos]:
            new_contents.append('\n')
        new_contents += contents[insertion_pos:]

        return new_contents

    def test_additions(self, contents,insertion_pos):

        # If it works as is, we can simply skip the additions.
        if self.test_compile( self.inserted_contents(contents,insertion_pos) ):
            return True

        if len(self.possible_additions) == 0:
            if DEBUG: print("File does not compile, and there's no remaining additions we can make.")
            return False # Skip additional trials, we're at the maximum already.

        if not self.test_compile(self.inserted_contents( contents, insertion_pos, self.possible_additions )):
            if DEBUG: print("File does not compile, even with all remaining additions")
            return False # We can't even compile with all the contents added, don't bother minimizing

        possible_additions = self.possible_additions
        unneeded_additions = []

        # We know it works with all of them.
        # Now try to remove the added ones one-by-one to come up with a minimal set.
        # possible_additions should be in a most-encompasing to most-specific order.
        progress = True
        while progress:
            progress = False
            for ii in range(len(possible_additions)):
                trial_additions = possible_additions[:ii] + possible_additions[ii+1:]
                if self.test_compile(self.inserted_contents( contents, insertion_pos, trial_additions)):
                    # Header ii is unneeded, so we can remove it and continue working on the reduced set.
                    if DEBUG: print("Header file", possible_additions[ii][1], "is not needed to compile - ignoring.")
                    unneeded_additions.append( possible_additions[ii] )
                    possible_additions = trial_additions
                    progress = True
                    break

        # possible_additions should now be whittled down to the minimal set.
        self.known_additions += possible_additions
        self.possible_additions = unneeded_additions
        self.modified = True

        return True

    def test_modifications(self):

        # First check to see if we're not compiling, but can if we add suggested additions.
        if not self.test_additions(self.contents,self.insert_pos):
            # We failed - no sense of testing any of the deletions, as we can't recover:
            self.failed = True
            if DEBUG: print("File", self.filename, "doesn't even compile without removals.")
            return

        # We go through in descending library level order, to try to keep the lower level headers around in preference to the upper ones.
        for fn in sorted( self.deletions.keys(), key=include_sort_key ):
            test_contents = self.contents[:] # Make a copy
            modified = False
            for lineno in self.deletions[fn]:
                lineno = int(lineno) -1  # Compensate for zero-based
                line = test_contents[lineno]
                splitline = line.split()
                if line.startswith("// AUTOREMOVED IWYU:"):
                    continue # Don't warn about re-removing line.
                if len(splitline) < 2:
                    print("ERROR: Line to remove does not match expected content --", self.filename)
                    continue
                if splitline[0] != '#include':
                    print("ERROR - Line to remove does not match expected content --", self.filename)
                    continue
                if splitline[1][1:-1] != fn:
                    print("ERROR:: Line to remove does not match expected content --", self.filename)
                    continue
                if "DO NOT AUTO-REMOVE" in line:
                    if DEBUG:
                        print("Skipping line removal for ", fn, " - comment tells it to stay --", self.filename)
                    continue
                test_contents[lineno] = '// AUTOREMOVED IWYU: ' + line.lstrip()
                modified = True
            if not modified:
                continue
            if DEBUG: print("Trying to remove",fn,"from",self.filename)
            if self.test_additions(test_contents,self.insert_pos):
                if DEBUG: print("... removal succeeded")
                self.contents = test_contents
                self.modified = True
            else:
                print("\tIssue removing",fn,"from",self.filename)

    def write_results(self):
        if os.path.isfile(self.testfilename):
            os.remove(self.testfilename)
        if self.modified and not self.failed:
            self.known_additions.sort() # Make sure we have a decent order to the included headers.
            with open(self.outfilename, 'w') as f:
                f.writelines(self.inserted_contents(self.contents,self.insert_pos))
            return True
        else:
            if DEBUG: print('\t', "No alterations to",self.filename,"were made")
            return False


def apply_changes(filename, instructions, options ):
    '''Apply the changes for the filename.
    Return true on success and false on failure.'''

    obj = CompileTest(filename,instructions,options)

    obj.test_modifications()
    return obj.write_results()

def process_file(filename, options):

    print("PROCESSING", filename)
    with open(filename) as f:
        instructions = json.load(f)

    success = apply_changes( filename[:-len('.riwyuf')], instructions, options )

    if success and not options.nodelete:
        os.remove( filename )

def process_dir(dirname):
    for item in os.listdir(dirname):
        name = os.path.join(dirname,item)
        if os.path.isdir(name):
            for fn in process_dir(name):
                yield fn
        elif os.path.isfile(name):
            if name.endswith(".riwyuf"):
                yield name

def find_files(pathlist):
    for name in pathlist:
        if os.path.isdir(name):
            for fn in process_dir(name):
                yield fn
        elif os.path.isfile(name) and name.endswith('.riwyuf'):
            yield name
        elif os.path.isfile( name+'.riwyuf' ):
            yield name+'.riwyuf'
        else:
            print("Cannot find file or directory: "+name)
            if os.path.isfile(name):
                print("Corresponding .riwyuf not found")

if __name__ == "__main__":
    if not os.path.exists( "./src" ):
        print( "Script must be run from within Rosetta/main/source/" )
        exit()

    parser = OptionParser(usage="usage: %prog [OPTIONS] [FILES|DIRECTORIES]")
    parser.set_description(__doc__)
    parser.add_option("--nodelete",
      default=False, action="store_true",
      help="Don't delete the .riwyuf file afterwards." )
    parser.add_option("--debug",
      default=False, action="store_true",
      help="Print extra debugging info." )
    parser.add_option("--overwrite",
      default=False, action="store_true",
      help="Overwrite the main file, rather than making a new trim file." )
    parser.add_option("--compiler",
      default=None,
      help="The command for the Clang compiler to use in compilation testing." )
    parser.add_option("-j",
      default=0, type=int,
      help="Number of processes to run." )

    (options, args) = parser.parse_args(args=sys.argv[1:])

    DEBUG = options.debug

    if options.j == 0:
        for fn in find_files(args):
            process_file(fn, options)
    else:
        import multiprocessing
        pool = multiprocessing.Pool(options.j)

        def callback(arg):
            print("Processing is done")

        for fn in find_files(args):
            pool.apply_async(process_file, (fn, options) )

        pool.close()
        pool.join()
