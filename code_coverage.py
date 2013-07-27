#!/usr/bin/env python

import sys, os, os.path
import glob
import subprocess
import math
import pickle
from optparse import OptionParser, IndentedHelpFormatter

def lowerwilsonbound(pos, n):
    # See http://www.evanmiller.org/how-not-to-sort-by-average-rating.html
    # Basically, we want 40/100 to be a "better" score than 4/10 (or even 1/2)
    # We're arbitrarily choosing a 95% confindence interval to simplify things
    if n == 0:
        return 0
    z = 1.96
    phat = 1.0*pos/n
    return (phat + z*z/(2*n) - z * math.sqrt((phat*(1-phat)+z*z/(4*n))/n))/(1+z*z/n)

class CoverageInfo:
    def __init__(self, options):
        self.files = {} # Dictionary of filename:(dictionary of lineno:ex/nonex status)
        self.userinfo = {} # Dictionary of filename:(dictionary of lineno:(user, rev) tuples)
        self.revlist = set() # Set of git revisions between options.rev and current HEAD
        self.options = options
        if options.load:
            f = open(options.load,'rb')
            try:
                self.files, self.userinfo = pickle.load(f)
            finally:
                f.close()

    def save(self, filename):
        f = open(options.save, "wb")
        try:
            pickle.dump( (self.files, self.userinfo), f )
        finally:
            f.close()

    def process(self, filename):
        parentfile = filename[:-5].replace("#","/") # get rid of ".gcov", convert "#" to "/"
        f = open(filename)
        try:
            for line in f:
                l = line.split(':')
                count = l[0].strip()
                if count == '-': # non-coding line
                    continue
                try:
                    lineno = int(l[1])
                except (IndexError, ValueError):
                    sys.stderr.write("Badly formated line in file %s: '%s'\n" % (filename, line) )
                    sys.exit()
                if not self.options.exits and ("assert" in l[2] or "utility_exit" in l[2]):
                    continue
                if l[2].startswith("/*EOF*/"):
                    # gcov occasionally puts execution counts on the end of files - ignore them
                    #(Not only are they irrelevant, they also mess up user statistics generation).
                    continue
                if count == "#####" or count == "====": #non-executed line
                    #If it was found in a previous run, keep the same status
                    if not lineno in self.files.setdefault(parentfile,{}):
                        self.files[parentfile][lineno] = False
                else:
                    self.files.setdefault(parentfile,{})[lineno] = True
        finally:
            f.close()

    def genuserinfo_svn(self):
        "Produce userinfo statistics."
        for filename in self.files:
            #print "Genuser for", filename #DEBUG
            p = subprocess.Popen(["svn", "blame", filename], stdout=subprocess.PIPE)
            svnout, err = p.communicate()
            if p.returncode != 0:
                sys.stdout.write("Skipping user statistics for "+filename+"\n")
                continue
            svnout = svnout.split('\n')
            fileuserinfo = self.userinfo.setdefault(filename,{})
            for lineno in self.files[filename]:
                rev, user = svnout[lineno-1].split()[0:2]
                fileuserinfo[lineno] = (user, rev)

    def genuserinfo(self):
        "Produce userinfo statistics."
        for filename in self.files:
            #print "Genuser for", filename #DEBUG
            # -c is tab seperated format, -l is full revision number, -w is ignore whitespace
            p = subprocess.Popen(["git", "blame", "-clw", filename], stdout=subprocess.PIPE)
            gitout, err = p.communicate()
            if p.returncode != 0:
                sys.stdout.write("Skipping user statistics for "+filename+"\n")
                continue
            gitout = gitout.split('\n')
            fileuserinfo = self.userinfo.setdefault(filename,{})
            for lineno in self.files[filename]:
                rev, user = gitout[lineno-1].split('\t')[0:2]
                user = user[1:].strip() # Remove the leading left paren and extraneous whitespace
                fileuserinfo[lineno] = (user, rev)

    def genrevlist(self):
        "Produce a set of revisions over which to operate."
        if self.options.rev is None:
            self.revlist = None
            return
        #Get a list of revisions between the given revision and the current head
        p = subprocess.Popen(["git", "rev-list", self.options.rev+"^..HEAD"], stdout=subprocess.PIPE)
        gitout, err = p.communicate()
        if p.returncode != 0:
            sys.stdout.write("Can't find "+self.options.rev+" in history of current commit - ignoring.\n")
            self.revlist = set()
        else:
            self.revlist = set(gitout.split("\n"))

    def filestats(self, username=None):
        "Produce a per-file summary of code coverage."
        totalex = totalnonex = 0
        filestats = {} # Dictionary of filename:[ex, noex] lists
        for fn in self.files:
            for line in self.files[fn]:
                if username is not None or self.options.rev is not None:
                    if fn not in self.userinfo:
                        sys.stdout.write("Skipping user/revision statistics for "+fn+"\n")
                        continue
                    user, rev = self.userinfo[fn][line]
                    if username is not None and user != username:
                        continue
                    if self.revlist and rev not in self.revlist:
                        continue
                if self.files[fn][line]:
                    #Executed
                    filestats.setdefault(fn, [0,0])[0] += 1
                else:
                    #Unexecuted
                    filestats.setdefault(fn, [0,0])[1] += 1
            if fn not in filestats:
                continue
            totalex += filestats[fn][0]
            totalnonex += filestats[fn][1]

        if totalnonex + totalex == 0:
           sys.stdout.write("Skipping user "+str(username)+" - no relevant lines found \n")
           return
        if username is None:
            outfilename = "byfile_coverage"
        else:
            outfilename = username.replace(" ","_") + "_coverage"
        if self.options.rev is not None:
            outfilename += "_r_" + str(self.options.rev) + ".txt"
        else:
            outfilename += ".txt"
        output = open(outfilename,"w")
        try:
            width = 80
            if self.options.invert:
                output.write("%-*s\tLinesExecuted\tTotalLines\tPercent\n" % (width, "Filename") )
                output.write("%-*s\t%d\t%d\t%0.2f\n" % ( width, "Total", totalex, totalnonex + totalex, 100.0 * totalnonex / (totalnonex + totalex) ) )
            else:
                output.write("%-*s\tLinesNotExecuted\tTotalLines\tPercent\n" % (width, "Filename") )
                output.write("%-*s\t%d\t%d\t%0.2f\n" % ( width, "Total", totalnonex, totalnonex + totalex, 100.0 * totalnonex / (totalnonex + totalex) ) )
            output.write( self.sort_and_format(filestats, width) )
        finally:
            output.close()

    def userstats(self):
        "Produce per-user statistics."
        userstats = {} # Dictionary of username:[ex, noex] lists
        for fn in self.userinfo:
            for ln in self.userinfo[fn]:
                user, rev = self.userinfo[fn][ln]
                if self.revlist and rev not in self.revlist:
                    continue
                if self.files[fn][ln]:
                    #Executed
                    userstats.setdefault(user,[0,0])[0] += 1
                else:
                    #Unexecuted
                    userstats.setdefault(user,[0,0])[1] += 1

        outfilename = "byuser_coverage"
        if self.options.rev is not None:
            outfilename += "_r_" + str(self.options.rev) + ".txt"
        else:
            outfilename += ".txt"
        output = open(outfilename,"w")
        try:
            width=10
            if self.options.invert:
                output.write("%-*s\tLinesExecuted\tTotalLines\tPercent\n" % (width, "User") )
            else:
                output.write("%-*s\tLinesNotExecuted\tTotalLines\tPercent\n" % (width, "User") )
            output.write( self.sort_and_format(userstats, width) )
        finally:
            output.close()

    def sort_and_format(self, indict, width=0):
        """Takes a dictionary of {item:[ex, noex]} and returns a formatted, sorted string with
        item   noex   ex+noex    percent
        on each line."""

        itemlist = []
        rsort = True
        for item, (ex, nonex) in indict.items():
            nlines = ex + nonex
            if nlines == 0:
                continue
            frac = float(nonex) / nlines
            if self.options.alphabetical and self.options.invert:
                itemlist.append( (item, float(ex)/nlines, ex, nlines, item ) )
                rsort = False
            elif self.options.alphabetical:
                itemlist.append( (item, float(nonex)/nlines, nonex, nlines, item ) )
                rsort = False
            elif self.options.bylines and self.options.invert:
                itemlist.append( (ex, float(ex)/nlines, ex, nlines, item ) )
            elif self.options.bylines:
                itemlist.append( (nonex, frac, float(nonex)/nlines, nlines, item ) )
            elif self.options.invert:
                itemlist.append( (lowerwilsonbound(ex, nlines), float(ex)/nlines, ex, nlines, item) )
            else:
                itemlist.append( (lowerwilsonbound(nonex, nlines), float(nonex)/nlines, nonex, nlines, item) )

        itemlist.sort(reverse = rsort)

        output = []
        for sort, frac, num, total, item in itemlist:
            output.append("%-*s\t%d\t%d\t%0.2f" % ( width, item, num, total, 100*frac) )
        return '\n'.join(output)

def findbuilddir():
    #Find the appropriate build directory
    build_ds = [ d for (d, dn, fn) in os.walk("build/src") if os.path.basename(d) == "gcov" ]
    if len(build_ds) != 1:
        sys.stderr.write("Found %d coverage build directories - can currently only handle one.\n" % len(build_ds) )
        sys.exit()
    return build_ds[0]

def main(options):

    coverageinfo = CoverageInfo(options)

    # While multiple invocations of coverage-augmented programs sum the statistics, multiple invocations of the gcov program itself do not.
    # Unfortunately, older gcov implementations (such as the one on the Baker lab digs) have limited support for out-of-tree compiling, and as such must
    # be processed directory by directory, which necessitates collecting the statistics for each directory before rerunning, as subsequent invocations will clobber them

    ### Non-workable "all at once" run
    ##ccfiles = [ os.path.join(d,f) for (d, dn, fn) in os.walk("src") for f in fn if os.path.splitext(f)[1] == ".cc"]
    ##subprocess.call(["gcov", "-p", "-o", build_d] + ccfiles)

    if options.file:
        if not options.file.startswith("src"):
            sys.stderr.write("Filename to process must be a relative path starting with src/ \n")
            sys.exit()

        #Check that we aren't going to clobber existing gcov files
        if len(glob.glob("*.gcov")):
            sys.stderr.write("There already exists *.gcov files in this directory. Please delete them first. \nThis script doesn't play well with others ;).\n")
            sys.exit()

        build_d = findbuilddir()

        stem = os.path.dirname(options.file)[4:] # remove 'src/'
        gcovcall = ["gcov", "-p", "-o", os.path.join(build_d,stem), options.file]
        subprocess.call(gcovcall)
        sys.stdout.write("\n\nCalled: " + ' '.join(gcovcall) + "\n\n" )
        sys.exit()

    if not coverageinfo.files:
        #Check that we aren't going to clobber existing gcov files
        if len(glob.glob("*.gcov")):
            sys.stderr.write("There already exists *.gcov files in this directory. Please delete them first. \nThis script doesn't play well with others ;).\n")
            sys.exit()

        build_d = findbuilddir()

        for dirpath, dirnames, filenames in os.walk("src"):
            ccfiles = [ f for f in filenames if os.path.splitext(f)[1] == ".cc" ]
            if not ccfiles:
                continue
            stem = dirpath[4:] # get rid of src and path separator
            gcovcall = ["gcov", "-p", "-o", os.path.join(build_d,stem)] + ccfiles
            sys.stdout.write("Calling: " + ' '.join(gcovcall) + "\n" )
            subprocess.call(gcovcall)
            for filename in glob.glob("src*.gcov"):
                coverageinfo.process(filename)
            #Cleanup (in case the next invocation doesn't overwrite the files
            for filename in glob.glob("*.gcov"):
                os.remove(filename)

    if options.genuserinfo and not coverageinfo.userinfo:
        coverageinfo.genuserinfo()

    if options.save:
        coverageinfo.save(options.save)

    coverageinfo.genrevlist()

    coverageinfo.filestats()

    if options.byuser:
        coverageinfo.userstats()

    if options.userlist:
        for username in options.userlist:
            coverageinfo.filestats(username)

# Better handle multiple paragraph descriptions.
class PreformattedDescFormatter (IndentedHelpFormatter):
    def format_description(self, description):
        return description.strip() + "\n" # Remove leading/trailing whitespace

if __name__ == "__main__":
    parser = OptionParser(usage="usage: %prog [options] [usernames]",
        description="""This program runs the gcc gcov code coverage tool on the Rosetta codebase and
outputs summary statistics. It assumes that the running statistic have already
been produced.

To gather coverage statistics, first you have to compile Rosetta with the gcc
compiler and "extras=gcov". (This places the appropriate *.gcno files in the
build directory). Note that the statistics are a little better if you compile
without optimization, though this is probably not strictly necessary.

$ scons -j8 extras=gcov bin && scons -j8 cat=test extras=gcov

(Note that if you need to recompile, you'll need to clear out the build
directory and recompile - unfortunately gcov has internal checks that don't
work well with incremental compiles.)

Then you should run whatever conditions you wish to test. (For example, unit tests
and/or integration tests). The number of calls will sum over all subsequent invocation.
If you want to reset, just remove all the *.gcda files from the build directory.

$ test/run.py -j8 --extras=gcov --mute all
$ ./integration.py -j8 --extras=gcov

Then you can run this script to get the statistics for the code.

By default, this script ignores lines with asserts (regular, runtime or Py)
or utility exits in it (but not all lines in blocks that will utility exit).
The philosophy being that normal runs shouldn't execute those lines anyway,
so they shouldn't count against coverage statistics. This check is done with
a simple substring match, so any line with the strings "assert" or
"utility_exit" will not be counted. Use the -x flag to change this behavior.

You can also output files with per user statistics with the -u option, or
statistics for just those lines attributed to particular users by specifying
their name according to git on the commandline (This is their fullname,
rather than the github or email address).
Note this uses git utilities, so this should be run from a git checkout directory.

Generating the running report and the user statistics is rather computationally
intensive, so you can cache the results of the preload with the -s/-l options
to save time later, if you want to do more detailed analysis. (e.g. with the
username or -r options). Note the git user statistics are only computed and
saved if they're needed.

Finally, for ease of use, you can invoke the relevant gcov with the -f option
(if present this will not run anything else), which will produce a breakdown
of the executed/non-executed lines in a given file.

""",formatter=PreformattedDescFormatter())
    parser.add_option("-x", "--exits", action="store_true", help="don't ignore lines with asserts or utility_exits")
    parser.add_option("-a", "--alphabetical", action="store_true", help="sort output alphabetically by full path, rather than by percentage unexecuted")
    parser.add_option("-n", "--bylines", action="store_true", help="sort output by number of unexecuted lines, rather than by percentage")
    parser.add_option("-i", "--invert", action="store_true", help="Print executed rather than non-executed statistics")
    parser.add_option("-u", "--byuser", action="store_true", help="also output per-committer (via git blame) statistics")
    parser.add_option("-r", "--rev", type="string", help="Only count lines from commits more recent than the given revision.")
    parser.add_option("-s", "--save", help="Save the interpreted data state to FILE, so subsequent runs don't have to rerun gcov/git", metavar="FILE")
    parser.add_option("-l", "--load", help="Load the interpreted data state from FILE, so you don't have to rerun gcov/git", metavar="FILE")
    parser.add_option("-f", "--file", help="Don't run any other option, but generate the gcov output files for FILE and associated headers", metavar="FILE")

    options, args = parser.parse_args(sys.argv[1:])
    options.userlist = args
    options.genuserinfo = options.byuser or (len(options.userlist) != 0) or (options.rev is not None)

    #Check to make sure we're in the rosetta_source/ directory.
    if not os.path.isdir("src") or not os.path.isdir("build"):
        sys.stderr.write("Script must be invoked from the rosetta_source/ directory!\n")
        sys.exit()

    #Check that the gcov program is availble
    devnull = open(os.devnull, "w")
    try:
        subprocess.call(["gcov", "-h"], stdout = devnull, stderr = devnull)
    except Exception:
        sys.stderr.write("Issue calling the gcov executable - make sure it's installed and on your path.\n")
        sys.exit()
    devnull.close()

    #Check that git blame is availible, if we're doing user
    if( options.genuserinfo ):
        devnull = open(os.devnull, "w")
        try:
            if subprocess.call(["git", "status", "-u", "no"], stdout = devnull, stderr = devnull):
                sys.stderr.write("The current directory doesn't look to be a git directory. git needed for user statistics.\n")
                sys.exit()
        except Exception:
            sys.stderr.write("Issue calling git. Need it to do user or revision statistics.\n")
            sys.exit()
        devnull.close()

    main(options)

