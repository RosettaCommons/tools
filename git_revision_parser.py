#!/usr/bin/env python
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

# For help, run with -h

import sys, os, subprocess
from os import path
from optparse import OptionParser
import urllib2
from HTMLParser import HTMLParser
DEVNULL = open(os.devnull, 'wb')

class RevMapParser(HTMLParser):
    def __init__(self):
        HTMLParser.__init__(self)
        self.rev_map = {}
        self.reverse_revmap = {}
        self.revnum = None
    def handle_data(self, data):
        dat = data.strip()
        if len(data) == 40 and self.revnum is not None:
            #Is SHA1 hash
            self.rev_map[ self.revnum] = data
            self.reverse_revmap[ data ] = self.revnum
            self.revnum = None
        elif 5 <= len(data) <= 39: #
            try:
                val = int(data)
                #If it works, it's likely a revision number.
                self.revnum = val
            except ValueError:
                pass #Not a revision number, ignore.

def get_revmap():
    parser = RevMapParser();
    parser.feed( urllib2.urlopen("http://test.rosettacommons.org/revisions_map").read() )
    return parser.rev_map, parser.reverse_revmap

def get_mergebase(revision, target, message):
    process = subprocess.Popen( ["git","merge-base", revision, target],  stdout=subprocess.PIPE, stderr=subprocess.PIPE )
    out, err = process.communicate()
    if process.returncode != 0:
        print "ERROR: Can't find merge base for", message
        print "Did you run 'git fetch' first?"
        raise IOError
    return out.strip()

def get_revision(target_SHA, revmap, reverse_revmap):
    # Find merge base with earliest web-listed commit
    earliest_web_merge_base = get_mergebase(revmap[ min(revmap.keys()) ] , target_SHA, "earliest web-map revision")

    if( earliest_web_merge_base == revmap[ min(revmap.keys()) ] ):
        #We're somewhere in the midst of the web revisions.
        #Do a bisection sort
        webkeys = list(revmap.keys())
        webkeys.sort()
        last_good = 0
        first_bad = len(webkeys)
        while( first_bad - last_good > 1 ):
            test_commit = (first_bad + last_good)/2 #At least 2 difference - won't be coincident on last_good
            if( revmap[ webkeys[test_commit] ] == target_SHA ):
                #Exact match
                return (webkeys[test_commit], webkeys[test_commit])
            mergebase = get_mergebase( revmap[ webkeys[test_commit] ], target_SHA, "r%d" % webkeys[test_commit] )
            if mergebase == revmap[ webkeys[test_commit] ]:
                last_good = test_commit
            else:
                first_bad = test_commit
        if last_good == len(webkeys)-1:
            return ( webkeys[last_good], None )
        else:
            return ( webkeys[last_good], webkeys[first_bad] )
    else:
        #We're like pre-web commits. Get the log message and see if it's a svn-transitioned commit
        process = subprocess.Popen( ["git","log", "-1", target_SHA],  stdout=subprocess.PIPE, stderr=DEVNULL )
        log_message, err = process.communicate()
        if process.returncode != 0:
            print "ERROR: Can't find log message for desired commit."
            raise IOError
        location = log_message.find("https://svn.rosettacommons.org/source/trunk/rosetta@")
        if( location != -1 ):
            rev = log_message[location:].split('@')[1].split()[0]
            return (int(rev), int(rev)) #Exactly

        #We're possibly in the strange mists of pre-web, post svn.
        # 1acd87ab39e2f8d2960ded74e5412dbc641880a7 is the last svn commit I know of
        last_svn_merge_base = get_mergebase("1acd87ab39e2f8d2960ded74e5412dbc641880a7", target_SHA, "last known svn revision")
        if last_svn_merge_base == "1acd87ab39e2f8d2960ded74e5412dbc641880a7":
            return (55111, min(revmap.keys()) )

        #One last shot - we're on a branch off of the svn versions
        process = subprocess.Popen( ["git","log", "-1", last_svn_merge_base],  stdout=subprocess.PIPE, stderr=DEVNULL )
        log_message, err = process.communicate()
        if process.returncode != 0:
            print "ERROR: Can't find log message for putative svn commit."
            raise IOError
        location = log_message.find("https://svn.rosettacommons.org/source/trunk/rosetta@")
        if( location != -1 ):
            rev = log_message[location:].split('@')[1].split()[0]
            return (rev, None) # Is there a good way to get the next commit after?

    return (None, None)

def get_release(target):
    process = subprocess.Popen( ["git","branch", "-a"],  stdout=subprocess.PIPE, stderr=DEVNULL )
    branches, err = process.communicate()
    if process.returncode != 0:
        print "ERROR: Can't find branch information"
        raise IOError
    #The earliest merge base of the weeklies
    earliest_base = "161b1ce83a51cf4767e2d609652fef8893cce790"
    weeklies = []
    releases = [] #ignore for now
    for line in branches.split():
        if line.find("weekly_releases/2") != -1:
            #This is the current situation - change if things change
            desig = line.strip().split("/")[-1]
            if( desig[0] != "2" ):
                continue #management, etxc.
            desig = desig.replace("-wk", "_")
            year, week = desig.split('_')[:2]
            weeklies.append( ( int(year), int(week), line.strip() ) )
        if line.find("releases/rosetta") != -1:
            releases.append( line.strip() )

    #This all assumes that weeklies are simple (not interwoven) branches off of the master trunk
    weeklies.sort()
    if( len(weeklies) < 2 ):
        return "UNKNOWN"
    earliest_weekly_branch_point = get_mergebase( weeklies[0][2], weeklies[-1][2], "earliest weekly branch point" )

    with_earliest_mergebase = get_mergebase( weeklies[0][2], target, "earliest weekly" )
    with_latest_mergebase = get_mergebase( weeklies[-1][2], target, "latest weekly" )

    if with_earliest_mergebase != earliest_weekly_branch_point:
        return "Before "+str(weeklies[0][0])+"wk"+str(weeklies[0][1])

    last_before = 0
    first_after = len(weeklies)
    while( first_after - last_before > 1 ):
        test_commit = (first_after + last_before)/2 #At least 2 difference - won't be coincident on last_good
        mergebase = get_mergebase( weeklies[ test_commit ][2], target, str(weeklies[ test_commit ][0])+"wk"+str(weeklies[ test_commit ][1]) )
        if mergebase == with_latest_mergebase:
            #print weeklies[ test_commit ][0], weeklies[ test_commit ][1], "After"
            first_after = test_commit
        else:
            #print weeklies[ test_commit ][0], weeklies[ test_commit ][1], "Before"
            last_before = test_commit
    if first_after == len(weeklies)-1:
        return "After "+str(weeklies[-1][0])+"wk"+str(weeklies[-1][1])
    else:
        return ("Before "+str(weeklies[first_after][0])+"wk"+str(weeklies[first_after][1])+
                " After "+str(weeklies[last_before][0])+"wk"+str(weeklies[last_before][1]))

def main(argv):
    '''
USAGE: Print relevant information about a given SHA1 hash.

This program should be run within the "main" Rosetta repository.
Either pass it a revision designation (e.g. a SHA1 hash or a branch name)
or run it without arguments to use the current HEAD revision.

It will look at various references to determine where the commit stands with regards to tested and released revisions.

Note: It is best to do a "git fetch" prior to running the program to update all the references.
It also needs internet access to fetch publically availible web pages (but not GitHub info).
    '''
    parser = OptionParser(usage="usage: %prog [commit reference] ",)
    parser.set_description(main.__doc__)

    (options, args) = parser.parse_args(args=argv)

    if len(args) > 1:
        print main.__doc__
        return -1

    if len(args) == 1:
        target = args[0]
    else:
        target = "HEAD"

    # Check if we're in the main directory. 835245c1dac6178bab6934bbe2564eff1af30cd4 is a very early commit.
    # By omitting the last digit we're prompting git to fill it in.
    process = subprocess.Popen( ["git","rev-parse", "835245c1dac6178bab6934bbe2564eff1af30cd"],  stdout=subprocess.PIPE, stderr=DEVNULL )
    out, err = process.communicate()
    if process.returncode != 0 or out.strip() != "835245c1dac6178bab6934bbe2564eff1af30cd4":
        print "ERROR: Script must be run under the Rosetta main directory"
        return -1

    process = subprocess.Popen( ["git","rev-parse",target],  stdout=subprocess.PIPE, stderr=DEVNULL )
    target_SHA, err = process.communicate()
    target_SHA = target_SHA.strip()
    if process.returncode != 0:
        print "ERROR: Cannot find commit for designation", target
        return -1

    # If we can't find a merge base with it, we likely have a bad revision number.
    process = subprocess.Popen( ["git","merge-base", "835245c1dac6178bab6934bbe2564eff1af30cd4", target],  stdout=subprocess.PIPE, stderr=DEVNULL )
    out, err = process.communicate()
    if process.returncode != 0 or out.strip() != "835245c1dac6178bab6934bbe2564eff1af30cd4":
        print "ERROR: Commit designation", target, "does not appear to be a good Rosetta main commit."
        return -1

    # We're a decent commit - find the revision map.
    revmap, reverse_revmap = get_revmap();

    try:
        rev_min, rev_max = get_revision( target_SHA, revmap, reverse_revmap )
    except IOError:
        return -1

    # Now try to find the corresponding release
    try:
        release_designation = get_release( target_SHA )
    except IOError:
        return -1

    #Now print out result
    print "For commit:", target
    print "SHA1:", target_SHA
    print "Test server revision:",
    if rev_min == rev_max:
        if rev_min is None:
            print "UNKNOWN"
        else:
            print "r%d" % rev_min
    else:
        if rev_max is None:
            print "After r%d" % rev_min
        elif rev_max is None:
            print "Before r%d" % rev_max
        else:
            print "Before r%d, After r%d" % (rev_max, rev_min)
    print "Weekly releases:", release_designation

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
