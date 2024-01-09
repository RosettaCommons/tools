#!/usr/bin/env python3

'''
This script will find large files in the repository, and then attempt to print the branches on which they can be reached.
This includes branches on which the files have been added and then subsequently deleted. (As a simplification, the labels MAIN and FOLDIT will be printed for files which are in the main Rosetta branch and the main Foldit branch, respectively.)

The files will be printed out in reverse order of size. It may take a while to complete, so feel free to cancel the script when you get what you need.

NOTE: This script will run a git prune & git gc to make its job easier

TO TRACK DOWN COMMITS FOR A BLOB: `git whatchanged --all --find-object=${blob_SHA1}`
'''

import os, sys
import argparse
import subprocess


# For simplicity, I keep things mostly as raw bitstrings throughout.

def cleanup():
    subprocess.run( "git remote prune origin", shell=True )
    subprocess.run( "git gc", shell=True )

def find_blobs():
    '''Returns list of (bitstring) hashes of blobs in the repo'''
    run = subprocess.run( "git rev-list --objects --all", stdout=subprocess.PIPE, shell=True )
    retval = []
    for line in run.stdout.split(b'\n'):
        sl = line.split()
        if len(sl) == 0:
            continue
        retval.append(sl[0])

    return retval

def blob_sizes(blobs, cutoff):
    '''Returns a sorted list of (size, hash) for the input blob values.
    (size in bytes, hash in bytestring.)
    Won't include anything smaller than cutoff.
    '''
    inputdata = b'\n'.join(blobs)
    run = subprocess.run( b"git cat-file --batch-check='%(objecttype) %(objectsize) %(objectname)'", stdout=subprocess.PIPE, shell=True, input=inputdata )

    retval = []

    for line in run.stdout.split(b'\n'):
        ls = line.split()
        if len(ls) == 0:
            continue
        if len(ls) != 3:
            print("ERROR: can't interpret blob size line ", line)
            continue

        t, s, n = ls
        if t != b'blob':
            continue

        size = int( s.decode() );

        if size < cutoff:
            continue

        retval.append( (size, n) )

    retval.sort(reverse=True)

    return retval

def find_commits(blob):
    '''Find the commit(s) which involve a particular blob
    Returns list of tuples containing the (bitstring) hashs of the commit and the corresponding name.'''
    command = b"git whatchanged --all --oneline --no-renames --find-object=" + blob
    run = subprocess.run( command, stdout=subprocess.PIPE, shell=True )
    commits = []
    current_commit = None
    for line in run.stdout.split(b'\n'):
        if line.startswith(b':'):
            sl = line.split()
            if len(sl) <= 5:
                print("ERROR: With command ", command)
                print("ERROR: Line too short: ", line )
            elif current_commit is None:
                print("ERROR: With command ", command)
                print("ERROR: Couldn't find commit for line: ", line )
            elif sl[4] != 'D': # Don't bother with delete commits
                commits.append( (current_commit, sl[5] ) )
                continue # We may have an add on the same commit
            current_commit = None
        else:
            sl = line.split()
            if len(sl) >= 1:
                current_commit = line.split()[0]
            #Ignore empty lines.

    return commits

def find_branches(commit):
    '''Find a list of branches (returned as bitstring names) which contain a particular commit'''
    run = subprocess.run( b"git branch -r --contains " + commit, stdout=subprocess.PIPE, shell=True )
    return run.stdout.split()

def find_tags(commit):
    '''Find a list of tags (returned as bitstrting names) which contain a particular commit'''
    run = subprocess.run( b"git tag --contains " + commit, stdout=subprocess.PIPE, shell=True )
    return run.stdout.split()

def pprint_size(size):
    if size < 1024:
        return str(size)+"B"
    size //= 1024
    if size < 1024:
        return str(size)+"KiB"
    size //= 1024
    if size < 1024:
        return str(size)+"MiB"
    size //= 1024
    return str(size)+"GiB"

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-c', '--noclean', action="store_true", help="Don't clean up repo first")
    parser.add_argument('-s', '--size', type=int, default=1024,
                        help='The size (in KB) under which not to bother with files.')
    parser.add_argument('-m', '--nomain', action="store_true", help="Don't print out blobs labeled MAIN/FOLDIT")
    args = parser.parse_args()

    if not args.noclean:
        cleanup()

    print("Getting list of files")
    blobs = find_blobs()
    #blobs = [b'45b5e296da6b0c6a44c53e41819955c7bb9fd69a', b'afe4de164e1b34959ac7b7a272ee94ba6df31f30' ]

    print("Getting file sizes")
    sizes = blob_sizes( blobs, args.size*1024 )

    print("Finding branches")
    for size, blob in sizes:
        commits = find_commits(blob)

        branches = set()
        tags = set()
        fnames = set()
        for chash, fname in commits:
            fnames.add(fname)
            branches.update( find_branches(chash) )
            tags.update( find_tags(chash) )
            if b'origin/main' in branches:
                branches = {b'MAIN'}
                tags = set()
                break
            if b'origin/interactive/develop' in branches:
                branches = {b'FOLDIT'}
                tags = set()
                break

        if args.nomain and (b'MAIN' in branches or b'FOLDIT' in branches):
            continue

        print( "; ".join(f.decode() for f in fnames), "  ", blob.decode(), "--", pprint_size(size) )
        for branch in branches:
            print("\t",branch.decode())
        for tag in tags:
            print("\tTag:",tag.decode())

        print()




