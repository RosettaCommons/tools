#!/usr/bin/env python3

'''
This script will find large files that are in a particular branch, but not in a reference branch (usually master).
It will also sum the size of all files (including deleted/modified ones) for the branch.
Note that the total space taken in the repo is likely to be much, much smaller than the printed size,
due to the compression which Git uses.
'''

import os, sys
import argparse
import subprocess
import json

# For simplicity, I keep things mostly as raw bytestrings throughout.

def get_blobs(branchname, refname):
    '''Returns dictionary of (bytestring) hashes of blobs to filename for the given (unicode string) branch name that are not also in refname.'''
    run = subprocess.run( "git rev-list --objects "+refname+".."+branchname , stdout=subprocess.PIPE, shell=True )
    retval = {}
    for line in run.stdout.split(b'\n'):
        sl = line.split()
        if len(sl) == 0:
            continue
        blob = sl[0]
        if len(sl) >=2:
            filename = sl[1]
        else:
            filename = b"<N/A>"
        retval[ blob ] = filename

    return retval

def blob_sizes(blobs):
    '''Returns a sorted list of (size, hash, filename ) for the input blob values.
    (size in bytes, hash in bytestring.)
    input is {hash: filename}
    '''
    inputdata = b'\n'.join( blobs.keys() )
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

        fn = blobs.get(n, "<N./A.>")
        retval.append( (size, n, fn) )

    retval.sort(reverse=True)

    return retval

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
    parser.add_argument('-b', '--branch', help="The branch to use (default:current branch)", default="HEAD")
    parser.add_argument('-r', '--ref', help="The reference branch (default:origin/master", default="origin/master")
    parser.add_argument('-j', '--json', help="If present, write results to the JSON file parameter.", default=None)
    parser.add_argument('-s', '--size', type=int, default=1024, help="Don't individually list files under the given size (in KB)" )
    args = parser.parse_args()

    blobs = get_blobs( args.branch, args.ref )
    sizes = blob_sizes( blobs )

    total_size = sum( s for (s,b,fn) in sizes )
    nblobs = len(sizes)

    json_output = {}
    json_output["total"] = total_size
    json_output["nentries"] = nblobs

    print("TOTAL SIZE:", pprint_size(total_size) )
    print("NENTRIES:", nblobs )
    print()

    files = {}
    for s, b, fn in sizes:
        if s <= args.size*1024:
            continue
        print( fn.decode(), " --", pprint_size(s) )

        files[ (fn+b':'+b).decode() ] = s

    json_output["files"] = files

    if args.json is not None:
        with open(args.json, 'w') as f:
            json.dump(json_output, f)
