#!/usr/bin/env python3

'''
This script will find large files that are in a particular branch, but not in a reference branch (usually main).
It will also sum the size of all files (including deleted/modified ones) for the branch.
Note that the total space taken in the repo is likely to be much, much smaller than the printed size,
due to the compression which Git uses.
'''

import os, sys, zlib
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


def calculate_size_using_blobs(args):
    blobs = get_blobs( args.branch, args.ref )
    sizes = blob_sizes( blobs )

    total_size = sum( s for (s,b,fn) in sizes )
    nblobs = len(sizes)

    json_output = {}
    json_output["total"] = total_size
    json_output["nentries"] = nblobs

    json_output['merge-head'] = subprocess.check_output(f'git rev-parse {args.branch}', shell=True).decode(encoding='utf-8', errors='backslashreplace').replace('\n', '')
    json_output['merge-base'] = subprocess.check_output(f'git rev-parse {args.ref}', shell=True).decode(encoding='utf-8', errors='backslashreplace').replace('\n', '')

    if args.verbose: print(f'Merge: {json_output["merge-base"]} ← {json_output["merge-head"]}')

    print("TOTAL SIZE:", pprint_size(total_size) )
    print("NENTRIES:", nblobs )
    print()

    files = {}
    for s, b, fn in sizes:
        if s <= args.size*1024:
            continue
        print( fn.decode(), ":", b.decode(), "--", pprint_size(s) )

        files[ (fn+b':'+b).decode() ] = s

    json_output["files"] = files

    if args.json:
        with open(args.json, 'w') as f:
            json.dump(json_output, f)


def calculate_size_using_diffs(args):
    merge_head = subprocess.check_output(f'git rev-parse {args.branch}', shell=True).decode(encoding='utf-8', errors='backslashreplace').replace('\n', '')
    merge_base = subprocess.check_output(f'git rev-parse {args.ref}', shell=True).decode(encoding='utf-8', errors='backslashreplace').replace('\n', '')

    print(f'merge: {merge_base} ← {merge_head}')

    commits = subprocess.check_output(f'git rev-list {merge_head} ^{merge_base}', shell=True).decode(encoding='utf-8', errors='backslashreplace')

    #print('Calculating commit-diff size for following commits:\n', commits)
    commits = commits.split()

    results = []
    for c in commits:
        diff = subprocess.check_output(f'git format-patch --stdout --binary -1 {c}', shell=True).decode(encoding='utf-8', errors='backslashreplace')
        compressed_diff = zlib.compress( diff.encode("utf-8") )

        results.append( dict(commit=c, diff=len(diff), compressed_diff=len(compressed_diff)) )


    print("Raw and compressed size of diff's for each commit:")
    results.sort(key = lambda c: -c['compressed_diff'] )
    for c in results: print(f'{c["commit"]} → raw: {pprint_size(c["diff"])}, compressed: {pprint_size(c["compressed_diff"])}')


    json_output = dict(
        commits = results,
        raw_diff_size = sum( ( c['diff'] for c in results ) ),
        compressed_diff_size = sum( ( c['compressed_diff'] for c in results ) ),
    )

    print(f'\ntotal-compressed-diff-size: {pprint_size(json_output["compressed_diff_size"])}, total-raw-diff-size: {pprint_size(json_output["raw_diff_size"])}')
    with open(args.json, 'w') as f:
         json.dump(json_output, f, sort_keys=True, indent=2)




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-b', '--branch', help="The branch to use (default:current branch)", default="HEAD")
    parser.add_argument('-r', '--ref', help="The reference branch (default:origin/main", default="origin/main")
    parser.add_argument('-j', '--json', help="If present, write results to the JSON file parameter.", default=None)
    parser.add_argument('-s', '--size', type=int, default=1024, help="Don't individually list files under the given size (in KB)" )
    parser.add_argument('-v', '--verbose', action="store_true", help="increase output verbosity")
    parser.add_argument('-d', '--diff', action="store_true", help="calculate size of a diff's")
    args = parser.parse_args()

    if args.diff: calculate_size_using_diffs(args)
    else: calculate_size_using_blobs(args)
