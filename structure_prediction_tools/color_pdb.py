#!/usr/local/bin/python2.7

from argparse import ArgumentParser
from sys import exit, stderr, stdout
from os.path import basename

aalist=['ALA', 'ARG', 'ASN', 'ASP',
        'CYS', 'GLU', 'GLN', 'GLY',
        'HIS', 'ILE', 'LEU', 'LYS',
        'MET', 'PHE', 'PRO', 'SER',
        'THR', 'TRP', 'TYR', 'VAL']

def get_score_dict( lines, score_term ):
    """ """
    # get index first
    col = ""
    for l in lines.split("\n"):
        #print l
        if l.startswith("label"):
            ls = l.split()
            try :
                col = ls.index( score_term )
            except:
                stderr.write("ERROR: couldn't find [%s]. Check the following score terms.\n" % args.score_term )
                ctr=1
                for score_term in ls:
                    print ctr, score_term
                    ctr+=1 
                return False

    if not col:
        stderr.write("ERROR: couldn't find the score header - use score_jd2 -out:pdb to dump a pdb with per residue energies\n")
        return False

    dict = {}
    for l in lines.split("\n"):
        if not l.strip(): continue
        ls = l.split()
        if ls[0][:3] in aalist: # first 3 characters are amino acid long name

            try:
                if "_" in ls[0]:
                    pos = ls[0].split("_")[-1] # HIS_D_2
                else:
                    pos = ls[1]

                pos   = int( pos )
                score = float(ls[col])
                dict[ pos ] = score
            except:
                stderr.write("ERROR: something is wrong here - pos: %s; score: %s\n"%( pos, ls[col])) 
                exit()

    return dict

def parse_score_file( scorefile, pos_col, score_col, autofill, autofill_score, total_residue ):
    dict = {}

    # this is for the scorefile that only contains scores at certain positions 
    if autofill:
        for pos in range(1,total_residue+1):
            dict[ pos ] = autofill_score

    with open( scorefile, "r" ) as f:
        for l in f:
            ls = l.strip().split()
            try:
                pos   = int(ls[ pos_col-1 ])
                score = float(ls[ score_col-1 ])
                dict[ pos ] = score
            except:
                stderr.write("WARNING: unparsable line: %s; %s %s\n" %( l.strip(), ls[ pos_col - 1 ], ls[ score_col - 1 ] ) )
                continue

    return dict 


def color_pdb( args ):
    pdblines = ""
    scorelines = ""
    with open( args.pdb, "r" ) as f:
        for l in f:
            if l.startswith("ATOM"):
                pdblines += l
            else:
                scorelines += l

    # get score_dict from two sources
    if args.scorefile:
        score_dict = parse_score_file( args.scorefile, args.pos_col, args.score_col, args.autofill, args.autofill_score, args.total_residue )

    elif scorelines:
        score_dict = get_score_dict( scorelines, args.score_term )


    if not score_dict:
        return False, False
    stderr.write("Reading in %s scores\n" %len(score_dict.keys()))

    minscore = min( score_dict.values() )
    # only normalize score when minscore < 0
    if minscore >= 0:
        minscore = 0

    outpdblines   = "REMARK B-factor column as %s\n" % args.score_term
    outscorelines = "rsn %s %s_add_%s\n" %( args.score_term, args.score_term, minscore )
    prev_rsn = ""
    for l in pdblines.split("\n"):
        if l.startswith("ATOM"):
            rsn = int(l[22:26].strip())
            outpdblines += l[:59] + "%6.2f\n" % (score_dict[rsn]-minscore)
            if rsn != prev_rsn:
                outscorelines += "%3s %5.4f %5.4f\n" %( rsn, score_dict[rsn], score_dict[rsn]-minscore)
            prev_rsn = rsn

    return outpdblines, outscorelines


if __name__=="__main__":
    parser = ArgumentParser("")
    parser.add_argument("pdb", help="")
    parser.add_argument("-c", "--score_term", type=str, help="")
    parser.add_argument("--scorefile", type=str, help="")
    parser.add_argument("--pos_col", default=1, type=int, help="")
    parser.add_argument("--score_col", default=2, type=int, help="")
    parser.add_argument("--autofill", action="store_true", help="")
    parser.add_argument("--autofill_score", type=int, default=0, help="")
    parser.add_argument("--total_residue", type=int, default=100000, help="")
    args = parser.parse_args()

    pdblines, scorelines = color_pdb( args )
    if not pdblines:
        exit()

    if args.score_term:
        tag = args.score_term
    else:
        tag = args.scorefile


    outfile = open( basename( args.pdb.split(".pdb")[0] ) +"_cb_"+tag+".pdb", "w")
    outfile.write(pdblines)
    outfile.close()

    #if args.dump_per_rsd_score:
    outfile = open( basename( args.pdb.split(".pdb")[0] ) +"_per_rsd_"+tag+".sc", "w")
    outfile.write(scorelines)
    outfile.close()


