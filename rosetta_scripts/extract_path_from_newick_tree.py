#!/usr/bin/env python2

#See https://www.rosettacommons.org/docs/wiki/scripting_documentation/RosettaScripts/multistage/TimeMachineExample
#for an example of how this is used

import sys

if len( sys.argv ) != 3:
    print "The first argument should be the element (leaf) you want to print the lineage for"
    print "The second argument should be the tree.\n"
    print "Example:"
    print "python extract_path_from_pewick_tree.py JR_105_1 '((((JR_113_1)JR_107_1,(JR_114_1)JR_108_1)JR_54_1,((JR_111_1)JR_101_1,(JR_112_1)JR_102_1)JR_77_1,(JR_105_1)JR_97_1)input_source_1)all'"

def get_lineage( leaf, tree ):
    scope = 0
    next_scope = 0
    #Add whitespace to help with parsing
    tree = tree.replace( "(", " ( " )
    tree = tree.replace( ")", " ) " )
    tree = tree.replace( ",", " , " )
    elements = tree.split()
    for i in range( 0, len( elements ) ):
        if elements[ i ] == "(":
            scope += 1
        elif elements[ i ] == ")":
            scope -= 1
            if scope == next_scope and scope != 0:
                print elements[ i+1 ]
                next_scope = scope - 1
        elif elements[ i ] == leaf:
            print leaf
            next_scope = scope - 1



leaf = sys.argv[ 1 ]
tree = sys.argv[ 2 ]

get_lineage( leaf, tree )
