import os, sys, re
import make_serialize_templates

sys.path.insert( 0, os.path.realpath(__file__).rpartition("/")[0]+"/../python_cc_reader" )
sys.path.insert( 0, os.path.realpath(__file__).rpartition("/")[0]+"/../external" )
#print( sys.path )

import blargs
from python_cc_reader.cpp_parser import code_reader
from python_cc_reader.beauty import beautifier


if __name__ == '__main__':
    with blargs.Parser(locals()) as p :
        p.str( "definitions" ).shorthand( "d" ).required().described_as( "The file with the concatenated output of the class names and fields" ).required()
        p.str( "creator_file" ).described_as( "The file with each line being the namespace-scoped list of all the widget creator classes for the Widgets to which XML Schema routines must be added and the string key" ).required()
        p.str( "class_type" ).required()

    all_ok = True
    defs = make_serialize_templates.load_definitions( definitions )
    print("definitions loaded2")

    creator_classes = []
    class_keys = []
    for line in open( creator_file ).readlines() :
        if len(line) == 0 or line[0] == "#" : continue
        cols = line.split();
        assert( len(cols) == 2 )
        creator_classes.append( cols[0].strip() )
        class_keys.append( cols[1].strip() )

        if creator_classes[-1] not in defs.classes :
            print("Could not find class", creator_classes[-1], "in definitions file")
            all_ok = False


    widget_classes = []
    all_classes = list( creator_classes )
    for classname in creator_classes :
        widget_classname = classname.rpartition("Creator")[0]
        if widget_classname not in defs.classes :

            got_a_sub = False
            #print widget_classname, "and end-of-name", widget_classname[ (-1*len(class_type)): ]
            if len( widget_classname ) > len( class_type ) and widget_classname[ (-1*len(class_type)): ] == class_type :
                temp = widget_classname[:(-1*len(class_type))]
                #print temp
                if temp in defs.classes :
                    widget_classname = temp
                    #print "got it"
                    got_a_sub = True
            if not got_a_sub :
                widget_classname += class_type # e.g. "Mover" or "Filter"
                if widget_classname not in defs.classes :
                    print("Could not find class", widget_classname, "in definitions file")
                    all_ok = False
                    widget_classname = "FIXME"
        widget_classes.append( widget_classname )

    for ii in range(len(widget_classes)) :
        print(creator_classes[ ii ], class_keys[ ii ], widget_classes[ ii ])


        
