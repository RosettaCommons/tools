#!/usr/bin/env python
# Really simple function to create blank .cc, .hh, and .fwd.hh files
#   for movers. Adapted from create_rosetta_class.py. Probably
#   should not be copying code, but oh well.
# -- rhiju, june 2015

from sys import argv
import string
import argparse
from os import getcwd,system
from os.path import *

def Help():
    print argv[0], " <full_class_name>"
    print "Example of a full_class_name could be protocols/sample_stream/SampleStream"
    print "Must run this from a directory within rosetta/main..."
    exit( 0 )

args = argv
FORCE = 0

for i in range(len(argv)):
    if argv[i] == '-f' or argv[i]=='-force' or argv[i]=='--force':
        FORCE = 1
        del( argv[i] )

full_class_name = args[1]
full_class_name = full_class_name.replace( '::', ':' ).replace( ':', '/' ).replace( '.cc', '')
cols = string.split( full_class_name, '/' )

class_name = basename( full_class_name )

ok_namespaces = ['basic','core','devel','numeric','protocols','utility']
if not cols[0] in ok_namespaces:
    print cols[0], " not in ", ok_namespaces
    Help()

CWD = getcwd()
pos = CWD.find( 'main/' )
if pos < 0:    pos = CWD.find( 'tools' )
if pos < 0:    Help()
srcdir = newdir = CWD[:pos] + 'main/source/src/'
newdir = srcdir + dirname(full_class_name)

if not exists( newdir ):
    print "Making new directory: ", newdir
    system( 'mkdir -p '+newdir )
else:
    "Directory exists: ", newdir

author = 'Rhiju Das, rhiju@stanford.edu'

cc_file = full_class_name + '.cc'
hh_file = full_class_name + '.hh'
fwd_hh_file = full_class_name + '.fwd.hh'

abs_cc_file = srcdir + '/' + cc_file
abs_hh_file = srcdir + '/' + hh_file
abs_fwd_hh_file = srcdir + '/' + fwd_hh_file

def write_shared_header( fid, file_name ):
    fid.write( '// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-\n')
    fid.write( '// vi: set ts=2 noet:\n')
    fid.write( '//\n')
    fid.write( '// (c) Copyright Rosetta Commons Member Institutions.\n')
    fid.write( '// (c) This file is part of the Rosetta software suite and is made available under license.\n')
    fid.write( '// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.\n')
    fid.write( '// (c) For more information, see http://www.rosettacommons.org. Questions about this can be\n')
    fid.write( '// (c) addressed to University of Washington CoMotion, email: license@uw.edu.\n')
    fid.write( '\n')
    fid.write( '/// @file %s\n' % file_name)
    fid.write( '/// @brief \n' )
    fid.write( '/// @detailed\n')
    fid.write( '/// @author %s\n' % author )
    fid.write( '\n')
    fid.write( '\n')

namespace_cols = string.split( dirname( full_class_name ), '/' )

tracer_name = ''
for col in namespace_cols: tracer_name += col+'.'
tracer_name += class_name

if exists( abs_cc_file ) and not FORCE:
    print abs_cc_file, ' exists already. Use -force to override.'
else:
    fid = open( abs_cc_file, 'w' )
    write_shared_header( fid, cc_file )
    fid.write( '#include <%s>\n' % hh_file )
    fid.write( '\n')
    fid.write( '#include <basic/Tracer.hh>\n' )
    fid.write( '\n')
    fid.write( 'static basic::Tracer TR( "%s" );\n' % tracer_name )
    fid.write( '\n')
    fid.write( 'using namespace core;\n')
    fid.write( '\n')
    for col in namespace_cols: fid.write( 'namespace %s {\n' % col )
    fid.write( '\n')
    fid.write( '\t//Constructor\n' )
    fid.write( '\t%s::%s()\n' % (class_name, class_name) )
    fid.write( '\t{}\n')
    fid.write( '\n')
    fid.write( '\t//Destructor\n' )
    fid.write( '\t%s::~%s()\n' % (class_name, class_name) )
    fid.write( '\t{}\n')
    fid.write( '\n')
    fid.write( '\tvoid\n')
    fid.write( '\t%s::apply( core::pose::Pose & pose )\n' % (class_name) )
    fid.write( '\t{\n')
    fid.write( '\t}\n')
    fid.write('\n')
    for col in namespace_cols[::-1]: fid.write( '} //%s \n' % col )

    fid.close()
    print 'Created: ', abs_cc_file
    print 'Add to %s.src.settings if you dare... ' % namespace_cols[0]

if exists( abs_hh_file ) and not FORCE:
    print abs_hh_file, ' exists already. Use -force to override.'
else:
    fid = open( abs_hh_file, 'w' )
    write_shared_header( fid, hh_file )
    included_tag = 'INCLUDED_%s_HH'  %  full_class_name.replace( '/', '_' )
    fid.write( '#ifndef %s\n' % included_tag )
    fid.write( '#define %s\n' % included_tag )
    fid.write('\n')
    fid.write('#include <protocols/moves/Mover.hh>\n')
    fid.write('#include <%s>\n' % fwd_hh_file )
    fid.write('\n')
    for col in namespace_cols: fid.write( 'namespace %s {\n' % col )
    fid.write('\n')
    fid.write('\tclass %s: public protocols::moves::Mover {\n' % class_name)
    fid.write('\t\n')
    fid.write('\tpublic:\n')
    fid.write('\t\n')
    fid.write('\t\t//constructor\n')
    fid.write('\t\t%s();\n' % class_name)
    fid.write('\t\n')
    fid.write('\t\t//destructor\n')
    fid.write('\t\t~%s();\n' % class_name)
    fid.write('\t\n')
    fid.write('\tpublic:\n')
    fid.write('\t\n')
    fid.write('\t\n')
    fid.write('\tvirtual void apply( core::pose::Pose & pose );\n')
    fid.write('\t\n')
    fid.write('\tvirtual std::string get_name() const{ return "%s"; }\n' % class_name)
    fid.write('\t\n')
    fid.write('\tprivate:\n')
    fid.write('\t\n')
    fid.write('\t};\n')
    fid.write('\n')
    for col in namespace_cols[::-1]: fid.write( '} //%s \n' % col )
    fid.write('\n')
    fid.write( '#endif\n' )
    fid.close()
    print 'Created: ', abs_hh_file


if exists( abs_fwd_hh_file ) and not FORCE:
    print abs_fwd_hh_file, ' exists already. Use -force to override.'
else:
    fid = open( abs_fwd_hh_file, 'w' )
    write_shared_header( fid, fwd_hh_file )
    included_tag = 'INCLUDED_%s_FWD_HH'  %  full_class_name.replace( '/', '_' )
    fid.write( '#ifndef %s\n' % included_tag )
    fid.write( '#define %s\n' % included_tag )
    fid.write('\n')
    fid.write('#include <utility/pointer/owning_ptr.hh>\n')
    fid.write('\n')
    for col in namespace_cols: fid.write( 'namespace %s {\n' % col )
    fid.write('\t\n')
    fid.write('\tclass %s;\n' % class_name)
    fid.write('\ttypedef utility::pointer::shared_ptr< %s > %sOP;\n' % (class_name, class_name) )
    fid.write('\ttypedef utility::pointer::shared_ptr< %s const > %sCOP;\n' % (class_name, class_name) )
    fid.write('\t\n')
    for col in namespace_cols[::-1]: fid.write( '} //%s \n' % col )
    fid.write('\n')
    fid.write( '#endif\n' )
    fid.close()
    print 'Created: ', abs_fwd_hh_file











