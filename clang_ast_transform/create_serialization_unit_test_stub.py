import blargs
import re
import os
import os.path

def cereal_includes() :
    return """
#ifdef    SERIALIZATION
#include <core/id/AtomID.hh>
#include <core/scoring/func/HarmonicFunc.hh>

// Cereal headers
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#endif
""".replace( "    ","\t" )

def find_last_include( flines ) :
    include_re = re.compile( "\s*#include" )
    for count, line in reversed(list(enumerate( flines ))) :
        if include_re.match( line ) :
            return count
    print( "did not find any #includes" )

def insert_needed_includes( flines ) :
    last = find_last_include( flines )
    flines = flines[:(last+1)] + cereal_includes().splitlines(True) + flines[(last+1):]
    return flines

def serialization_stub() :
    return """
    void test_serialize_CLASSNAME() {
        TS_ASSERT( true ); // for non-serialization builds
#ifdef SERIALIZATION
        using namespace NAMESPACE;
        using namespace core::scoring::func;
        using namespace core::id;

        FuncOP some_func( new HarmonicFunc( 1, 2 ));
        ConstraintOP instance( new CLASSNAME() ); // serialize this through a pointer to the base class

        std::ostringstream oss;
        {
            cereal::BinaryOutputArchive arc( oss );
            arc( instance );
        }

        ConstraintOP instance2; // deserialize also through a pointer to the base class
        std::istringstream iss( oss.str() );
        {
            cereal::BinaryInputArchive arc( iss );
            arc( instance2 );
        }

        // make sure the deserialized base class pointer points to a CLASSNAME
        TS_ASSERT( utility::pointer::dynamic_pointer_cast< CLASSNAME > ( instance2 ));
        TS_ASSERT( *instance == *instance2 );
#endif // SERIALIZATION
    }

""".replace( "    ","\t" )

def insert_serialization_stub_into_test_suite( flines, namespace_name, class_name ) :
    class_end = find_class_end( flines )
    stub = serialization_stub()
    stub = stub.replace( "NAMESPACE", namespace_name ).replace( "CLASSNAME", class_name )
    flines = flines[:class_end] + stub.splitlines(True) + flines[class_end:]
    return flines

def find_class_end( flines ) :
    classdec_re = re.compile( "\s*class\s" )
    suite_re = re.compile( ".*CxxTest::TestSuite" )
    lcb_re = re.compile( ".*{")
    suite_begun = False
    class_begun = False
    for count, line in enumerate( flines ) :
        #print( "find_class_end:", count, suite_begun, class_begun )
        if suite_begun :
            if line == "};\n" :
                return count
        elif class_begun :
            if suite_re.match( line ) :
                suite_begun = True
            elif lcb_re.match( line ) :
                class_begun = False # wrong class
        elif classdec_re.match( line ) :
            class_begun = True
            if suite_re.match( line ) :
                suite_begun = True

    if suite_begun :
        print( "Could not locate class ending" )
    else :
        print( "Could not locate suite declaration!" )
    return len(flines)

def empty_cxxtest_file() :
    return """// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/constaints/FILENAME.cxxtest.hh
/// @brief  test suite for CLASSNAME
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>

// basic headers
#include <basic/Tracer.hh>

// unit headers
#include <PATH/FILENAME>


class CLASSNAMETests : public CxxTest::TestSuite
{
public:

};

""".replace( "    ","\t" )

def add_serialization_stub_to_existing( test_file, namespace, classname ) :
    flines = open( test_file ).readlines()
    flines = insert_needed_includes( flines )
    flines = insert_serialization_stub_into_test_suite( flines, namespace, classname )
    print(( "writing to", test_file ))
    open( test_file, "w" ).writelines( flines )

def create_new_test_suite_w_serialization_stub( test_file, namespace, classname, hhfile ) :
    flines = empty_cxxtest_file()
    flines = flines.replace( "NAMESPACE", namespace ).replace( "PATH", namespace.replace("::","/"))
    flines = flines.replace( "CLASSNAME", classname ).replace( "FILENAME", hhfile )
    flines = flines.splitlines(True)
    flines = insert_needed_includes( flines )
    flines = insert_serialization_stub_into_test_suite( flines, namespace, classname )
    print(( "writing to", test_file ))
    open( test_file, "w" ).writelines( flines )


def create_test_suite_for_file( namespace, hhfilename, classname ) :
    cwd = os.getcwd()
    path_to_source = cwd.partition("source")[0]
    filepath_relative = namespace.replace("::","/")
    filename_dir = path_to_source + "source/src/" + filepath_relative
    test_dir = path_to_source + "source/test/" + filepath_relative
    test_file = test_dir + "/" + ( hhfilename[:-3] if hhfilename.find( "/" ) == -1 else hhfilename.rpartition("/")[2][:-3] ) + ".cxxtest.hh"

    #print( path_to_source, filepath_relative, filename_dir, test_dir, test_file )

    if os.path.isfile( test_file ) :
        add_serialization_stub_to_existing( test_file, namespace, classname )
    else :
        create_new_test_suite_w_serialization_stub( test_file, namespace, classname,
            hhfilename if hhfilename.find("/") == -1 else hhfilename.rpartition("/")[2] )


if __name__ == "__main__" :
    with blargs.Parser(locals()) as p :
        p.str( "namespace" ).required()
        p.str( "filename" ).shorthand("f").required()
        p.str( "classname" ).shorthand("c").required()


    create_test_suite_for_file( namespace, filename, classname )
