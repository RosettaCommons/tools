from code_utilities import *
from code_reader import *
import re

OP_regex1 = re.compile( "OP" )
OP_regex2 = re.compile( "owning_ptr" )

def class_needs_a_cc_dstor( cl ) :
    has_dstor  = False
    dstor_inhh = False
    has_ownptr = False
    for func in cl.functions :
        if func.name == "Destructor" :
            has_dstor = True
            if func.definition_present :
                dstor_inhh = True
                break
    for dmemb in cl.data_members :
        if OP_regex1.search( dmemb[ 0 ] ) or OP_regex2.search( dmemb[ 0 ] ) :
            print "   Found suspected OP", dmemb[ 0 ], "in class ", cl.name
            has_ownptr = True
    return has_ownptr and ( not has_dstor or dstor_inhh )

def class_needs_a_cc_copy_dstor( cl ) :


def find_all_classes_that_require_cc_dstors() :
    compilable_files, all_includes, file_contents = load_source_tree()

    hh_files = non_fwd_hh_subset( compilable_files )

    classes_needing_ccdstors = []

    for hh_file in hh_files :
        if hh_file == "devel/simple_options/option.hh":
            continue
        if hh_file.partition( "utility" )[ 2 ] :
            continue
        if hh_file.partition( "numeric" )[ 2 ] :
            continue
        if hh_file.partition( "ObjexxFCL" )[ 2 ] :
            continue
        print "Processing header", hh_file
        classes = read_classes_from_header( hh_file, file_contents[ hh_file ] )

        for cl in classes :
            print  " Processing class", cl.name
            if class_needs_a_cc_dstor( cl ) :
                classes_needing_ccdstors.append( cl )

    return classes_needing_ccdstors

if __name__ == "__main__" :
    classes_needing_ccdstors = find_all_classes_that_require_cc_dstors()
    for cl in classes_needing_ccdstors:
        print "CCDstor needed for class", cl.name, "declared in file", cl.file_declared_within

