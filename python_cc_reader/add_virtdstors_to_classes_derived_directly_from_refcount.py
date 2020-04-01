from python_cc_reader.cpp_parser.code_utilities import *
from python_cc_reader.cpp_parser.code_reader import *
import re

OP_regex1 = re.compile( "OP" )
OP_regex2 = re.compile( "owning_ptr" )

def class_derives_directly_from_refcount( cl ) :
    if not cl.base : return False
    if cl.base == "utility::pointer::ReferenceCount" : return True
    if cl.base == "pointer::ReferenceCount" : return True;
    if cl.base == "ReferenceCount" : return True;
    return False

def find_all_classes_derived_directly_from_refcount() :
    compilable_files, all_includes, file_contents = load_source_tree()

    hh_files = non_fwd_hh_subset( compilable_files )

    classes_deriving_from_refcount = []

    for hh_file in hh_files :
        if hh_file == "devel/simple_options/option.hh":
            continue
        if hh_file.partition( "utility" )[ 2 ].partition( "pointer" )[ 2 ] :
            continue
        if hh_file.partition( "ObjexxFCL" )[ 2 ] :
            continue
        #print "Processing header", hh_file
        classes = read_classes_from_header( hh_file, file_contents[ hh_file ] )

        for cl in classes :
            #print  " Processing class", cl.name
            if class_derives_directly_from_refcount( cl ) :
                classes_deriving_from_refcount.append( cl )

    return classes_deriving_from_refcount

def dstor_for_class( cl ) :
    for func in cl.functions :
        if func.name == "Destructor" :
            return func
    return None

def class_has_virtual_dstor( cl ) :
    func = dstor_for_class( cl )
    if func :
        return func.is_explicitly_virtual
    return False

def class_has_dstor( cl ) :
    return not dstor_for_class( cl ) == None

def class_has_dstor_in_header( cl ) :
    func = dstor_for_class( cl )
    if func :
        return func.definition_present
    return False

def rewrite_header_to_declare_dstor_virtual( cl ) :
    cr = HeavyCodeReader()
    cr.push_new_file( cl.file_declared_within )
    flines = open( cl.file_declared_within ).readlines()
    newlines = []
    done = False
    for line in flines :
        print(line, end=' ')
        cr.examine_line( line )
        newlines.append( line )
        if done : continue
        if cr.at_class_scope() and cr.current_class_name() == cl.name :
            if cr.just_parsed_a_class_function_declaration() and cr.function_name == "Destructor" :
                assert ( not cr.function_is_explicitly_virtual )
                linecount = len(newlines)
                while linecount >= 0 :
                    linecount -= 1
                    tildepos = newlines[ linecount ].find( "~" )
                    if tildepos != -1 :
                        newline = newlines[ linecount ][ : tildepos ]  + "virtual " + newlines[ linecount ][ tildepos : ]
                        print(newline)
                        newlines[ linecount ] = newline
                        done = True
                        break
                if linecount < 0 :
                    print("Error looking for destructor in class " + cl.name + " declared within file " + cl.file_declared_within)
                    sys.exit(1)
    open( cl.file_declared_within, "w" ).writelines( newlines )

def rewrite_header_to_include_virtdstor( cl ) :
    cr = HeavyCodeReader()
    cr.push_new_file( cl.file_declared_within )
    flines = open( cl.file_declared_within ).readlines()
    newlines = []
    done = False
    for line in flines :
        newlines.append( line )
        if done : continue
        cr.examine_line( line )
        if cr.at_class_scope() and cr.current_class_name() == cl.name and cr.current_privacy_level() == "public" and cr.statement_string == "" :
            newlines.append( "\t///@brief Automatically generated virtual destructor for class deriving directly from ReferenceCount\n")
            newlines.append( "\tvirtual ~" + cl.name + "();\n" )
            done = True
    if not done : 
        print("Error: unable to find an appropriate place to add a virtual destructor to class", cl.name, "from", cl.file_declared_within)
        return False
    open( cl.file_declared_within, "w" ).writelines( newlines )
    return True
            

def rewrite_header_to_replace_dstor_implementation_w_dstor_declaration( cl ) :
    cr = HeavyCodeReader()
    cr.push_new_file( cl.file_declared_within )
    flines = open( cl.file_declared_within ).readlines()
    newlines = []
    dstor = dstor_for_class( cl )
    if not dstor :
        print("ERROR: could not find destructor for class " + cl + " declared within '" + cl.file_declared_within + "'")
        return
    if dstor.definition_begin_line == -1 :
        print("ERROR: dstor for class " + cl + " is not defined within the header file '" + cl.file_declared_within + "'")
    done = False
    middle_of_dstor_removal = False
    dstorstring = "~" + cl.name
    for line in flines :
        cr.examine_line( line )
        if not middle_of_dstor_removal :
            newlines.append( line )
        else :
            if cr.curr_line() > dstor.definition_end_line :
                middle_of_dstor_removal = False
                newlines.append( line )
                done = True
            else :
                newlines.append( "// " + line )
        if done : continue
        #print cr.current_class_name(), cl.name, line,
        if not cr.current_class_name() == cl.name : continue
        #print dstorstring, line, cr.commentless_line,
        if cr.commentless_line.find( dstorstring ) != -1 :
            #print "Found destructor line", cr.curr_line(), dstor.definition_begin_line
            if cr.curr_line() == dstor.definition_begin_line :
                # ok: we are on the line where the definition begins
                # if it also ends on this line, then we should just replace the "{...}" with ";// {...}"
                # if it spans multiple lines, then lets remove the entire declaration
                if dstor.definition_begin_line == dstor.definition_end_line :
                    newlines[ cr.curr_line() - 1 ] = line[:line.find("~")] + "virtual " + line[ line.find("~") : line.find( "{" ) ] + "; // auto-removing definition from header" + line[ line.find( "{" ) : ]
                    done = True;
                    continue
                else :
                    newlines[ cr.curr_line() - 1] = "\t///@brief Autogenerated destructor definition\n"
                    newlines.append( "\tvirtual ~" + cl.name + "();\n" )
                    newlines.append( "\n" )
                    newlines.append( "//" + line )
                    middle_of_dstor_removal = True
                    continue
        if cr.curr_line == dstor.definition_end_line :
            # somehow we missed the declaration of the dstor, maybe because it was declared "~ ClassName"
            # go back until we find a line that has a "~"
            countback = cr.curr_line
            while countback > 0 :
                --countback
                if newlines[ countback ].find( "~" ) != -1 :
                    break
            if countback == 0 :
                print("ERROR: could not find beginning of destructor for class " + cl.name + " in file " + cl.file_declared_within)
                return
            for linenum in range( countback, cr.curr_line ) :
                newlines[ linenum ] = "// " + newlines[ linenum ]
            newlines.append( "\t///@brief Autogenerated destructor definition\n" )
            newlines.append( "\tvirtual ~" + cl.name + "();\n" )
            done = True
    if not done :
        print("ERROR: failed to replace destructor implemenation with a (virtual) destructor declaration in class", cl.name, "in header file", cl.file_declared_within)
        return
    open( cl.file_declared_within, "w" ).writelines( newlines )


# just put the dstor at the top of the ccfile
def add_empty_dstor_to_ccfile( cl ) :
    cr = HeavyCodeReader()
    ccfile = cl.file_declared_within[:-2]+"cc"
    if not os.path.isfile( ccfile ) :
        print("Unable to locate ccfile '" + ccfile + "' to place explicit virtual dstor within for class " + cl.name)
        return False
    cr.push_new_file( ccfile )
    flines = open( ccfile ).readlines()
    newlines = []
    done = False
    target_scope = cl.scope.rpartition( "::" )[ 0 ]
    if target_scope == cl.name :
        target_scope = "Global"
    for line in flines :
        cr.examine_line( line )
        newlines.append( line )
        if done : continue
        #print cr.full_scope_name(), target_scope, line,
        if cr.full_scope_name() == target_scope :
            newlines.append( "\n" )
            newlines.append( "/// @details Auto-generated virtual destructor\n" )
            newlines.append( cl.name + "::~" + cl.name + "() {}\n" )
            done = True
    if not done :
        print("ERROR: Failed to find the appropriate place to place the destructor for class", cl.name, "in", ccfile)
        print("seeking scope: '" + target_scope + "'")
        return False
    open( ccfile, "w" ).writelines( newlines )
    return True

if __name__ == "__main__" :

    if len(sys.argv) > 1 :
        for focused_header in sys.argv[1:] :
            classes_in_focused_header = read_classes_from_header( focused_header, open( focused_header ).readlines() )
            classes_in_focused_header_deriving_from_refcount = [x for x in classes_in_focused_header if class_derives_directly_from_refcount(x)]
            classes_needing_virtccdstors = [x for x in classes_in_focused_header_deriving_from_refcount if not class_has_virtual_dstor(x)]
            classes_w_nonvirtdstors = [x for x in classes_needing_virtccdstors if class_has_dstor( x )]
            classes_wo_dstors = [x for x in classes_needing_virtccdstors if not class_has_dstor( x )]
            classes_w_nonvirtdstors_in_header = [x for x in classes_w_nonvirtdstors if class_has_dstor_in_header(x)]
            classes_w_nonvirtdstors_in_ccfile = [x for x in classes_w_nonvirtdstors if not class_has_dstor_in_header(x)]

            #for cl in classes_w_nonvirtdstors_in_ccfile:
            #    print "Class", cl.name, "declared in file", cl.file_declared_within, "requires its destructor, declared within its header, to become virtual"
            #    rewrite_header_to_declare_dstor_virtual( cl )
            #
            #sys.exit(0)

            for cl in classes_wo_dstors :
                print("Processing class", cl.name)
                added_new_dstor_to_cc = add_empty_dstor_to_ccfile( cl )
                if added_new_dstor_to_cc :
                    rewrite_header_to_include_virtdstor( cl )

            #for cl in classes_w_nonvirtdstors_in_header:
            #    print "processing class", cl.name, "in file", cl.file_declared_within
            #    moved_dstor_to_cc = add_empty_dstor_to_ccfile( cl )
            #    if moved_dstor_to_cc :
            #        print "attempting to replace dstor implementation in ", cl.name
            #        rewrite_header_to_replace_dstor_implementation_w_dstor_declaration( cl )
        sys.exit(0)

    classes_deriving_from_refcount = find_all_classes_derived_directly_from_refcount()
    classes_needing_virtccdstors = [x for x in classes_deriving_from_refcount if not class_has_virtual_dstor(x)]
    classes_wo_dstors = [x for x in classes_needing_virtccdstors if not class_has_dstor( x )]
    classes_w_nonvirtdstors = [x for x in classes_needing_virtccdstors if class_has_dstor( x )]
    classes_w_nonvirtdstors_in_header = [x for x in classes_w_nonvirtdstors if class_has_dstor_in_header(x)]
    classes_w_nonvirtdstors_in_ccfile = [x for x in classes_w_nonvirtdstors if not class_has_dstor_in_header(x)]

    for cl in classes_wo_dstors :
        print("Processing class", cl.name, "from", cl.file_declared_within)
        added_new_dstor_to_cc = add_empty_dstor_to_ccfile( cl )
        if added_new_dstor_to_cc :
            rewrite_header_to_include_virtdstor( cl )

    #for cl in classes_needing_virtccdstors :
    #    print "Class", cl.name, "needs a virtual destructor"
    #sys.exit(0)

    #for cl in classes_w_nonvirtdstors_in_ccfile:
    #    print "Class", cl.name, "declared in file", cl.file_declared_within, "requires its destructor, declared within its header, to become virtual"
    #    #rewrite_header_to_declare_dstor_virtual( cl )
    #
    #for cl in classes_w_nonvirtdstors_in_header:
    #    print "Class", cl.name, "declared in file", cl.file_declared_within, "requries a virtual destructor, currently non-virtual and declared within its header, to become virtual and for its implementation to move to its .cc file"
    #    moved_dstor_to_cc = move_dstor_to_ccfile( cl )
    #    if moved_dstor_to_cc :
    #        rewrite_header_to_replace_dstor_implementation_w_dstor_declaration( cl )
