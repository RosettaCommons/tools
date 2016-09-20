import os, sys, re
import make_serialize_templates
import blargs
import beautifier

sys.path.insert( 0, os.path.realpath(__file__).rpartition("/")[0]+"/../" )
#print( sys.path )
import python_cc_reader.code_reader

def class_hh_stub() :
    return [
        "\tstd::string\n",
        "\tget_name() const override;\n",
        "\n",
        "\tstatic\n",
        "\tstd::string\n",
        "\tmover_name();\n",
        "\n",
        "\tstatic\n",
        "\tvoid\n",
        "\tprovide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );\n" ]

def creator_hh_stub() :
    return [
        "\tprotocols::moves::MoverOP create_mover() const override;\n",
        "\tstd::string keyname() const override;\n",
        "\tvoid provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;\n" ]

def cc_stub() :
    return [
        "std::string %(classname)s::get_name() const {\n",
        "\treturn mover_name();\n",
        "}\n",
        "\n",
        "std::string %(classname)s::mover_name() {\n",
        "\treturn \"%(classkey)s\";\n",
        "}\n",
        "\n",
        "void %(classname)s::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )\n",
        "{\n",
        "// TO DO!\n",
        "\tusing namespace utility::tag;\n",
        "\tAttributeList attlist; // TO DO: add attributes to this list\n",
        "\t// TO DO: perhaps this is not the right function to call? -- also, delete this comment\n",
        "\t%(defaultschemafunc)s( xsd, mover_name(), attlist );\n",
        "}\n",
        "\n",
        "std::string %(creatorname)s::keyname() const {\n",
        "\treturn %(classname)s::mover_name();\n",
        "}\n",
        "\n",
        "protocols::moves::MoverOP\n",
        "%(creatorname)s::create_mover() const {\n",
        "\treturn protocols::moves::MoverOP( new %(classname)s );\n",
        "}\n",
        "\n",
        "void %(creatorname)s::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const\n",
        "{\n",
        "\t%(classname)s::provide_xml_schema( xsd );\n",
        "}\n" ]

class FunctionDecl :
    def __init__( self ) :
        self.line_begin = 0
        self.line_end = 0
        self.name = ""
        self.scope = ""
        self.privacy = ""


def reblob_tokens( toks ) :
    blobs = []
    working_blob = []
    for i in xrange(len( toks )) :
        working_blob.append( toks[i].spelling )
        if i+1 != len(toks) :

            if toks[i].start + len(toks[i].spelling) == toks[i+1].start :
                pass
            else :
                blobs.append( "".join( working_blob ))
                working_blob = []
    blobs.append( "".join( working_blob ))
    return blobs

def tokens_last_descendant( tok ) :
    if not tok.children : return tok
    else : return tokens_last_descendant( tok.children[-1] )

def turn_function_prefix_into_name( toks ) :
    funcname = []
    blobs = reblob_tokens( toks )
    #print "blobs", blobs
    return_type = ""
    func_name = ""
    ignore = set([ "static", "virtual" ]);
    for blob in blobs :
        if blob in ignore : continue
        if not return_type :
            return_type = blob
        elif not func_name and ( blob == "const" or blob == "&" or blob == "*" ) :
            return_type = return_type + " " + blob
        elif not func_name :
            func_name = blob
        else :
            func_name = func_name + " " + blob
    if not func_name :
        return return_type
    return func_name

def functions_for_file( file_lines, file_name ) :
    print "reading functions from", file_name
    beaut = beautifier.Beautifier()
    beaut.filename = file_name
    for line in file_lines :
        beaut.tokenize_line( line )
    beaut.minimally_parse_file()

    privacy = "public"
    ns_stack = []
    funcs = []
    #print classname, "ntoks: ", len(beaut.all_tokens)

    for tok in beaut.all_tokens :
        if not tok.is_visible : continue
        #print "tok:", tok.spelling, tok.type, tok.context(), ns_stack
        if tok.type == "namespace" :
            for child in tok.children :
                if child.type != "namespace-scope" :
                    #print "namespace child:", child.spelling
                    ns_stack.append( child.spelling )
            continue
        if tok.type == "class" :
            privacy = "private"
            for child in tok.children :
                if child.type != "class-scope" and child.spelling != ";" and child.type != "class-inheritance-list" :
                    #print "class child: ", child.spelling
                    ns_stack.append( child.spelling )
            continue
        if tok.type == "namespace-scope" :
            continue
        if tok.context() == "class" and tok.spelling == ";" :
            ns_stack.pop()
            continue
        if tok.context() == "namespace-scope" or tok.context() == "class-scope" :
            if tok.spelling == "typdef" : continue
            #print "token:", tok.type, tok.spelling, tok.line_number
            if tok.type == "class-privacy" :
                privacy = tok.spelling
                continue

            # at namespace scope
            # read until we get the function name
            child_funcname_toks = [ tok ]

            if tok.children :
                for child in tok.children :
                    if child.type != "function-decl-argument-list" :
                        child_funcname_toks.append( child )
                    else :
                        break
                ns_scope =  "::".join( ns_stack )
                unscoped_funcname = turn_function_prefix_into_name( child_funcname_toks )
                fname = ns_scope + "::" + unscoped_funcname
                func = FunctionDecl()
                func.line_begin = tok.line_number
                func.line_end = tokens_last_descendant(tok).line_number
                func.name = fname
                func.scope = "::".join( ns_stack ) + ( "" if ( unscoped_funcname.find( "::" ) == -1 ) else "::" + unscoped_funcname.rpartition( "::" )[0] )
                func.privacy = privacy
                #print "func: ", func.name, func.privacy, func.line_begin, func.line_end
                funcs.append( func )

            if tok.context() == "namespace-scope" and tok.parents_last_child() :
                ns_stack.pop()
    return funcs

def comment_out_function( func, lines ) :
    for ii in xrange( func.line_begin, func.line_end+1 ) :
        lines[ ii ] = "// XRW TEMP " + lines[ ii ]
    return lines

def add_include_at_bottom_of_includes( lines, include_lines ) :

    beaut = beautifier.Beautifier()
    beaut.pound_if_setting = "take_neither"
    for line in lines :
        beaut.tokenize_line( line )
    beaut.minimally_parse_file()

    last_include_line = -1
    for ii in xrange(len(lines)) :
        if beaut.line_tokens[ii] and beaut.line_tokens[ii][ 0 ].invisible_by_macro : continue
        if len(lines[ ii ]) > 7 and lines[ ii ][0:8] == "#include" :
            last_include_line = ii
    assert( last_include_line != -1 )
    lines = lines[:(last_include_line+1)] + include_lines + lines[ (last_include_line+1): ]
    return lines    

if __name__ == '__main__':
    with blargs.Parser(locals()) as p :
        p.str( "definitions" ).shorthand( "d" ).required().described_as( "The file with the concatenated output of the class names and fields" ).required()
        p.str( "triple_file" ).described_as( "The file with each line being the namespace-scoped list of all the widget creator classes for the Widgets to which XML Schema routines must be added, the string key, and the name-space scoped name of the Mover itself" ).required()
        p.str( "default_helper_func" ).required()
        p.multiword( "extra_cc_includes" ).cast( lambda x : [ y for y in x.split() ] )

    defs = make_serialize_templates.load_definitions( definitions )
    print "definitions loaded2"

    creator_classes = []
    class_keys = []
    widget_classes = []
    all_ok = True
    for line in open( triple_file ).readlines() :
        if len(line) == 0 or line[0] == "#" : continue
        cols = line.split();
        assert( len(cols) == 3 )
        creator_classes.append( cols[0].strip() )
        class_keys.append( cols[1].strip() )
        widget_classes.append( cols[2].strip() )

        if creator_classes[-1] not in defs.classes :
            print "Could not find class", creator_classes[-1], "in definitions file"
            all_ok = False
        if widget_classes[-1] not in defs.classes :
            print "Could not find class", widget_classes[-1], "in definitions file"
            all_ok = False

    if not all_ok :
        sys.exit(1)

    for ii in xrange(len(creator_classes)) :
        creator_name = creator_classes[ ii ]
        widget_name = widget_classes[ ii ]
        class_key = class_keys[ ii ]
        print "processing", widget_name, creator_name

        creator_hh_fname = defs.classes[ creator_name ].filename
        creator_cc_fname = creator_hh_fname[ :-2 ] + "cc"
        widget_hh_fname = defs.classes[ widget_name ].filename
        widget_cc_fname = widget_hh_fname[ :-2 ] + "cc"
        #print "files: ", creator_hh_fname, widget_hh_fname, widget_cc_fname

        widget_hh_lines = open( widget_hh_fname ).readlines()
        widget_cc_lines = open( widget_cc_fname ).readlines()
        creator_hh_lines = open( creator_hh_fname ).readlines()
        creator_cc_lines = None
        if ( os.path.isfile( creator_cc_fname ) ) :
            creator_cc_lines = open( creator_cc_fname ).readlines()

        widget_hh_funcs  = functions_for_file( widget_hh_lines, widget_hh_fname )
        widget_cc_funcs  = functions_for_file( widget_cc_lines, widget_cc_fname )
        creator_hh_funcs = functions_for_file( creator_hh_lines, creator_hh_fname )
        creator_cc_funcs = None
        if creator_cc_lines : creator_cc_funcs = functions_for_file( creator_cc_lines, creator_cc_fname )

        # tasks:
        # 1. comment out everything in creator.hh
        # 2. comment out existing get_name function
        # 3. add #inclusion to appropriate files

        widget_func_names_to_remove = []
        widget_func_names_to_remove.append( widget_name + "::get_name" )
        widget_func_names_to_remove.append( widget_name + "::class_name" )
        widget_func_names_to_remove.append( widget_name + "::mover_name" )


        widget_cc_funcs_to_remove = []
        for func in widget_cc_funcs :
            if func.scope == creator_name :
                widget_cc_funcs_to_remove.append( func )
            if func.name in widget_func_names_to_remove :
                widget_cc_funcs_to_remove.append( func )

        widget_hh_funcs_to_remove = []
        for func in widget_hh_funcs :
            if func.name in widget_func_names_to_remove :
                widget_hh_funcs_to_remove.append( func )

        creator_hh_funcs_to_remove = []
        for func in creator_hh_funcs :
            if func.scope == creator_name :
                creator_hh_funcs_to_remove.append( func )

        creator_cc_funcs_to_remove = []
        if creator_cc_funcs :
            for func in creator_cc_funcs :
                if func.scope == creator_name :
                    creator_cc_funcs_to_remove.append( func )


        to_be_replaced = creator_name.rpartition("::")[2] + "::mover_name"
        replaced_with  = widget_name.rpartition("::")[2] + "::mover_name"
        for ii,line in enumerate(widget_cc_lines) :
            if len(line) >= 2 and line[0:2] == "//" : continue
            widget_cc_lines[ ii ] = line.replace( to_be_replaced, replaced_with )

        to_be_replaced = "class_name"
        replaced_with  = "mover_name"
        for ii,line in enumerate(widget_cc_lines) :
            if len(line) >= 2 and line[0:2] == "//" : continue
            widget_cc_lines[ ii ] = line.replace( to_be_replaced, replaced_with )

        creator_hh_insertion_line = -1
        widget_cc_insertion_line = -1
        widget_hh_insertion_line = -1

        for func in creator_hh_funcs :
            #print "creator func: ", func.scope, func.name
            if func.scope == creator_name :
                creator_hh_insertion_line = func.line_end
        for func in widget_cc_funcs :
            #print "cc func: ", func.scope, func.name
            if func.scope == creator_name or func.scope == widget_name :
                widget_cc_insertion_line = func.line_end
        for func in widget_hh_funcs :
            if func.scope == widget_name and func.privacy == "public" :
                widget_hh_insertion_line = func.line_end

        assert( creator_hh_insertion_line != -1 )
        assert( widget_cc_insertion_line != -1 )
        assert( widget_hh_insertion_line != -1 )

        for func in widget_cc_funcs_to_remove :
            widget_cc_lines = comment_out_function( func, widget_cc_lines )
        for func in widget_hh_funcs_to_remove :
            widget_hh_lines = comment_out_function( func, widget_hh_lines )
        for func in creator_hh_funcs_to_remove :
            creator_hh_lines = comment_out_function( func, creator_hh_lines )
        if creator_cc_funcs :
            for func in creator_cc_funcs_to_remove :
                creator_cc_lines = comment_out_function( func, creator_cc_lines )

        cc_replace_dict = {
            "classname" : widget_name.rpartition( "::" )[ 2 ],
            "defaultschemafunc" : default_helper_func,
            "classkey" : class_key,
            "creatorname" : creator_name.rpartition( "::" )[ 2 ] }
        cc_newlines = [ line % cc_replace_dict for line in cc_stub() ]
        widget_cc_lines = widget_cc_lines[:(widget_cc_insertion_line+1)] + [ "\n" ] + cc_newlines + [ "\n" ] + widget_cc_lines[(widget_cc_insertion_line+1):]

        creator_hh_lines = creator_hh_lines[:(creator_hh_insertion_line+1)] + creator_hh_stub() + creator_hh_lines[(creator_hh_insertion_line+1):]

        widget_hh_newlines = [ line % cc_replace_dict for line in class_hh_stub() ]
        widget_hh_lines = widget_hh_lines[:(widget_hh_insertion_line+1)] + ["\n"] + widget_hh_newlines + ["\n"] + widget_hh_lines[(widget_hh_insertion_line+1):]

        cc_includes = ["// XSD XRW Includes\n", "#include <utility/tag/XMLSchemaGeneration.hh>\n"]
        for inc in extra_cc_includes :
            cc_includes.append( "#include <%s>\n" % inc )
        if creator_cc_lines :
            cc_includes.append( "#include <%s>\n" % ( creator_hh_fname.partition( "src/" )[2] ) )
        widget_cc_lines = add_include_at_bottom_of_includes( widget_cc_lines, cc_includes )

        open( creator_hh_fname, "w" ).writelines( creator_hh_lines )
        open( widget_hh_fname, "w" ).writelines( widget_hh_lines )
        open( widget_cc_fname, "w" ).writelines( widget_cc_lines )
        if creator_cc_lines : open( creator_cc_fname, "w" ).writelines( creator_cc_lines )
        
