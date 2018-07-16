#!/usr/bin/python
from __future__ import print_function
import sys
import blargs

# This file defines the Beautifier class which parses Rosetta into an AST
# and then determines the indentation for each line based on the AST.
#
# This file also can be also be run as a stand-alone script to beautify
# a single source file.  The following two options are useful:
# --filename <fname> : which file should be beautified
# --overwrite : overwrite the input file with the beautified version of the file.
#               This triggers additional safety checks to ensure the logical
#               structure of the beautified file is identical to that of the original.
#               In the absence of this flag, the beautified output of X.cc is written
#               to X.cc.beaut


debug = False
#debug = True

debug_equiv = False
#debug_equiv = True

token_types = [ "top-level",
                "namespace",
                "namespace-scope",
                "scope",
                "for",
                "for-declaration",
                "for-scope",
                "do-while",
                "do-while-condition",
                "do-while-scope",
                "if",
                "if-condition",
                "if-scope",
                "else",
                "else-scope",
                "while",
                "while-condition",
                "while-scope",
                "switch",
                "switch-expression",
                "switch-scope",
                "case",
                "case-block",
                "case-block-scope",
                "class",
                "class-inheritance-list",
                "class-scope",
                "class-privacy",
                "struct",
                "struct-inheritance-list",
                "struct-scope",
                "union",
                "union-scope",
                "function",
                "function-decl-argument-list",
                "function-scope",
                "try",
                "try-scope",
                "catch",
                "catch-arg",
                "catch-scope",
                "enum",
                "enum-scope",
                "ctor-initializer-list",
                "block-comment",
                "statement",
                "statement-scope", # new for lambda functions that are defined within statements
                "substatement",
                "template",
                "template-arg-list",
                "empty-statement" ]

white_listed_macros = set([
    "OPT_KEY",
    "OPT_1GRP_KEY",
    "OPT_2GRP_KEY",
    "OPT_3GRP_KEY",
    "EXTERN_OPT_KEY",
    "EXTERN_OPT_1GRP_KEY",
    "EXTERN_OPT_2GRP_KEY",
    "EXTERN_OPT_3GRP_KEY",
    "ASSERT_ONLY",
    "CEREAL_FORCE_DYNAMIC_INIT",
    "CEREAL_REGISTER_DYNAMIC_INIT",
    "CEREAL_REGISTER_TYPE",
    "SPECIAL_COP_SERIALIZATION_HANDLING",
])

class Token :
    def __init__( self ) :
        self.spelling = ""
        self.start = 0
        self.one_past_end = 0
        self.line_number = 0
        self.orig_line_number = 0
        self.parent = None
        self.children = []
        self.is_reserve = False
        self.is_commented = False
        self.is_visible = True
        self.is_inside_string = False
        self.is_preprocessor_directive = False
        self.invisible_by_macro = False
        self.is_macro = False
        self.type = "substatement" # default is a vanilla "I'm part of another statement" type.
        self.depth = 0
        self.index = -1
    def context( self ) :
        if not self.parent :
            return "namespace-scope"
        else :
            return self.parent.type
    def set_parent( self, parent ) :
        self.parent = parent
        parent.children.append( self )
    def equivalent( self, other ) :
        if self.spelling != other.spelling :
            if not ( self.spelling == "\t" or self.spelling == "\\t" and other.spelling == "\t" or self.spelling == "\\t" ) :
                return False
        elif self.type != other.type : return False
        elif self.invisible_by_macro != other.invisible_by_macro : return False
        elif self.is_visible != other.is_visible : return False
        elif self.is_inside_string != other.is_inside_string : return False
        return True
    def parents_last_child( self ) :
        if not self.parent : return False
        return self is self.parent.children[-1]
    def regular_rcb( self ) :
        # is this a non-string, visible, left-curly brace?
        return not self.is_inside_string and self.is_visible and self.spelling == "}"


class Beautifier :
    def __init__( self ) :
        self.filename = ""
        self.line = ""
        self.all_lines = [] # this doesn't change after the file is read
        self.new_lines = []
        self.line_number = -1 # so that the first line ends up line 0
        self.this_line_tokens = []
        self.all_tokens = [] # tokens for the entire file
        self.line_tokens = [] # the tokens for each line
        self.new_line_tokens = [] # the set of tokens for each line as the lines are getting rewritten
        self.namespace_stack = []
        self.line_indentations = []
        self.context_stack = []
        self.nested_ifdefs  = [ ( "None", "None", True ) ]
        self.defined_macros = []
        self.macros_inquired_about = []
        self.reserve_words = set( [ "using", "namespace", "class", "for", "while", "do", \
            "repeat", "public", "private", "protected", "template", "typedef", "typename", \
            "operator", "inline", "explicit", "static", "mutable", "virtual", "friend", \
            "unsigned", "struct", "union", "try", "catch", "override", "default", "delete" ] )
        self.macros_at_zero_indentation = set( [ "#if", "#ifdef", "#ifndef", "#else", "#elif", "#endif" ] )
        self.whitespace = set([" ", "\t", "\n"])
        self.dividers = set([";",":",",","(",")","{","}","=","[","]","<",">","&","|","\\",
                             '"',"'","?","!","+","-"]) # *, and / are not dividers, but are treated like one
        self.privacy_types = set([ "public", "protected", "private"])
        self.scope_types = set( ["namespace-scope", "for-scope", "do-while-scope", "if-scope",
                                 "else-scope", "while-scope", "switch-scope",
                                 "case-block-scope", "class-scope", "struct-scope", "union-scope",
                                 "function-scope", "try-scope", "catch-scope", "scope", "statement-scope" ] )
        self.for_types = set( [ "for", "BOOST_FOREACH", "foreach", "foreach_", "boost_foreach", "FEM_DO",
                                "FEM_DO_SAFE", "FEM_DOSTEP",
                                "FORVC", "FORTAGS" ] )
        self.pound_if_setting = "" # can be "take_if" or "take_else" to broadly specify to take one or the other for a file
                                   # without providing a full-fledged preprocessor directive reader, or "take_neither" to skip everything
        self.binary_math_symbols = set(["*","-","+","/"])
        self.in_comment_block = False
        self.in_string = False
        self.in_char_quote = False
        self.in_macro_definition = False

    def macro_definitions_say_line_visible( self ) :
        return self.nested_ifdefs[-1][2]

    # look at the contents of a line and chop it up into a set of tokens
    # then decide what tokens are parts of strings and what parts are
    # parts of comments and what parts are invisible due to macro definitions.
    # Save all the tokens seen in all lines so that they can be "parsed"
    # and understood.
    def tokenize_line( self, line ) :
        self.line_number += 1
        self.breakup_line_into_tokens( line )
        self.process_strings_and_comments()

        # copy the list, shallow-copy the list elements
        self.line_tokens.append( list( self.this_line_tokens ) )

        # steal the elements of the list and add them to the
        # set of all the tokens for the file
        self.all_tokens.extend( self.this_line_tokens )


    def breakup_line_into_tokens( self, line ) :
        self.line = line
        self.all_lines.append( line )
        self.this_line_tokens = []
        start = -1
        i = 0
        while ( i < len(self.line)) :
            #print(len(self.all_lines), i, self.line[i], start)
            if self.line[i] in self.whitespace :
                if start >= 0 :
                    self.take_token_for_range( start, i )
                start = -1
                if self.line[i] == "\t" :
                    self.take_token_for_range(i,i+1)
                i+=1
                continue
            if self.line[ i ] in self.dividers :
                if start >= 0 :
                    self.take_token_for_range( start, i )
                if i+1 < len(self.line) and self.line[i] == "\\" and self.line[i+1] == "t" :
                    self.take_token_for_range( i, i+2 )
                    start = -1
                    i+=2
                elif i+1 < len(self.line) and ( self.line[i] == "<" or self.line[i] == ">" ) and self.line[i+1] == self.line[i] :
                    self.take_token_for_range( i, i+2 )
                    start = -1
                    i+=2
                else :
                    self.take_token_for_range( i, i+1 )
                    start = -1
                    i+=1
                continue
            # make to-end-of-line comments and block-comment starting symbols their own tokens
            if self.line[ i ] == "/" :
                if i+1 < len(self.line ) and ( self.line[i+1] == "/" or self.line[i+1] == "*" ) :
                    if start >= 0  :
                        self.take_token_for_range( start, i )
                    self.take_token_for_range( i, i+2 )
                    start = -1
                    i+=2
                    continue
                else :
                    if start >= 0 :
                        self.take_token_for_range( start, i )
                    self.take_token_for_range(i,i+1)
                    start = -1
                    i+=1
                    continue
            # make block-comment ending symbols their own tokens
            if self.line[ i ] == "*" :
                if i+1 < len( self.line ) and self.line[i+1] == "/" :
                    if start >= 0 :
                        self.take_token_for_range( start, i )
                    self.take_token_for_range( i, i+2 )
                    start = -1
                    i+=2
                    continue
                if start >= 0 :
                    self.take_token_for_range( start, i )
                self.take_token_for_range(i,i+1)
                start = -1
                i+=1
                continue
            # make string opening or closing symbols their own tokens
            if self.line[ i ] == '"' :
                if i > 0 and ( self.line[i-1] != "\\" ) :
                    if start >= 0 :
                        self.take_token_for_range( start, i )
                    self.take_token_for_range( i,i+1)
                    start = -1
                    i += 1
                    continue
            # if its not whitespace
            if start < 0 :
                start = i
            i+=1
        if start >= 0 :
            self.take_token_for_range( start, len(self.line) )

    def take_token_for_range( self, start, end_plus_one ) :
        tok = Token()
        tok.spelling = self.line[start:end_plus_one]
        tok.start = start
        tok.one_past_end = end_plus_one
        tok.line_number = self.line_number
        tok.orig_line_number = self.line_number
        self.this_line_tokens.append( tok )

    def process_strings_and_comments( self ) :
        # first check to see if we're in a comment block
        # and then check to see if the first token begins with a #, if it does, then
        # the whole line is invisible, and we have to decide if the next lines
        # are visible or not

        process_as_macro_line = False
        for tok in self.this_line_tokens :
            if tok.spelling == "\t" : continue
            if not self.in_comment_block and tok.spelling[0] == "#" :
                # print("process as macro line", self.line.strip())
                process_as_macro_line = True
            else :
                break;

        if self.in_macro_definition :
            # just don't mess with macro definitions; they shouldn't even be there!
            if debug: print("still in macro definition", self.line_number)
            if self.this_line_tokens :
                for tok in self.this_line_tokens :
                    tok.is_visible = False
                    tok.invisible_by_macro = True
                assert( len( self.this_line_tokens ) != 0 )
                if self.this_line_tokens[-1].spelling != "\\" :
                    self.in_macro_definition = False
                    if debug: print("final line in macro definition", self.line_number)
                return
            else :
                # empty line
                self.in_macro_definition = False
                if debug: print("final line in macro definition", self.line_number)
                return

        escape_tok = None
        single_char = False
        i = 0
        while i < len( self.this_line_tokens ) :
            # for i, tok in enumerate( self.this_line_tokens ) :
            tok = self.this_line_tokens[i]
            if not self.macro_definitions_say_line_visible() :
                tok.is_visible = False
                tok.invisible_by_macro = True
                # print("token is invisible by macro:", tok.spelling, tok.line_number+1)
            if not self.in_string and tok.spelling == "\t" :
                # keep around tab tokens only so long as they are inside strings
                # they'll be replaced by escaped tabs when the file is written out
                self.this_line_tokens.remove( tok )
                continue; # no increment
            elif self.in_comment_block :
                tok.is_visible = False
                tok.is_commented = True
                if tok.spelling == "*/" and self.macro_definitions_say_line_visible() :
                    self.in_comment_block = False
            elif self.in_string :
                if debug : print("in string", tok.spelling, tok.line_number)
                tok.is_inside_string = True
                if tok.spelling == "\\" :
                    if escape_tok and tok.start == escape_tok.one_past_end :
                        escape_tok = None
                    else :
                        escape_tok = tok
                elif tok.spelling == '"':
                    if escape_tok and tok.start == escape_tok.one_past_end :
                        escape_tok = None
                    else :
                        self.in_string = False
                else :
                    escape_tok = None
            elif self.in_char_quote :
                if debug : print("in char quote", tok.spelling, tok.line_number)
                tok.is_inside_string = True
                if tok.spelling == "\\" :
                    if escape_tok and tok.start == escape_tok.one_past_end :
                        escape_tok = None
                    else :
                        escape_tok = tok
                elif tok.spelling == "'":
                    if escape_tok and tok.start == escape_tok.one_past_end :
                        escape_tok = None
                    else :
                        self.in_char_quote = False
                else :
                    escape_tok = None
            elif len(tok.spelling) == 2 and ( tok.spelling == "//" or ( tok.spelling == "/*" and self.macro_definitions_say_line_visible() ) ) :
                # comment begin
                tok.is_visible = False
                tok.is_commented = True
                if tok.spelling == "//" :
                    #the rest of the line is a comment
                    j = i
                    while j < len(self.this_line_tokens) :
                        if self.this_line_tokens[j].spelling == "\t" :
                            self.this_line_tokens.remove(self.this_line_tokens[j])
                        else :
                            self.this_line_tokens[j].is_visible = False
                            self.this_line_tokens[j].is_commented = True
                            j+=1
                    break
                if tok.spelling == "/*" :
                    #we're in a comment block
                    self.in_comment_block = True
            elif tok.spelling == '"' :
                # double check we're not inside the single character " held inside two single quotes
                # should now be unreachable if i > 0 and self.this_line_tokens[i-1].spelling == "'" :
                # should now be unreachable     i+=1
                # should now be unreachable     continue
                if debug : print("in string", tok.spelling, tok.line_number)
                tok.is_inside_string = True
                self.in_string = True
            elif tok.spelling == "'" :
                if debug : print("in char quote", tok.spelling, tok.line_number)
                tok.is_inside_string = True
                self.in_char_quote = True
            # elif tok.spelling == "'" and not self.in_string :
            #     if single_char :
            #         single_char = False
            #     else :
            #         single_char = True
            #     if debug : print("in string", tok.spelling, tok.line_number)
            #     tok.is_inside_string = True
            # elif single_char :
            #     if debug : print("in string", tok.spelling, tok.line_number)
            #     tok.is_inside_string = True
            i += 1

        if process_as_macro_line : self.process_macro_line()

    def process_macro_line( self ) :
        # this line has been shown to start with "#" so the entire line is being treated as a macro definition
        # print("process macro line: ", self.all_lines[-1].strip())

        for tok in self.this_line_tokens :
            tok.is_visible = False
            tok.is_preprocessor_directive = True
            # print("tok: ", tok.spelling, "is now invisible")
        tok0 = self.this_line_tokens[ 0 ]
        if tok0.spelling == "#include" or tok0.spelling == "#pragma" :
            return

        for tok in self.this_line_tokens :
            if tok.spelling not in self.macros_inquired_about :
                self.macros_inquired_about.append( tok.spelling )

        if tok0.spelling in self.macros_at_zero_indentation :
            # treat certain macros that may be invisible as if they are visible for
            # the sake of indenting them to level 0.
            for tok in self.this_line_tokens :
                tok.invisible_by_macro = False
        if tok0.spelling == "#define" :
            self.defined_macros.append( self.this_line_tokens[1].spelling )
            if debug: print("#defined:", self.this_line_tokens[1].spelling)
            if self.this_line_tokens[-1].spelling == "\\" :
                self.in_macro_definition = True
                if debug: print("in macro definition", self.line_number)
        elif tok0.spelling == "#ifdef" :
            if len(self.this_line_tokens) > 1 :
                name = self.this_line_tokens[1].spelling
                if name in self.defined_macros :
                    # say "true" if the name is in the set of defined macros, and
                    # the parent macro visibility was also "true"
                    self.nested_ifdefs.append( (name, "ifdef", self.nested_ifdefs[-1][2] ) )
                else :
                    self.nested_ifdefs.append( (name, "ifdef", False ))
        elif tok0.spelling == "#ifndef" :
            if len(self.this_line_tokens) > 1 :
                name = self.this_line_tokens[1].spelling
                if name not in self.defined_macros :
                    self.nested_ifdefs.append( (name, "ifndef", self.nested_ifdefs[-1][2] ) )
                else :
                    self.nested_ifdefs.append( (name, "ifndef", False ))
        elif tok0.spelling == "#if" :
            name = " ".join( [ x.spelling for x in self.this_line_tokens[1:] if not x.is_commented ] )
            #print("hit #if", self.pound_if_setting, self.line_number)
            if self.pound_if_setting == "take_if" :
                #print("Take if, arived at #if:", self.all_lines[-1].strip())
                self.nested_ifdefs.append( (name, "if", self.nested_ifdefs[-1][2] ) )
            elif self.pound_if_setting == "take_else" :
                #print("Take else, arived at #if:", self.all_lines[-1].strip())
                self.nested_ifdefs.append( (name, "if", False ) )
            elif self.pound_if_setting == "take_neither" :
                self.nested_ifdefs.append( (name, "if", False ) )
            else :
                # conservative -- don't try to parse or logically evaluate this line
                # just say "yeah, it's probably true" and if you hit an else later, then
                # instead of negating it, say the else is also true.
                # this could cause problems
                #print("take both")
                self.nested_ifdefs.append( (name, "if", self.nested_ifdefs[-1][2] ))
        elif tok0.spelling == "#else" :
            name, iftype, truthval = self.nested_ifdefs.pop()
            #print("pop!", name, iftype, truthval)
            if iftype == "else" :
                # we have hit a problem
                print("ERROR: Encountered '#else' inappropriately.  Ifdef stack: ")
                for n2, it2, tv2 in reversed( self.nested_ifdefs ) :
                    print("   ", n2, it2, tv2)
                assert( iftype != "else" )
            if self.nested_ifdefs[-1][2] :
                if iftype == "ifdef" or iftype == "ifndef" :
                    self.nested_ifdefs.append( (name, "else", not truthval) )
                else :
                    assert( iftype == "if" )
                    if self.pound_if_setting == "take_if" :
                        #print("Take if, arived at #else:", self.all_lines[-1].strip())
                        self.nested_ifdefs.append( (name, "else", False ) )
                    elif self.pound_if_setting == "take_else" :
                        #print("Take else, arived at #else:", self.all_lines[-1].strip())
                        self.nested_ifdefs.append( (name, "else", True ) )
                    elif self.pound_if_setting == "take_neither" :
                        self.nested_ifdefs.append( (name, "else", False ) )
                    else :
                        #print("conservative -- take both")
                        self.nested_ifdefs.append( (name, "else", True) )
            else :
                self.nested_ifdefs.append( (name, "else", False ))
        elif tok0.spelling == "#endif" :
            assert( self.nested_ifdefs[-1][1] in set( ["if", "ifdef", "ifndef", "else" ]) )
            self.nested_ifdefs.pop()
        else :
            pass
            #print("WARNING: Unhandled preprocessor directive:", tok0.spelling, "on", self.line_number+1, "of", self.filename, ":", self.line,)
        #print("processing macros", self.line_number+1, ", ".join( [ "%s %s" % ( x[1], "true" if x[2] else "false" ) for x in self.nested_ifdefs ] ))


    def find_next_visible_token( self, i, stack=None ) :
        if debug : self.print_entry( "find_next_visible_token", i, stack )
        while i < len(self.all_tokens) :
            if self.all_tokens[i].is_visible : return i
            if stack != None :
                self.set_parent( i, stack )
            i += 1
        return i

    def print_entry( self, fname, i, stack ) :
        if stack :
            if i < len(self.all_tokens) :
                print((" "*len(stack)), fname, i, self.all_tokens[i].spelling, "line number", self.all_tokens[i].line_number+1, self.all_tokens[i].is_visible)
            else :
                print((" "*len(stack)),fname,i,"reached last token")
        else :
            if i < len(self.all_tokens) :
                # print("stack-less", fname, "token", i, ",", self.all_tokens[i].spelling, "line number", self.all_tokens[i].line_number+1)
                pass
            else :
                # print("stack-less", fname, "token", i, "beyond last token")
                pass

    def renumber_tokens( self ) :
        for i,tok in enumerate(self.all_tokens) :
            tok.index = i

    def minimally_parse_file( self ) :
        # at this point, the file has been split into a large set of tokens
        # and we're going to go through the tokens and, to a limited extent,
        # interpret the tokens into a structure
        self.renumber_tokens()
        i = 0
        self.toplevel_token = Token()
        self.toplevel_token.type = "namespace-scope"
        stack = [ self.toplevel_token ] # you start at namespace scope
        while i < len( self.all_tokens ) :
            i = self.process_statement( i, stack )

        for i in xrange(len(self.all_tokens)) :
            if self.all_tokens[ i ].parent is self.all_tokens[i] :
                itok = self.all_tokens[i]
                print("Token", i, "is it's own parent:", itok.spelling, itok.line_number, itok.type)
                assert( not self.all_tokens[ i ].parent is self.all_tokens[i] )
            for tok in self.all_tokens[ i ].children :
                assert( not tok is self.all_tokens[i] )

        if debug : self.print_depth_stack( stack )



    def process_statement( self, i, stack ) :
        if debug : self.print_entry("process_statement",i,stack)
        i = self.find_next_visible_token(i,stack)
        if i == len( self.all_tokens ) : return i
        if debug : self.print_entry("process_statement",i,stack)
        i_spelling = self.all_tokens[i].spelling
        if i_spelling in self.for_types :
            return self.process_for(i,stack)
        elif i_spelling == "if" :
            return self.process_if(i,stack)
        elif i_spelling == "do" :
            return self.process_do_while(i,stack)
        elif i_spelling == "while" :
            return self.process_while(i,stack)
        elif i_spelling == "namespace" :
            return self.process_namespace(i,stack)
        elif i_spelling == "class" or i_spelling == "struct" :
            # treat classes and structs identically
            return self.process_class_decl(i,stack)
        elif i_spelling == "union" :
            return self.process_union(i,stack)
        elif i_spelling == "switch" :
            return self.process_switch(i,stack)
        elif i_spelling == "case" or i_spelling == "default" :
            return self.process_case(i,stack)
        elif i_spelling == "try" :
            return self.process_try(i,stack)
        elif i_spelling == "{" :
            return self.process_scope(i,stack)
        elif i_spelling == ";" :
            self.set_parent(i,stack,"empty-statement")
            return i+1;
        elif i_spelling == "using" :
            return self.process_simple_statement(i,stack)
        elif i_spelling == "typedef" :
            j = self.find_next_visible_token(i+1)
            #print("encountered typedef", self.all_tokens[j].spelling)
            if self.all_tokens[j].spelling == "struct" :
                self.set_parent(i,stack,"statement")
                i = self.find_next_visible_token(i+1,stack)
                return self.process_class_decl(i,stack)
            else :
                return self.process_simple_statement(i,stack)
        elif i_spelling == "public" or i_spelling == "protected" or i_spelling == "private" :
            return self.process_privacy_declaration(i,stack)
        elif i_spelling == "template" :
            # just ignore templates, treat them like statements that don't end in ;'s
            i = self.process_template(i,stack)
            return self.process_statement(i,stack)
        elif i_spelling == "enum" :
            return self.process_enum(i,stack)
        elif i_spelling in white_listed_macros :
            return self.process_simple_statement(i,stack)
        elif stack[-1].type == "namespace-scope" or stack[-1].type == "class-scope" :
            return self.process_function_preamble_or_variable(i,stack)
        else :
            return self.process_simple_statement(i,stack)

    def set_parent( self, i, stack, token_type = "substatement" ) :
        #if debug : print(" " * len(stack), " set parent of ", i, stack[-1].index, self.all_tokens[i].spelling, stack[-1].spelling)
        self.all_tokens[i].type   = token_type
        self.all_tokens[i].parent = stack[-1]
        stack[-1].children.append( self.all_tokens[i] )

    def handle_read_all_tokens_error( self, function_name, stack ) :
        print("ERROR: Ran out of tokens in function", function_name, "in file", self.filename)
        self.print_depth_stack( stack )
        sys.exit(1)

    def print_depth_stack( self, stack ) :
        print("Depth stack: ",)
        for elem in reversed( stack ) :
            if elem :
                print("  ", elem.spelling, (  "%d %s" % ( elem.line_number+1, self.all_lines[ elem.line_number ].strip() ) ), "vis" if elem.is_visible else "invis")

    def read_to_end_paren( self, i, stack ) :
        # the starting "(" should have already been read
        paren_depth = 1
        while i < len( self.all_tokens ) :
            if debug : print((" "*len(stack)), "read to end paren", i, self.all_tokens[i].spelling, self.all_tokens[i].line_number+1, self.all_tokens[i].is_visible)
            self.set_parent(i,stack)
            if not self.all_tokens[i].is_visible or self.all_tokens[i].is_inside_string :
                pass
            elif self.all_tokens[i].spelling == ")" :
                paren_depth -= 1
                if paren_depth == 0 :
                    return i+1
            elif self.all_tokens[i].spelling == "(" :
                paren_depth += 1
            i += 1
        self.handle_read_all_tokens_error("read_to_end_paren", stack)

    def process_template(self,i,stack) :
        if debug : self.print_entry("process_template",i,stack)
        self.set_parent(i,stack,"template")
        stack.append(self.all_tokens[i])
        i+=1
        template_params = "not-begun"
        template_depth = 0
        while i < len(self.all_tokens) :
            #print("process template ", i, self.all_tokens[i].spelling)
            if self.all_tokens[i].spelling == "<" :
                template_depth += 1
                if template_params == "not-begun" :
                    self.set_parent(i,stack,"template-arg-list")
                    stack.append(self.all_tokens[i])
                    template_params = "started"
                    i+=1
                    continue
            elif self.all_tokens[i].spelling == ">" :
                assert( template_params == "started" )
                template_depth -= 1
                if template_depth == 0 :
                    self.set_parent(i,stack)
                    stack.pop() #pop template-arg-list
                    stack.pop() #pop template
                    if debug : self.print_entry("exitting process_template",i,stack)
                    return i+1
            elif self.all_tokens[i].spelling == ";" and template_depth == 0 :
                # perhaps we'll never find the template params
                assert( template_params == "not-begun" )
                self.set_parent(i,stack)
                stack.pop() # pop template
                return i+1
            self.set_parent(i,stack)
            i+=1
        self.handle_read_all_tokens_error( "process_template", stack )

    def var_in_class_dec_initializer( self, stack ) :
        # print("stack[-1].type", stack[-1].type)
        # print("stack[-2].type", stack[-2].type)
        # print("stack[-3].type", stack[-3].type)
        # print("stack[-4].type", stack[-4].type)
        return len(stack) >= 4 and \
            stack[-1].type == "statement" and \
            stack[-2].type == "function-scope" and \
            stack[-3].type == "function" and \
            stack[-4].type == "class-scope"

    def check_lambda_compatability( self, token_list ) :
        if not self.could_be_lambda( token_list ) and len( token_list ) != 0 :
            del token_list[:]

    def could_be_lambda( self, token_list ) :
        if len(token_list) == 0 : return False
        seen_lparen = False
        last_was_minus = False
        seen_arrow = False
        seen_throw = False
        seen_throw_lparen = False
        template_depth = 0; # how many "<"s have we seen?
        for tok_ind in xrange( len(token_list) ) :
            tok = token_list[ tok_ind ]
            if tok_ind == 0 :
                if tok != "[" : return False
            elif tok_ind == 1 :
                if tok == "{" : return True
                if tok != "(" : return False
                seen_lparen = True
            elif last_was_minus and tok != ">" :
                return False
            elif last_was_minus and tok == ">" :
                if seen_arrow : return False
                seen_arrow = True
            elif tok == "(" :
                if seen_lparen :
                    if not seen_throw or ( seen_throw and seen_throw_lparen ) :
                        return False
                    seen_throw_lparen = True
                else :
                    seen_lparen = True
            elif tok == "-" :
                if seen_lparen and template_depth == 0 :
                    last_was_minus = True
                else :
                    return False
            elif seen_lparen and not seen_arrow :
                if tok == "throw" :
                    if seen_throw or seen_nothrow: return False
                    seen_throw = True
                elif tok == "nothrow" :
                    if seen_nothrow or seen_throw: return False
                    seen_nothrow = True
            #elif tok == "{" :
            #    # at this point; tok_ind should point to the last element in the list
            #    return True
            elif seen_arrow :
                # this could be a return type and I don't want to figure out if it is, so, let's assume it is
                # if we see a "<" treat it like the start of a template parameter
                if tok == "<" :
                    template_depth += 1
                if tok == ">" :
                    if template_depth > 0 : template_depth -= 1
                elif tok == "," :
                    if template_depth == 0 : return False
            else :
                # basically any other thing that we find prevents this from being a lambda
                return False
        return True

    def process_simple_statement(self,i,stack) :
        if debug : self.print_entry("process_simple_statement",i,stack)

        open_scope = "{"
        close_scope = "}"
        open_brackets = set( [ "[", "(", "{" ] )
        close_brackets = set( [ "]", ")", "}" ] )
        closing_bracket_for_opener = { "[" : "]", "(" : ")", "{" : "}" }
        bracket_stack = []
        lambda_detection_stack = [ [] ] # what symbols have we found consistent with a lambda function declaration given our current bracket stack?

        if self.all_tokens[i].spelling in open_brackets and not self.all_tokens[i].is_inside_string :
            bracket_stack.append( self.all_tokens[i].spelling )
            lambda_detection_stack[-1].append( self.all_tokens[i].spelling )
            lambda_detection_stack.append( [] )

        macro_without_ending_semicolon = False
        if self.all_tokens[i].spelling in white_listed_macros :
            if debug : print("Encountered white-listed macro", self.all_tokens[i].spelling)
            macro_without_ending_semicolon = True
            encountered_first_lparen = False

        self.set_parent(i,stack,"statement")
        stack.append(self.all_tokens[i])

        # global-scope function call prefix test
        # e.g. "::time()"
        if self.all_tokens[i].spelling == ":" and i+1 < len(self.all_tokens) and self.all_tokens[i+1].spelling == ":" :
            self.set_parent(i+1,stack)
            i+=2
        else :
            i+=1;

        # read all the tokes from i until the next semicolon
        last = i
        n_question_marks = 0
        while i < len(self.all_tokens) :

            if not self.all_tokens[i].is_visible or self.all_tokens[i].is_inside_string :
                self.set_parent(i,stack)
                i+=1
                continue

            lambda_detection_stack[-1].append( self.all_tokens[i].spelling )
            #print(self.all_tokens[i].spelling, len( lambda_detection_stack ))
            self.check_lambda_compatability( lambda_detection_stack[-1] )

            if self.all_tokens[i].spelling == ";" and len(bracket_stack) == 0 :
                # print("statement:", " ".join( [ x.spelling for x in self.all_tokens[i:i+1] ] ))
                self.set_parent(i,stack)
                stack.pop()
                if debug : self.print_entry("exitting process_simple_statement",i,stack)
                return i+1
            elif self.all_tokens[i].spelling == ":" :
                if self.all_tokens[i+1].spelling != ":" :
                    if debug : print("found colon:", i, self.all_tokens[i].line_number, self.all_tokens[i].is_inside_string, self.all_tokens[i].start, "n?s", n_question_marks, "i+1", self.all_tokens[i+1].spelling)
                    if n_question_marks == 0 :
                        # goto handling; you're allowed to define statements on a line and end them with a colon
                        # I hope I'm not missing anything important here
                        self.set_parent(i,stack)
                        stack.pop()
                        i = self.find_next_visible_token(i+1,stack)
                        if debug : self.print_entry("label declaration; soldiering onwards",i,stack)
                        return self.process_statement(i,stack)
                    else :
                        n_question_marks -= 1;
                else :
                    self.set_parent(i,stack)
                    self.set_parent(i+1,stack)
                    i+=2
                    continue
            elif self.all_tokens[i].spelling == "?" :
                n_question_marks += 1
            elif self.all_tokens[i].spelling == open_scope and self.could_be_lambda( lambda_detection_stack[-1] ) :
                if debug : self.print_entry( "found open scope inside simple statement", i, stack )
                i = self.process_scope(i,stack)
                if debug : self.print_entry( "finished dealing with open scope inside simple statement", i, stack )
                del lambda_detection_stack[-1][:] # empty out the list of the last entry in the lambda-detection stack
                continue
            elif self.all_tokens[i].spelling in open_brackets :
                bracket_stack.append( self.all_tokens[i].spelling )
                lambda_detection_stack.append( [ ] )
                if macro_without_ending_semicolon and not encountered_first_lparen and self.all_tokens[i].spelling == "(" :
                    encountered_first_lparen = True
            elif self.all_tokens[i].spelling in close_brackets :
                if len(bracket_stack) == 0 :
                    # Perhaps this is a constant declared in a class?
                    if self.var_in_class_dec_initializer( stack ) and self.all_tokens[i].spelling == "}" :
                        if debug : self.print_entry("simple statement has turned out to be a variable!",i,stack )
                        # OK! -- then we found what we were looking for and can quit
                        self.set_parent(i,stack)
                        stack.pop()
                        return i
                    print("Closing bracket '"+self.all_tokens[i].spelling+"' found that does not match an opening bracket in", self.filename)
                    sys.exit(1)
                if closing_bracket_for_opener[ bracket_stack[-1] ] != self.all_tokens[i].spelling :
                    print("Closing bracket does not match opening bracket in", self.filename)
                    print("opener:", bracket_stack[-1])
                    print("closer:", self.all_tokens[i].spelling)
                    print("expected:", closing_bracket_for_opener[ bracket_stack[-1] ])
                    #assert( closing_bracket_for_opener[ bracket_stack[-1] ] == self.all_tokens[i].spelling  )
                    sys.exit(1)
                bracket_stack.pop()
                lambda_detection_stack.pop()
                if len(bracket_stack) == 0 and macro_without_ending_semicolon and encountered_first_lparen :
                    self.set_parent(i,stack)
                    stack.pop()
                    return i+1

            self.set_parent(i,stack)
            i+=1
        self.handle_read_all_tokens_error( "process_simple_statement", stack )

    def current_classname( self, i, stack ) :
        assert( stack[-1].type == "class-scope" )
        assert( len(stack) > 2 )
        # let's get the parent "class" token
        class_token = stack[-2]
        assert( class_token )
        # and then the next visible token should be the class name
        # class names are never namespace qualified are they? That'd make that annoying
        i = self.find_next_visible_token( class_token.index+1 )
        return self.all_tokens[i].spelling

    def function_is_constructor( self, functionname, classname ) :
        #print("function is constructor? ", functionname, classname)
        return functionname == classname

    def process_function_preamble_or_variable(self,i,stack) :
        # print("process_function_preamble_or_variable", i, depth)
        # treat function declarations and variable declarations interchangably.
        # so we should gobble up the tokens that make this a function
        # and see if we have a constructor, in which case, we should also
        # be watchful for an intializer list
        if debug : self.print_entry("process_function_preamble_or_variable",i,stack)
        classname = ""
        if stack[-1].type == "class-scope" :
            classname = self.current_classname( i, stack )
            #if debug: print("class name!", classname)
        arglist = "not-begun"
        is_ctor = False
        found_init_list = False
        found_square_brackets = False
        found_parens = False
        found_operator = False
        found_equals = False

        i_initial = i # save this in case we aren't actually entering a function; treat it like a statement in that case
        self.set_parent(i,stack,"function")
        stack.append( self.all_tokens[i] )
        i+=1
        while i < len(self.all_tokens) :
            # print((" "*len(stack)), "process function preamble or variable", i, self.all_tokens[i].spelling, self.all_tokens[i].line_number+1)
            if not self.all_tokens[i].is_visible or self.all_tokens[i].is_inside_string :
                # don't do anything special; just don't try and interpret this token as having meaning
                pass
            elif self.all_tokens[i].spelling == "[" :
                # ok -- maybe we're looking at an array declaration
                if not found_parens :
                    found_square_brackets = True
            elif self.all_tokens[i].spelling == "=" :
                # ok -- maybe we're looking at a variable declaration w/ assignment
                if not found_parens :
                    found_equals = True
            elif self.all_tokens[i].spelling == ":" :
                # print("found : -- ", classname, is_ctor, arglist, found_init_list)
                if is_ctor and arglist == "ended" and not found_init_list :
                    # we're now reading the initializer list
                    # stack.pop() # remove the function-decl-argument-list
                    self.set_parent( i,stack,"ctor-initializer-list")
                    stack.append(self.all_tokens[i])
                    i+=1
                    found_init_list = True
                    continue
                elif arglist == "not-begun" :
                    # this might be constructor, we don't know
                    # this might also be some namespace qualifier on a return type
                    if i+1 < len(self.all_tokens) and self.all_tokens[i+1].spelling == ":" :
                        if self.all_tokens[i-1].spelling == ">" :
                            # The class is templated, so we have to look backwards past the template arguments
                            # to figure out its name so we can figure out if this is a constructor definition.
                            #print("found a '<'; searching backwards for a '<'")
                            j = i-2
                            count_gts = 1
                            while j > 0 :
                                #print("looking at", j, self.all_tokens[j].spelling, count_gts)
                                if not self.all_tokens[j].is_visible :
                                    pass
                                elif self.all_tokens[j].spelling == "<" :
                                    count_gts -= 1
                                    if count_gts == 0 :
                                        break
                                elif self.all_tokens[j].spelling == ">" :
                                    count_gts += 1
                                j -= 1
                            classname = self.all_tokens[j-1].spelling
                        else :
                            classname = self.all_tokens[i-1].spelling
                        #print("classname:", classname)

            elif self.all_tokens[i].spelling == "(" and arglist == "not-begun" :
                found_parens = True
                funcname = self.all_tokens[i-1].spelling
                if funcname != "operator" :
                    if self.function_is_constructor( funcname, classname ) :
                        is_ctor = True
                    arglist = "started"
                    self.set_parent(i,stack,"function-decl-argument-list")
                    stack.append(self.all_tokens[i])
                    i+=1
                    i = self.read_to_end_paren( i, stack )
                    arglist = "ended"
                    stack.pop() # remove decl-argument-list from stack
                    continue
            elif self.all_tokens[i].spelling == ";" :
                # ok, there's no function body here, just a function name declaration
                # so we're not going to increase the stack depth
                # in fact, this may not have been a function at all -- it might have
                # been a variale!
                assert( not found_init_list )
                self.set_parent(i,stack)
                stack.pop(); # get rid of "function" on stack
                i+=1
                if arglist == "not-begun" :
                    # ok this is just a variable, go ahead and reinterpret what we've just processed
                    # as a simple statement and remove my place in my parent's list of its children
                    token_i = self.all_tokens[i_initial].parent.children.pop()
                    assert( token_i is self.all_tokens[i_initial] )
                    i = self.process_simple_statement(i_initial,stack)
                if debug : self.print_entry("exitting process_function_preamble_or_variable",i,stack)
                return i

            elif self.all_tokens[i].spelling == "{" :
                # print("found_square_brackets", found_square_brackets)
                # print("found_equals", found_equals)
                # print("found_operator", found_operator)
                if ( found_square_brackets or found_equals ) and not found_operator :
                    i = self.process_to_end_rcb(i,stack,"array-value-initializer")
                    i = self.find_next_visible_token(i,stack)
                    assert( self.all_tokens[i].spelling == ";" )
                    self.set_parent(i,stack)
                    stack.pop() # remove function
                    return i+1
                else :
                    if ( found_init_list ) :
                        stack.pop() # remove ctor-initializer-list
                    # ok! now descend into the function block
                    i = self.process_statement( i, stack )
                    stack.pop()
                    if debug : self.print_entry("exitting process_function_preamble_or_variable",i,stack)
                    return i
            elif self.all_tokens[i].spelling == "operator" :
                found_operator = True

            self.set_parent(i,stack)
            i+=1
        self.handle_read_all_tokens_error( "process_function_preamble_or_variable", stack )

    def process_for(self,i,stack) :
        if debug : self.print_entry("process_for",i,stack)

        assert( self.all_tokens[i].spelling in self.for_types )
        self.set_parent(i,stack,"for")
        stack.append(self.all_tokens[i] )
        i+=1
        arglist = "not-yet-begun"
        while i < len(self.all_tokens) :
            if not self.all_tokens[i].is_visible or self.all_tokens[i].is_inside_string :
                pass
            elif self.all_tokens[i].spelling == "(" :
                if arglist == "not-yet-begun" :
                    arglist = "started"
                    self.set_parent(i,stack,"for-declaration")
                    stack.append(self.all_tokens[i])
                    i+=1
                    i = self.read_to_end_paren( i, stack )
                    # now remove for-declaration from the stack and process the next statment
                    stack.pop()
                    i = self.process_statement( i, stack )
                    stack.pop()
                    if debug : self.print_entry("exitting process_for",i,stack)
                    return i
            self.set_parent(i,stack)
            i += 1
        # whoa! we shouldn't have reached here!
        self.handle_read_all_tokens_error( "process_for", stack )

    def process_if(self,i,stack) :
        assert( self.all_tokens[i].spelling == "if" )
        if debug : self.print_entry("process_if",i,stack)
        self.set_parent(i,stack,"if")
        stack.append(self.all_tokens[i])
        i+=1
        while i < len(self.all_tokens) :
            if not self.all_tokens[i].is_visible or self.all_tokens[i].is_inside_string :
                pass
            elif self.all_tokens[i].spelling == "(" :
                self.set_parent(i,stack,"if-condition")
                stack.append(self.all_tokens[i])
                i+=1
                i = self.read_to_end_paren(i,stack)
                stack.pop() #remove if-condition
                i = self.process_statement( i, stack )
                j = self.find_next_visible_token( i )
                if self.all_tokens[j].spelling == "else" :
                    i = self.find_next_visible_token(i,stack)
                    self.set_parent(i,stack,"else") # set the parent for the else as the if
                    stack.append(self.all_tokens[i])
                    i = self.process_statement(i+1,stack)
                    stack.pop() #remove else
                stack.pop() #remove if
                if debug : self.print_entry("exitting process_if",i,stack)
                return i
            self.set_parent(i,stack)
            i += 1

        # whoa! We did not expect to get here!
        self.handle_read_all_tokens_error( "process_if", stack )

    def process_do_while( self, i, stack ) :
        assert( self.all_tokens[i].spelling == "do" )
        if debug : self.print_entry("process_do_while",i,stack)

        self.set_parent(i,stack,"do-while")
        stack.append(self.all_tokens[i])
        i+=1
        i = self.process_statement(i,stack)
        while i < len(self.all_tokens ) :
            if not self.all_tokens[i].is_visible or self.all_tokens[i].is_inside_string:
                pass
            elif self.all_tokens[i].spelling == "(" :
                self.set_parent(i,stack,"do-while-condition")
                stack.append(self.all_tokens[i])
                i+=1
                i = self.read_to_end_paren( i, stack )
                stack.pop() # pop do-while-condition
                i = self.find_next_visible_token(i,stack)
                assert( self.all_tokens[i].spelling == ";" )
                self.set_parent(i,stack)
                stack.pop() # pop do-while
                if debug : self.print_entry("exitting process_do_while",i,stack)
                return i+1
            self.set_parent(i,stack)
            i+=1

        self.handle_read_all_tokens_error( "process_do_while", stack )

    def process_try( self, i, stack ) :
        assert( self.all_tokens[i].spelling == "try" )
        if debug : self.print_entry("process_try",i,stack)

        self.set_parent(i,stack,"try")
        stack.append(self.all_tokens[i])
        i+=1
        found_catch = False
        while i < len(self.all_tokens ) :
            if not self.all_tokens[i].is_visible or self.all_tokens[i].is_inside_string:
                pass
            elif self.all_tokens[i].spelling == "{" :
                i = self.process_statement(i,stack)
                break
            self.set_parent(i,stack)
            i+=1
        j = self.find_next_visible_token(i,stack)
        while j < len( self.all_tokens ):
            #print("next visible: is it a catch?", j, self.all_tokens[j].spelling)
            if self.all_tokens[j].spelling == "catch" :
                for k in xrange(i,j) :
                    self.set_parent(k,stack)
                i = self.process_catch(j,stack)
                j = self.find_next_visible_token(i,stack)
                continue
            else :
                stack.pop() # remove try
                return i

        self.handle_read_all_tokens_error( "process_try", stack )
    def process_catch( self, i, stack ) :
        assert( self.all_tokens[i].spelling == "catch" )
        self.set_parent(i,stack,"catch")
        stack.append(self.all_tokens[i])
        i+=1
        while i < len( self.all_tokens ) :
            if not self.all_tokens[i].is_visible or self.all_tokens[i].is_inside_string:
                pass
            elif self.all_tokens[i].spelling == "(" :
                self.set_parent(i,stack,"catch-arg")
                stack.append(self.all_tokens[i])
                i+=1
                i = self.read_to_end_paren(i,stack)
                stack.pop() # remove catch-arg
                i = self.process_statement(i,stack)
                stack.pop() # remove catch
                return i
            self.set_parent(i,stack)
            i+=1
        self.handle_read_all_tokens_error( "process-catch", stack )


    def process_while(self,i,stack) :
        assert( self.all_tokens[i].spelling == "while" )
        if debug : self.print_entry("process_while",i,stack)
        self.set_parent(i,stack,"while")
        stack.append(self.all_tokens[i])
        i+=1
        while i < len(self.all_tokens) :
            if not self.all_tokens[i].is_visible or self.all_tokens[i].is_inside_string:
                pass
            elif self.all_tokens[i].spelling == "(" :
                self.set_parent(i,stack,"while-condition")
                stack.append(self.all_tokens[i])
                i+=1
                i = self.read_to_end_paren( i, stack )
                # ok -- now process the body of the while loop
                stack.pop() # remove while-condition
                i = self.process_statement( i, stack )
                stack.pop() # remove while
                if debug : self.print_entry("exitting process_while",i,stack)
                return i
            self.set_parent(i,stack)
            i+=1
        # we should not have arrived here
        self.handle_read_all_tokens_error( "process_while", stack )

    def process_namespace(self,i,stack) :
        if debug : self.print_entry("process_namespace",i,stack)
        self.set_parent(i,stack,"namespace")
        stack.append(self.all_tokens[i])
        i+=1
        while i < len(self.all_tokens) :
            i = self.find_next_visible_token( i, stack)
            if self.all_tokens[i].spelling == "{" :
                i = self.process_statement( i, stack )
                stack.pop()
                if debug :
                    if i < len(self.all_tokens) : self.print_entry("exitting process_namespace",i,stack)
                    else : print("exiting process_namespace; all tokens processed")
                return i
            if self.all_tokens[i].spelling == ";" :
                self.set_parent(i,stack)
                stack.pop()
                return i+1
            self.set_parent(i,stack)
            i+=1
        # we should not have arrived here
        self.handle_read_all_tokens_error( "process_namespace", stack )

    def process_class_decl(self,i,stack) :
        assert( self.all_tokens[i].spelling == "class" or self.all_tokens[i].spelling == "struct" )
        if debug : self.print_entry("process_class_decl",i,stack)

        seen_inheritance_colon = False
        found_square_brackets = False
        found_parens = False

        self.set_parent(i,stack,"class")
        stack.append(self.all_tokens[i])
        i+=1
        template_bracket_count = 0

        while i < len( self.all_tokens ) :
            if not self.all_tokens[i].is_visible or self.all_tokens[i].is_inside_string:
                pass
            elif self.all_tokens[i].spelling == "<" :
                template_bracket_count += 1
            elif self.all_tokens[i].spelling == ">" :
                template_bracket_count -= 1
            elif template_bracket_count == 0 and self.all_tokens[i].spelling == "[" :
                # ok -- maybe we're looking at an array declaration; this happens with structs
                # where you say struct structname varname[] = {};
                found_square_brackets = True
            elif self.all_tokens[i].spelling == "{" :
                if found_square_brackets :
                    assert( not seen_inheritance_colon );
                    i = self.process_to_end_rcb(i,stack,"array-value-initializer")
                    i = self.find_next_visible_token(i,stack)
                    assert( self.all_tokens[i].spelling == ";" )
                    self.set_parent(i,stack)
                    stack.pop() # remove class
                    return i+1
                else :
                    if seen_inheritance_colon :
                        stack.pop() # remove class-inheritance-list
                    i = self.process_statement( i, stack )
                    while i < len( self.all_tokens ) :
                        i = self.find_next_visible_token( i, stack )
                        self.set_parent( i, stack )
                        if self.all_tokens[i].spelling == ";" : break; # case where struct is anon and struct variables declared next
                        i+=1
                    stack.pop()
                    if debug : self.print_entry("exitting process_class_decl",i,stack)
                    return i+1
            elif self.all_tokens[i].spelling == ":" :
                if not seen_inheritance_colon :
                    if debug: print(" " * len(stack), "encountered inheritance colon")
                    seen_inheritance_colon = True
                    self.set_parent(i,stack,"class-inheritance-list")
                    stack.append(self.all_tokens[i])
                    i+=1
                    continue
            elif self.all_tokens[i].spelling == ";" :
                if seen_inheritance_colon :
                    # ok, we didn't actually see an inheritance colon; what we saw was a colon and probably because this
                    # isn't a class or struct declaration so much as it's a struct-instance declaration
                    # e.g. "struct mystruct varname = std::somefunction();"
                    stack.pop() # remove class-inheritance-list
                # this is just a forward declaration
                self.set_parent(i,stack)
                stack.pop()
                if debug : self.print_entry("exitting process_class_decl",i,stack)
                return i+1
            self.set_parent(i,stack)
            i += 1
        # we should not have arrived here
        self.handle_read_all_tokens_error( "process_class_decl", stack )

    def process_union(self,i,stack) :
        assert( self.all_tokens[i].spelling == "union" )
        if debug : self.print_entry("process_union",i,stack)
        self.set_parent(i,stack,"union")
        stack.append(self.all_tokens[i])
        i+=1
        processed_types = False
        while i < len( self.all_tokens ) :
            if not self.all_tokens[i].is_visible or self.all_tokens[i].is_inside_string:
                pass
            elif self.all_tokens[i].spelling == "{" and not processed_types :
                i = self.process_statement(i,stack)
                processed_types = True;
                continue
            elif self.all_tokens[i].spelling == ";" :
                assert( processed_types )
                self.set_parent(i,stack)
                stack.pop() # remove union
                if debug : self.print_entry("exitting process_union",i,stack)
                return i+1
            self.set_parent(i,stack)
            i+=1
        self.handle_read_all_tokens_error( "process_union", stack );

    def process_switch(self,i,stack) :
        assert( self.all_tokens[i].spelling == "switch" )
        if debug : self.print_entry("process_switch",i,stack)
        self.set_parent(i,stack,"switch")
        stack.append(self.all_tokens[i])
        i+=1
        while i < len(self.all_tokens) :
            if not self.all_tokens[i].is_visible or self.all_tokens[i].is_inside_string :
                pass
            if self.all_tokens[i].spelling == "(" :
                self.set_parent(i,stack,"switch-expression")
                stack.append(self.all_tokens[i])
                i+=1
                i = self.read_to_end_paren( i, stack )
                stack.pop() # remove switch-expression
                i = self.process_statement( i, stack )
                stack.pop() # remove switch
                if debug : self.print_entry("exitting process_switch",i,stack)
                return i
            self.set_parent(i,stack)
            i+=1
        self.handle_read_all_tokens_error( "process_class_decl", stack )

    def process_case(self,i,stack) :
        assert( self.all_tokens[i].spelling == "case" or self.all_tokens[i].spelling == "default" )
        if debug : self.print_entry("process_case",i,stack)

        self.set_parent(i,stack,"case")
        stack.append(self.all_tokens[i])
        i+=1
        while i < len(self.all_tokens) :
            if not self.all_tokens[i].is_visible or self.all_tokens[i].is_inside_string :
                pass
            elif self.all_tokens[i].spelling == ":"  :
                if self.all_tokens[i+1].spelling != ":" :
                  self.set_parent(i,stack,"case-block")
                  stack.append(self.all_tokens[i])
                  i+=1
                  # now process a bunch of statements
                  while i < len( self.all_tokens ) :
                      i = self.find_next_visible_token(i,stack)
                      i_spelling = self.all_tokens[i].spelling
                      if i_spelling == "case" or i_spelling == "default" or i_spelling == "}" :
                          stack.pop() # pop case-block
                          stack.pop() # pop case
                          if debug : self.print_entry("exitting process_case",i,stack)
                          return i
                      else :
                          i = self.process_statement(i,stack)
                else :
                    # jump past the next colon
                    self.set_parent(i,stack)
                    i+=1
            self.set_parent(i,stack)
            i+=1
        # we should not have gotten here
        self.handle_read_all_tokens_error( "process_case_or_default_case", stack )

    def process_to_end_rcb(self,i,stack,scopename) :
        if debug : self.print_entry( "process_to_end_rcb", i, stack )
        found_initial_lcb = False
        cbdepth = 0
        while i < len( self.all_tokens ) :
            if self.all_tokens[i].is_inside_string or not self.all_tokens[i].is_visible :
                pass
            elif self.all_tokens[i].spelling == "{" :
                cbdepth += 1
                if not found_initial_lcb :
                    self.set_parent(i,stack,scopename)
                    found_initial_lcb = True
                    stack.append(self.all_tokens[i])
                    i+=1
                    continue
            elif self.all_tokens[i].spelling == "}" :
                assert( found_initial_lcb )
                cbdepth -= 1
                if cbdepth == 0 :
                    self.set_parent(i,stack)
                    stack.pop()
                    return i+1
            self.set_parent(i,stack)
            i+=1
        #whoops!
        self.handle_read_all_tokens_error("process_read_to_end_rcb",stack)


    def process_enum(self,i,stack) :
        if debug : self.print_entry("process_enum",i,stack)
        assert( self.all_tokens[i].spelling == "enum" )
        self.set_parent(i,stack,"enum")
        stack.append(self.all_tokens[i])
        i+=1
        while i < len( self.all_tokens ) :
            i = self.find_next_visible_token(i,stack)
            assert( i < len(self.all_tokens))
            if self.all_tokens[i].spelling == "{" :
                i = self.process_to_end_rcb(i,stack,"enum-scope")
                continue
            elif self.all_tokens[i].spelling == ";" :
                self.set_parent(i,stack)
                stack.pop()
                return i+1
            self.set_parent(i,stack)
            i+=1
        # we should not have gotten here
        self.handle_read_all_tokens_error( "process_enum", stack )

    def process_scope(self,i,stack) :
        if debug : self.print_entry("process_scope",i,stack)

        scopename = stack[-1].type + "-scope"
        if scopename not in token_types :
            scopename = "scope"
        self.set_parent(i,stack,scopename)
        stack.append(self.all_tokens[i])
        i+=1
        while i < len(self.all_tokens) :
            i = self.find_next_visible_token(i,stack)

            assert( i < len(self.all_tokens) )
            if self.all_tokens[i].spelling == "}" :
                # print((" "*len(stack)), "process_scope", i, "found end }", self.all_tokens[i].line_number+1)
                #ok, done processing at this depth
                self.set_parent(i,stack)
                stack.pop() # remove "scope" or "*-scope"
                if debug : self.print_entry("exitting process_scope",i,stack)
                return i+1
            else :
                # print((" "*len(stack)), "process_scope encountered non-{, non-} at", i, self.all_tokens[i].spelling, self.all_tokens[i].line_number+1)
                i = self.process_statement(i,stack)
        self.handle_read_all_tokens_error( "process_scope", stack )

    # def process_scope_end(self,i,stack) :
    #     # print("process_scope_end", i, depth)
    #     self.set_parent(i,depth-1,stack.pop())
    #     return i+1, depth-1
    def process_privacy_declaration(self,i,stack) :
        if debug : self.print_entry("process_privacy_declaration",i,stack)
        self.set_parent(i,stack,"class-privacy")
        stack.append(self.all_tokens[i])
        i+=1
        while i < len(self.all_tokens) :
            if not self.all_tokens[i].is_visible or self.all_tokens[i].is_inside_string :
                pass
            elif self.all_tokens[i].spelling == ":" :
                self.set_parent(i,stack)
                stack.pop()
                if debug : self.print_entry("exitting process_privacy_declaration",i,stack)
                return i+1
            self.set_parent(i,stack)
            i+=1
        self.handle_read_all_tokens_error( "process_privacy_declaration", stack )


    def whole_line_invisible( self, line_number ) :
        for tok in self.line_tokens[ line_number ] :
            if tok.is_visible: return False
        return True

    def next_visible_token_on_line( self, line_number, starting_token_ind = 0 ) :
        i = starting_token_ind
        while i < len( self.line_tokens[line_number] )  :
            if self.line_tokens[line_number][i].is_visible : return i
        return i

    def determine_indentation_level( self, line_number ) :

        if len( self.line_tokens[ line_number ] ) == 0 : return

        # look at the first token on this line
        first_token = self.line_tokens[ line_number ][ 0 ]
        # print("first token: ", first_token.spelling, "" if not first_token.parent else first_token.parent.spelling)

        if first_token.spelling in self.macros_at_zero_indentation :
            self.line_indentations[ line_number ] = 0
        elif first_token.context() == "namespace" :
            # declaring new namespace -- use your parent's indentation
            if line_number != 0 :
                self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number]
        elif first_token.context() == "namespace-scope" :
            # namespaces scopes don't indent
            if line_number != 0 and first_token.parent.parent :
                self.line_indentations[line_number] = self.line_indentations[first_token.parent.parent.line_number]
        elif first_token.context() == "scope" :
            if first_token.spelling == "}" and first_token.is_visible and not first_token.is_inside_string :
                self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number]
            else :
                # indent one from the parent's indentation
                self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number] + 1
        elif first_token.context() == "for" :
            # this case might happen if a for-scope isn't inserted around a
            # for statement with only a single sub-statement and that statement
            # is on a second line; that should not happen, though.
            self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number] + 1
        elif first_token.context() == "for-declaration" :
            # this is where where the three statements in a for loop take more than one line
            # indent twice
            self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number] + 2
        elif first_token.context() == "for-scope" :
            if first_token.spelling == "}" and first_token.is_visible and not first_token.is_inside_string :
                self.line_indentations[line_number] = self.line_indentations[first_token.parent.parent.line_number]
            else :
                self.line_indentations[line_number] = self.line_indentations[first_token.parent.parent.line_number]+1
        elif first_token.context() == "do-while" :
            # this is an odd case -- between the "do" and the "{"
            self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number]
        elif first_token.context() == "do-while-condition" :
            # the do-while condition has occupied more than one line, I guess.  Indent twice like for and if
            self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number]+2
        elif first_token.context() == "do-while-scope" :
            if first_token.spelling == "}" and first_token.is_visible and not first_token.is_inside_string :
                self.line_indentations[line_number] = self.line_indentations[first_token.parent.parent.line_number]
            else :
                self.line_indentations[line_number] = self.line_indentations[first_token.parent.parent.line_number]+1
        elif first_token.context() == "if" :
            if first_token.type == "else" :
                self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number]
            elif first_token.spelling == "(" and first_token.is_visible :
                # the "(" at the beginning of the if-condition or a comment
                # after the if-condition but before the {
                self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number]+2
            else :
                # a comment after the if declaration either before an else declaration
                # or before the opening "{"
                self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number]
        elif first_token.context() == "if-condition" :
            # indent twice if the if-condition has gone on so long that it wraps to a second line
            self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number]+2
        elif first_token.context() == "if-scope" :
            if first_token.spelling == "}" and first_token.is_visible and not first_token.is_inside_string :
                self.line_indentations[line_number] = self.line_indentations[first_token.parent.parent.line_number]
            elif first_token.parent.parent.context() == "else" :
                self.line_indentations[line_number] = self.line_indentations[first_token.parent.parent.parent.line_number ]+1
            else :
                self.line_indentations[line_number] = self.line_indentations[first_token.parent.parent.line_number]+1
        elif first_token.context() == "else" :
            # indent to the same level as the if
            self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number]
        elif first_token.context() == "else-scope" :
            if first_token.spelling == "}" and first_token.is_visible and not first_token.is_inside_string :
                # print("else-scope indentation:", first_token.parent.parent.spelling, self.line_indentations[first_token.parent.parent.line_number])
                self.line_indentations[line_number] = self.line_indentations[first_token.parent.parent.line_number]
            else :
                self.line_indentations[line_number] = self.line_indentations[first_token.parent.parent.line_number]+1
        elif first_token.context() == "while" :
            # donno, this shouldn't really happen
            self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number]+2
        elif first_token.context() == "while-condition" :
            self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number]+2
        elif first_token.context() == "while-scope" :
            if first_token.spelling == "}" and first_token.is_visible and not first_token.is_inside_string :
                self.line_indentations[line_number] = self.line_indentations[first_token.parent.parent.line_number]
            else :
                self.line_indentations[line_number] = self.line_indentations[first_token.parent.parent.line_number]+1
        elif first_token.context() == "switch" :
            self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number]+2
        elif first_token.context() == "switch-expression" :
            self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number]+1
        elif first_token.context() == "switch-scope" :
            # don't indent for case: and default: statements, only for case-block and case-block-scopes.
            self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number]
            #if first_token.spelling == "}" and first_token.is_visible and not first_token.is_inside_string :
            #    self.line_indentations[line_number] = self.line_indentations[first_token.parent.parent.line_number]
            #else :
            #    self.line_indentations[line_number] = self.line_indentations[first_token.parent.parent.line_number]+1
        elif first_token.context() == "case" :
            # dunno -- is this a case condition that went for more than one line?
            self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number]
        elif first_token.context() == "case-block" :
            self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number]+1
            pass
        elif first_token.context() == "case-block-scope" :
            if first_token.spelling == "}" and first_token.is_visible and not first_token.is_inside_string :
                self.line_indentations[line_number] = self.line_indentations[first_token.parent.parent.line_number]
            else :
                self.line_indentations[line_number] = self.line_indentations[first_token.parent.parent.line_number]+1
        elif first_token.context() == "class" :
            # kind of an odd duck
            # indent once if we're in the process of declaring a class unless we're a "{"
            if self.line_tokens[line_number][0].spelling == "{" and first_token.is_visible and not first_token.is_inside_string :
                self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number]
            else :
                self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number]+1
        elif first_token.context() == "class-inheritance-list" :
            self.line_indentations[line_number] = self.line_indentations[first_token.parent.parent.line_number]+1
        elif first_token.context() == "class-scope" :
            # print("class scope ", line_number, " and parent on ", first_token.parent.line_number)
            if ( first_token.spelling == "}" or first_token.type == "class-privacy" ) and first_token.is_visible and not first_token.is_inside_string :
                self.line_indentations[line_number] = self.line_indentations[first_token.parent.parent.line_number]
            else :
                self.line_indentations[line_number] = self.line_indentations[first_token.parent.parent.line_number]+1
        elif first_token.context() == "class-privacy" :
            # don't indent privacy declarations
            self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number]
        elif first_token.context() == "union" :
            pass
        elif first_token.context() == "union-scope" :
            self.line_indentations[line_number] = self.line_indentations[first_token.parent.parent.line_number]+1
        elif first_token.context() == "function" :
            # we're in the middle of a function declaration -- e.g. "inline \n void \n etc"
            # or we're declaring a variable -- it's hard to tell
            # or its the first "(" opening the argument list, or the "{" opening the function-scope.
            self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number]
        elif first_token.context() == "function-decl-argument-list" :
            if first_token.spelling == ")" and first_token is first_token.parent.children[-1] :
                # the final closing parenthesis in an argument list gets indented at it's grandparent's indentation
                self.line_indentations[line_number] = self.line_indentations[first_token.parent.parent.line_number]
            else :
                # mid-argument-list
                self.line_indentations[line_number] = self.line_indentations[first_token.parent.parent.line_number]+1
        elif first_token.context() == "function-scope" :
            # use grand-parent's indentation as the reference; "function-scope" is the parent, "function" is the grandparent
            if first_token.spelling == "}" and first_token.is_visible and not first_token.is_inside_string :
                self.line_indentations[line_number] = self.line_indentations[first_token.parent.parent.line_number]
            else :
                self.line_indentations[line_number] = self.line_indentations[first_token.parent.parent.line_number]+1
        elif first_token.context() == "ctor-initializer-list" :
            # print("indentation for ctor-initializer-list", first_token.parent.parent.type)
            self.line_indentations[line_number] = self.line_indentations[first_token.parent.parent.line_number]+1
        elif first_token.context() == "statement" :
            #print("indentation for statement line?", first_token.spelling, first_token.index, "and", first_token.parent.spelling, first_token.parent.index)
            if first_token.parent.context() == "statement-scope" :
                #print("self.line_tokens[ first_token.parent.line_number ][0] is not first_token.parent", self.line_tokens[ first_token.parent.line_number ][0] is not first_token.parent)
                if ( self.line_tokens[ first_token.parent.line_number ][0] is not first_token.parent ) :
                    self.line_indentations[ line_number ] = self.line_indentations[first_token.parent.line_number] + 1
                else :
                    self.line_indentations[ line_number ] = self.line_indentations[first_token.parent.line_number]
            # this statement has run on to a second line
            elif first_token.is_visible and first_token.spelling == ")" and first_token.parent.children[-2] is first_token and \
                first_token.parent.children[-1].spelling == ";" :
                # the terminating ) in a function call, if it's the last token before the ";", should indent to
                # the same level as the parent.
                self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number]
            else :
                self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number]+1
        elif first_token.context() == "statement-scope" :
            # this statement has run on to a second line
            #print("indentation for statement scope line?", first_token.spelling, first_token.index, "and", first_token.parent.spelling, first_token.parent.index)
            if first_token.spelling == "}" and first_token.is_visible and not first_token.is_inside_string :
                self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number]
            else :
                self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number]+1
        elif first_token.context() == "template" :
            #print("indentation for template line?", first_token.spelling, first_token.index, "and", first_token.parent.spelling, first_token.parent.index)
            # this statement has run on to a second line
            self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number]
        elif first_token.context() == "template-arg-list" :
            #print("indentation for template-arg-list line?", first_token.spelling, first_token.index, "and", first_token.parent.spelling, first_token.parent.index)
            # this statement has run on to a second line
            if first_token.is_visible and first_token.spelling == ">" and first_token.parent.children[-1] is first_token :
                # the terminating > in a list of template arguments, should indent to
                # the same level as the parent.
                self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number]
            else :
                self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number]+1

        elif first_token.context() == "try-scope" :
            if first_token.spelling == "}" and first_token.is_visible and not first_token.is_inside_string :
                self.line_indentations[line_number] = self.line_indentations[first_token.parent.parent.line_number]
            else :
                self.line_indentations[line_number] = self.line_indentations[first_token.parent.parent.line_number]+1
        elif first_token.context() == "catch-scope" :
            if first_token.spelling == "}" and first_token.is_visible and not first_token.is_inside_string :
                self.line_indentations[line_number] = self.line_indentations[first_token.parent.parent.line_number]
            else :
                self.line_indentations[line_number] = self.line_indentations[first_token.parent.parent.line_number]+1
        elif first_token.context() == "enum" :
            if first_token.spelling == "{" and first_token.is_visible and not first_token.is_inside_string :
                self.line_indentations[ line_number ] = self.line_indentations[ first_token.parent.line_number]
            else :
                self.line_indentations[ line_number ] = self.line_indentations[ first_token.parent.line_number]+1
        elif first_token.context() == "enum-scope" :
            if first_token.spelling == "}" and first_token.is_visible and not first_token.is_inside_string :
                self.line_indentations[ line_number ] = self.line_indentations[ first_token.parent.parent.line_number]
            else :
                self.line_indentations[ line_number ] = self.line_indentations[ first_token.parent.parent.line_number]+1



    def line_from_line_tokens( self, line_number ) :
        toks = self.line_tokens[line_number]
        if len(toks) == 0 or not toks[0].invisible_by_macro :
            # print(line_number, [ x.spelling for x in toks ])
            spellings = ["\t" * self.line_indentations[line_number] ]
            for i, tok in enumerate( toks ) :
                if tok.spelling == "\t" :
                    assert( tok.is_inside_string )
                    spellings.append( "\\t" ) # output tabs inside strings as escaped tab sequences
                else :
                    spellings.append( tok.spelling )
                if i+1 != len( toks ) :
                    spellings.append( " " * ( toks[i+1].start - tok.one_past_end ) ) #preserve spaces
            spellings.append("\n")
            return "".join( spellings )
        else :
            # do not modify lines that are invisible
            return self.all_lines[ toks[0].orig_line_number ]

    def insert_space_after( self, tok ) :
        line_toks = self.line_tokens[ tok.line_number ]
        i = 0
        while i < len(line_toks) :
            if line_toks[i] is tok :
                i+=1
                break
            i+=1
        if i == len(line_toks) :
            # whoa!
            print("Failed to find tok:", tok.spelling, "on line", tok.line_number, "in file", self.filename)
            sys.exit(1)
        while i < len(line_toks) :
            line_toks[i].start = line_toks[i].start+1
            line_toks[i].one_past_end = line_toks[i].one_past_end+1
            i+=1

    def enforce_space_between( self, tok1, tok2 ) :
        # print("tok1", tok1.spelling, tok1.line_number, tok1.start, tok1.one_past_end)
        # print("tok2", tok2.spelling, tok2.line_number, tok2.start, tok2.one_past_end)
        if tok1.line_number == tok2.line_number and tok1.one_past_end == tok2.start :
            self.insert_space_after( tok1 )

    def correct_spacing( self ) :
        for i,tok in enumerate(self.all_tokens) :
            if not tok.is_visible or tok.is_inside_string : continue

            if i > 0 :
                prevtok = self.all_tokens[i-1]
            else :
                prevtok = None
            if i+1 != len(self.all_tokens) :
                nexttok = self.all_tokens[i+1]
            else :
                nexttok = None

            if tok.context() == "if" :
                if tok.spelling == "(" :
                    # ok, let's check to make sure the "if" statement before this is properly spaced
                    assert( tok.parent )
                    self.enforce_space_between( tok.parent, tok )
                    assert( nexttok )
                    self.enforce_space_between( tok, nexttok )
                elif tok.spelling == "{" :
                    # ok, let's make sure there's a space after the ")" that proceeded this
                    assert( prevtok )
                    self.enforce_space_between( prevtok, tok )
                    assert( nexttok )
                    if nexttok.spelling != "}" :
                        self.enforce_space_between( tok, nexttok )
                elif tok.spelling == "else" :
                    # make sure there's a space between the if's "}" and the next "{"
                    assert( prevtok )
                    if prevtok.spelling == "}":
                        self.enforce_space_between( prevtok, tok )
                    if nexttok.spelling == "{" :
                        self.enforce_space_between( tok, nexttok )
            elif tok.context() == "if-condition" :
                if tok.spelling == ")" and tok.parent.children[-1] is tok :
                    self.enforce_space_between( prevtok, tok )
                    self.enforce_space_between( tok, nexttok )
            elif tok.context() == "else" :
                if tok.spelling == "{" and nexttok.spelling != "}" :
                    self.enforce_space_between( tok, nexttok )
            elif tok.context() == "else-scope" :
                if tok.spelling == "}" and prevtok.spelling != "}" :
                    self.enforce_space_between( prevtok, tok )
            elif tok.context() == "for" :
                if tok.spelling == "(" :
                    self.enforce_space_between( prevtok, tok )
                    self.enforce_space_between( tok, nexttok )
                elif tok.spelling == "{" :
                    self.enforce_space_between( prevtok, tok ) # for (..)_{
            elif tok.context() == "for-declaration" :
                if tok.spelling == ";" :
                    self.enforce_space_between( tok, nexttok ) # for ( int i = 1;_i<=..), for (..;_++i )
                elif tok.spelling == ")" and tok.parent.children[-1] is tok :
                    self.enforce_space_between( prevtok, tok ) # for (..++i_)
                    if nexttok.type != "empty-statement" :
                        self.enforce_space_between( tok, nexttok )
            elif tok.context() == "while" :
                if tok.spelling == "(" :
                    self.enforce_space_between( prevtok, tok ) # while_(
                    self.enforce_space_between( tok, nexttok ) # while (_condition
            elif tok.context() == "while-condition" :
                if tok.spelling == ")" and tok.parent.children[-1] is tok :
                    self.enforce_space_between( prevtok, tok ) # while ( condition_)
                    self.enforce_space_between( tok, nexttok ) # while ( .. )_{
            elif tok.context() == "case" :
                if tok.spelling == ":" and nexttok.context() == "case-block" :
                    self.enforce_space_between( prevtok, tok )
                    self.enforce_space_between( tok, nexttok )
            # elif tok.context() == "function" :
            #     if tok.spelling == "(" and nexttok.spelling != ")" :
            #         self.enforce_space_between( tok, nexttok )
            # elif tok.context() == "function-decl-argument-list" :
            #     if tok.spelling == ")" and tok.parent.children[-1] is tok and prevtok.spelling != "(" :
            #         self.enforce_space_between( prevtok, tok )
            elif tok.context() == "namespace" :
                if tok.spelling == "{" :
                    self.enforce_space_between( prevtok, tok )

    def last_descendent( self, token ) :
        # return the last descendent of a given token
        # does not require that the index of the token be up-to-date
        # print(token.spelling, token.line_number)
        if not token.children :
            return token
        else :
            return self.last_descendent( token.children[-1] )

    #def excise_left_curly_brace_and_move_after( self, lcb_tok, dest_tok ) :
    def excise_token_and_move_upwards_and_after( self, moving_token, dest_tok, adjust_parentage=True, space_from_prev=1 ) :
        # move moving_token up to the line that dest_tok is on.
        # dest_tok == the token after which the moving token is to be placed.
        # if adjust_parentage is true, then if there were tokens between moving_token and
        # dest_tok, adjust their parent to be the moving_token -- this must make sense.
        # i.e. a for statement may only have one child besides the for-declaration child.

        orig_line_number = moving_token.line_number
        orig_line = self.line_tokens[ orig_line_number ]
        dest_line_number = dest_tok.line_number
        dest_line = self.line_tokens[ dest_line_number ]
        dest_tok_ind = dest_line.index( dest_tok )
        moving_token.line_number = dest_line_number
        moving_token_length = moving_token.one_past_end - moving_token.start
        moving_token.start = dest_tok.one_past_end + space_from_prev
        moving_token.one_past_end = moving_token.start + moving_token_length
        # shift the rest of the tokens on this line right by two spaces
        for i in xrange( dest_tok_ind+1, len(dest_line) ) :
            tok = dest_line[i]
            tok.start += space_from_prev+moving_token_length
            tok.one_past_end += space_from_prev+moving_token_length
        dest_line.insert( dest_tok_ind+1, moving_token )
        if adjust_parentage :
            new_children = []
            for i in xrange( dest_line_number, orig_line_number+1 ) :
                found_moving_token = False
                for tok in self.line_tokens[ i ] :
                    if i == dest_line_number and not found_moving_token :
                        if tok is moving_token :
                            found_moving_token = True
                        continue
                    if i == orig_line_number :
                        if tok is moving_token :
                            break
                    orig_parent = tok.parent
                    orig_parent.children.remove( tok )
                    # print("move upward:", tok.spelling, "added as child of", moving_token.spelling)
                    tok.parent = moving_token
                    new_children.append( tok )

            if new_children :
                moving_token.children = new_children + moving_token.children
                # we also have to adjust self.all_tokens
                # i = self.all_tokens.index( moving_token )
                # while i > 0 :
                #     self.all_tokens[i] = self.all_tokens[i-1]
                #     if self.all_tokens[i] is new_children[0] :
                #         self.all_tokens[i-1] = moving_token
                #         break
                #     i-=1
                # assert( i != 0 )
        orig_line.remove( moving_token )

    def delete_empty_line( self, old_line ) :
        # iterate across all the tokens and decrement the line_number for any that are after
        # this old line number; also delete one line from the line_toks table
        for tok in self.all_tokens :
            if tok.line_number > old_line : tok.line_number -= 1
        self.line_tokens = self.line_tokens[:old_line] + self.line_tokens[old_line+1:]


    def insert_token_lines( self, tok_lines ) :
        # tokens must already know their new line numbers.
        # first and last lines of tok_lines can't be empty.
        # adjusts self.line_tokens, but not self.all_tokens.
        # print("insert token lines: ", [ " ".join( [ x.spelling for x in line ]) for line in tok_lines ])
        nlines = len( tok_lines )
        assert( len(tok_lines[-1]) > 0 and len(tok_lines[0]) > 0 )
        firstline = tok_lines[0][0].line_number
        self.line_tokens = self.line_tokens[:firstline] + tok_lines + self.line_tokens[firstline:]
        for i in xrange( tok_lines[-1][0].line_number + 1, len(self.line_tokens) ) :
            for tok in self.line_tokens[i] :
                tok.line_number += nlines

    def insert_line_with_single_token( self, tok ) :
        self.insert_token_lines( [[ tok ]] )

    def excise_right_curly_brace_and_move_to_next_line_on_own( self, rcb ) :
        rcb_line = self.line_tokens[ rcb.line_number ]

        if rcb is not rcb_line[ -1 ] :
            # ok, we have comments or something at the end of this line
            # if it's nothing but comments, then leave them on this line
            # otherwise, they should be on their own line after the }
            rcb_ind = rcb_line.index( rcb )

            any_visible = False
            for i in xrange( rcb_ind+1, len( rcb_line )) :
                if self.line_tokens[ rcb.line_number ][ i ].is_visible :
                    any_visible = True
                    break;
            if any_visible :
                # ok! then we need to move these statements to their own line
                # we don't have to ajust parentage or the all_tokens list
                toks_to_move = rcb_line[ rcb_ind: ]
                self.line_tokens[ rcb.line_number ] = rcb_line[:rcb_ind]
                rcb.line_number += 1
                newlines = []
                newlines.append( [ rcb ] )
                newlines.append( toks_to_move[1:] )
                for tok in newlines[1] :
                    tok.line_number += 2
                self.insert_token_lines( newlines )
                return;
            else :
                # case where there are comments at the end of the line after rcb.
                # let's not adjust the start and one_past_end positions of these tokens
                # but do adjust the tree so that these tokens are assigned rcb's parent
                for i in xrange(rcb_ind+1,len(rcb_line)) :
                    tok = rcb_line[i]
                    orig_parent = tok.parent
                    orig_parent.children.remove( tok )
                    tok.set_parent( rcb.parent )
                # update the all_tokens list
                # Dont bother i = self.all_tokens.index( rcb )
                # Dont bother while i < len(self.all_tokens)-1 :
                # Dont bother     self.all_tokens[i] = self.all_tokens[i+1]
                # Dont bother     if self.all_tokens[i] is rcb_line[-1] :
                # Dont bother         self.all_tokens[i+1] = rcb
                # Dont bother         break
                # Dont bother     i+=1
                # Dont bother assert( i < len(self.all_tokens)-1 )
        else :
            # the rcb is the last token on its line
            # no need to adjust parentage or the all_tokens list
            pass
        rcb_line.remove( rcb )
        rcb.line_number += 1
        self.insert_line_with_single_token( rcb )


    def add_left_curly_brace_after( self, tok_before ) :
        lcb = Token()
        lcb.spelling = "{"
        lcb.start = tok_before.one_past_end + 1
        lcb.one_past_end = lcb.start + 1
        lcb.line_number = tok_before.line_number
        line_toks = self.line_tokens[ tok_before.line_number ]
        tok_before_index = line_toks.index( tok_before )
        for i in xrange( tok_before_index+1, len( line_toks ) ) :
            line_toks[i].start += 2
            line_toks[i].one_past_end += 2
        line_toks.insert( tok_before_index+1, lcb )
        #self.all_tokens.insert( self.all_tokens.index( tok_before ) + 1, lcb )
        self.all_tokens.append( lcb )
        return lcb

    def add_right_curly_brace_on_own_line_after( self, tok_before, lcb_tok ) :
        rcb_tok = Token()
        rcb_tok.spelling = "}"
        rcb_tok.start = 0 # don't worry about the indentation
        rcb_tok.one_past_end = 1
        rcb_tok.line_number = tok_before.line_number + 1
        #self.all_tokens.insert( self.all_tokens.index( tok_before )+1, rcb_tok )
        self.all_tokens.append( rcb_tok )
        self.insert_line_with_single_token( rcb_tok )
        rcb_tok.set_parent( lcb_tok )

    def move_tokens_after_to_their_own_next_line( self, tok ) :
        tok_line = self.line_tokens[ tok.line_number ]
        tok_pos = tok_line.index( tok )
        toks_for_next_line = tok_line[ (tok_pos+1): ]
        # print("toks_for_next_line", [ x.spelling for x in toks_for_next_line ])
        for nltok in toks_for_next_line :
            nltok.line_number += 1
        self.line_tokens[ tok.line_number ] = tok_line[:(tok_pos+1)]
        self.insert_token_lines( [ toks_for_next_line ] )

    def adjust_token_tree( self, new_parent, child_begin, child_end ) :
        # takes an inclusive range of children between child_begin and child_end
        # of a particular node (old_parent) and makes all of them children
        # of new_parent; then new_parent is inserted as a child of new_parent
        # into the position that child_begin formerly occupied.
        assert( child_begin.parent is child_end.parent )
        old_parent = child_begin.parent
        children_copy = list( old_parent.children )
        child_begin_index = old_parent.children.index( child_begin )
        i = child_begin_index
        while i < len( children_copy ) :
            child = children_copy[i]
            old_parent.children.remove( child )
            child.set_parent( new_parent )
            if child is child_end :
                break
            i+=1
        assert( i < len( children_copy ))
        new_parent.parent = old_parent
        old_parent.children.insert( child_begin_index, new_parent )


    def put_empty_curly_braces_on_line_together( self, lcb, rcb ) :
        # move the rcb up to the lcb
        self.line_tokens[ rcb.line_number ].remove( rcb )
        old_line_number = rcb.line_number
        self.line_tokens[ lcb.line_number ].append( rcb )
        rcb.line_number = lcb.line_number
        rcb.start = lcb.one_past_end
        rcb.one_past_end = rcb.start+1
        # print("toks left on rcb line:", ", ".join( [x.spelling for x in self.line_tokens[old_line_number] ]))
        self.delete_empty_lines_between( lcb.line_number+1, old_line_number )

    def visible_tokens_on_line_after( self, tok ) :
        line = self.line_tokens[ tok.line_number ]
        ind = line.index( tok )
        i = ind+1
        while i < len( line ) :
            if line[i].is_visible :
                return True
            i+=1
        return False

    def next_visible_token_on_line_after( self, tok ) :
        line = self.line_tokens[ tok.line_number ]
        i = line.index( tok )+1
        while i < len( line ) :
            if line[i].is_visible :
                return line[i]
            i+=1
        return -1

    def token_is_scope_ender( self, tok ) :
        return tok.is_visible and not tok.is_inside_string and tok.spelling == "}" \
            and tok.parent.children[-1] is tok and tok.parent.type in self.scope_types

    def recurse_to_ifelse_block_parent_if( self, elsetok ) :
        ifparent = elsetok.parent
        if ifparent.context() == "else" :
            return self.recurse_to_ifelse_block_parent_if( ifparent.parent )
        else :
            return ifparent
    def find_last_visible_token_on_line( self, line_number ) :
        line = self.line_tokens[line_number]
        i = len(line)-1
        while i > 0 :
            if line[i].is_visible : return i
            i -=1
        return i

    def adjust_parentage_of_comments_after( self, line, first_to_adjust, comments_new_parent, prev_child_of_comment ) :
        # print("prev_child_of_comment:", "None" if not prev_child_of_comment else prev_child_of_comment.spelling)
        new_child_ind = comments_new_parent.children.index( prev_child_of_comment )+1 if prev_child_of_comment else 0
        # print("new_child_ind:", new_child_ind)
        for i in xrange( first_to_adjust, len(line) ) :
            itok = line[i]
            itok.parent.children.remove( itok )
            itok.parent = comments_new_parent
            comments_new_parent.children.insert( new_child_ind, itok )
            # print("inserting child", itok.spelling, "at index", new_child_ind, ": ", ", ".join( [x.spelling for x in comments_new_parent.children ]))
            new_child_ind += 1
        return new_child_ind

    def get_parent_for_eol_comments_for_else( self, else_tok, last_desc, line_number ) :
        # if we're moving the entirety of line_number upwards in front of existing comments,
        # then those comments need to be assigned a new parent to rework the tree.
        # Return last_visible_ind, the index of the last token on this line that is visible
        # (this line should have at least one visible token on it), the new parent for
        # the comments at the end of the line, and the child of this new parent that should be
        # the older child of the those comments

        last_visible_ind = self.find_last_visible_token_on_line( line_number )
        assert( last_visible_ind != -1 )
        last_visible = self.line_tokens[ line_number ][ last_visible_ind ]

        if line_number == last_desc.line_number :
            iftok = self.recurse_to_ifelse_block_parent_if(else_tok)
            comments_new_parent = iftok.parent
            prev_child = iftok # comments will be appended after iftok in the same context as iftok
        else :
            if len( last_visible.children ) != 0 :
                # if there are children of this node, then it can take children -- make the comment a child of this node
                comments_new_parent = last_visible
                prev_child = None
            else :
                # if there are no children of this node, then this node's parent can take children
                last_desc_of_last_vis_parent = self.last_descendent( last_visible.parent )
                if last_visible is not last_desc_of_last_vis_parent :
                    # we're in the middle of the list of descendents from the parent of this node.
                    # make the eol-comments descend form this node's parent.
                    comments_new_parent = last_visible.parent
                    prev_child = last_visible
                else :
                    # we've found the last child of this parent -- so look upwards to the grandparent
                    # or beyond, to figure out where we should insert these comments
                    parent = last_visible.parent
                    grandparent = last_visible.parent.parent
                    while grandparent.parent :
                        # if the grandparent is the top-level node, go ahead and stop recursing
                        if parent is not grandparent.children[-1] :
                            # we've found a node where the parent is not the last child of it's parent;
                            # we may stop
                            break
                        else :
                            parent = grandparent
                            grandparent = parent.parent

                    comments_new_parent = grandparent
                    prev_child = parent

        return last_visible_ind, comments_new_parent, prev_child

    def move_line_upwards_and_after( self, tok, tok_before, last_visible_ind, comments_new_parent, prev_child ) :
        # take everything after tok on its line and move it upwards to just after the tok_before token

        orig_line_number = tok.line_number
        orig_line = list( self.line_tokens[ tok.line_number ] )
        tok_line_index = orig_line.index( tok )
        i = tok_line_index
        before_tok = tok_before
        last_tok = self.line_tokens[ tok_before.line_number ][-1]
        if last_tok.is_visible :
            # if there were no end-of-line comments at the end of the destination line,
            # don't really change the way comments appear at the end of this line.
            last_tok = orig_line[ last_visible_ind ]
        prev_end = 0

        while i < len( orig_line ) :
            itok = orig_line[i]
            itok_end = itok.one_past_end
            if itok.is_visible or i <= last_visible_ind :
                self.excise_token_and_move_upwards_and_after( itok, before_tok, False, 1 if i==tok_line_index else itok.start-prev_end )
                before_tok = itok
            else :
                # move already-at-the-end-of-the-line comments to the end of the line
                self.excise_token_and_move_upwards_and_after( itok, last_tok, False, 1 if i+1==last_visible_ind else itok.start-prev_end )
                # itok.parent.children.remove(itok)
                # itok.set_parent(comments_new_parent)
                last_tok = itok
            prev_end = itok_end
            i+=1
        self.adjust_parentage_of_comments_after( orig_line, last_visible_ind+1, comments_new_parent, prev_child )

        # now adjust the parentage of all tokens on lines between the moving line and the destination line
        # the only thing between them
        #print("now adjust the parentage of all tokens on lines between the moving line and the destination line")
        #print("comments_new_parent:", comments_new_parent.spelling)
        #print("prev_child:", "None" if not prev_child else prev_child.spelling)
        for line_number in xrange( tok_before.line_number+1, orig_line_number+1 ) :
            for tok in self.line_tokens[line_number] :
                # print("tok on line", line_number, tok.spelling, tok.is_visible, tok.invisible_by_macro, tok.is_preprocessor_directive)
                assert( not tok.is_visible and not tok.invisible_by_macro and not tok.is_preprocessor_directive )
            new_child_ind_plus_one = self.adjust_parentage_of_comments_after( self.line_tokens[line_number], 0 , comments_new_parent, prev_child )
            #print("new_child_ind_plus_one:", new_child_ind_plus_one)
            if new_child_ind_plus_one != 0 :
                prev_child = comments_new_parent.children[ new_child_ind_plus_one - 1 ]

        self.delete_empty_lines_between( tok_before.line_number, orig_line_number )


    def find_rcb_child( self, tok ) :
        # search through the children on tok and look for the right-curly-brace child
        # (there should be only one)
        rcbs = filter( lambda x : x.is_visible and not x.is_inside_string and x.spelling == "}", tok.children )
        assert( len(rcbs) == 1 )
        return rcbs[0]

    def intervening_macro_tokens( self, tok1, tok2 ) :
        # return True if there are any tokens between tok1 and tok2
        # in the self.line_tokens lists which are either preprocessor directives
        # or that are invisible by macros.
        # tok1 must preceed tok2.
        line1 = tok1.line_number
        line2 = tok2.line_number
        for line_number in xrange( line1, line2+1 ) :
            toks = self.line_tokens[line_number]
            ind = 0
            if line_number == line1 :
                ind = toks.index( tok1 )+1
            while ind < len(toks) :
                tok = toks[ind]
                if tok is tok2 : break
                if tok.is_preprocessor_directive or tok.invisible_by_macro : return True
                ind+=1
        return False

    def delete_empty_lines_between( self, start, end ) :
        for line_number in xrange( end, start-1, -1 ) :
            if len(self.line_tokens[line_number]) == 0 :
                self.delete_empty_line( line_number )


    def adjust_else( self, tok ) :
        # print("adjust else", tok.line_number)
        last_desc = self.last_descendent( tok )
        else_body = None
        for child in tok.children :
            if child.is_visible :
                else_body = child
                break
        # ok -- first check if the previous child of the parent was an if-scope
        # and if so, and if the if-scope starts and ends on different lines, then
        # make sure the else is moved up to the line the if-scope ends on
        parent_if = tok.parent
        # for child in parent_if.children :
        #     print("parent if child:", child.spelling, child.type)
        else_ind = parent_if.children.index( tok )
        # print("else ind:", else_ind)
        prev_visible_child_ind = else_ind - 1
        assert( else_ind != 0 )
        while prev_visible_child_ind >= 0 :
            # print("prev_child:", parent_if.children[ prev_visible_child_ind ].spelling)
            if parent_if.children[ prev_visible_child_ind ].is_visible :
                break
            prev_visible_child_ind -= 1
        assert( prev_visible_child_ind != -1 )
        prev_child = parent_if.children[ prev_visible_child_ind ]
        if prev_child.type == "if-scope" :
            # if debug : print("else children:", ", ".join( [ x.spelling for x in tok.children ] ))
            # if debug : print("else.child[0].children:", ", ".join( [ x.spelling for x in tok.children[0].children ] ))
            # if debug : print("prev_child children:", ", ".join( [ x.spelling for x in prev_child.children ] ))
            if_lcb = prev_child #rename for clarity
            if_rcb = self.find_rcb_child( if_lcb )
            assert( if_rcb.spelling == "}" )
            if if_rcb.line_number != tok.line_number and if_lcb.line_number != if_rcb.line_number and not self.intervening_macro_tokens( if_rcb, tok ) :
                # ok -- move the "else" up to behind the rcb
                last_visible_ind, comments_new_parent, prev_child = self.get_parent_for_eol_comments_for_else( tok, last_desc, tok.line_number )
                # print("last_visible_ind", last_visible_ind, "comments_new_parent", comments_new_parent.spelling, comments_new_parent.type, "prev_child", prev_child.spelling if prev_child else "None")
                self.move_line_upwards_and_after( tok, if_rcb, last_visible_ind, comments_new_parent, prev_child )
                # print("comments new parent children:'%s'" % "', '".join( [ x.spelling for x in comments_new_parent.children ] ))

            #otherwise, leave the else in place.

        if else_body.type == "if" :
            if tok.line_number != else_body.line_number and not self.intervening_macro_tokens( tok, else_body ):
                # ok, make sure that the whole line that the "if" statement is on ends up on this line
                last_visible_ind, comments_new_parent, prev_child = self.get_parent_for_eol_comments_for_else( tok, last_desc, else_body.line_number )
                # print("moving if behind else upwards", last_visible_ind, comments_new_parent.spelling, "None" if not prev_child else prev_child.spelling)
                self.move_line_upwards_and_after( else_body, tok, last_visible_ind, comments_new_parent, prev_child )
        elif else_body.type == "else-scope" :
            else_lcb = else_body
            else_rcb = else_body.children[-1]
            assert( else_lcb.spelling == "{" and else_rcb.spelling == "}" )
            if tok.line_number != else_lcb.line_number and not self.intervening_macro_tokens(tok, else_body) :
                # move the else_lcb up to the line with the else
                orig_line_number = else_lcb.line_number
                self.excise_token_and_move_upwards_and_after( else_lcb, tok )
                if len( self.line_tokens[ orig_line_number ] ) == 0 : self.delete_empty_line( orig_line_number )

            if len( else_lcb.children ) != 1 :
                # make sure the rcb is not on the same line as the last token
                # of the second-to-last child of the else_lcb unless it's also on
                # the same line as the lcb
                last_before_rcb = self.last_descendent( else_lcb.children[-2] )
                if last_before_rcb.line_number == else_rcb.line_number and else_rcb.line_number != else_body.line_number :
                    # print('ok -- move rcb to its own line')
                    self.excise_right_curly_brace_and_move_to_next_line_on_own( else_rcb )
            else :
                # ok, then else_lcb has a single child, the else_rcb.
                # make sure the {}s are on the same line together
                if else_lcb.line_number != else_rcb.line_number :
                    self.put_empty_curly_braces_on_line_together( else_lcb, else_rcb )
        else :
            # check to see if we ought to have curly braces
            last_es_desc = self.last_descendent( else_body )
            if last_es_desc.line_number != tok.line_number and not self.intervening_macro_tokens(tok,else_body) :
                # print('we need to add a pair of curly braces')
                lcb_tok = self.add_left_curly_brace_after( tok )
                lcb_tok.type = "else-scope"
                self.adjust_token_tree( lcb_tok, tok.children[0], else_body ) # start with the first child of the else (possibly a comment)
                if else_body.line_number == lcb_tok.line_number :
                    # print('move these tokens down to their own line')
                    self.move_tokens_after_to_their_own_next_line( lcb_tok )
                if self.visible_tokens_on_line_after( last_es_desc ) :
                    print("Possible bug identified on line", last_es_desc.orig_line_number+1, "of file", self.filename, "1")
                    self.move_tokens_after_to_their_own_next_line( last_es_desc )
                self.add_right_curly_brace_on_own_line_after( last_es_desc, lcb_tok )
            else :
                # make sure that there are no other tokens on this line except a scope-ending "}"
                #if else_body.type != "if" and
                if self.line_tokens[ last_es_desc.line_number ][ -1 ] is not last_es_desc :
                    if self.visible_tokens_on_line_after( last_es_desc ) :
                        next_visible = self.next_visible_token_on_line_after( last_es_desc )
                        if not self.token_is_scope_ender( next_visible ) :
                            # print('move these tokens onto their own line')
                            print("Possible bug identified on line", last_es_desc.orig_line_number+1, "of file", self.filename, "2")
                            self.move_tokens_after_to_their_own_next_line( last_es_desc )
        self.adjust_lines_for_children( tok ) # recurse

    def adjust_if( self, tok ) :
        #print("adjust if", tok.line_number)
        # ok -- let's look at the descendents; they should be an if-condition and
        # a statement
        last_desc = self.last_descendent( tok )
        if_cond = None
        if_body = None
        for child in tok.children :
            if child.type == "if-condition" :
                if_cond = child
            elif child.is_visible :
                if_body = child
                break # the last descendent may be an else statement
        end_paren = if_cond.children[-1]
        assert( end_paren.spelling == ")" )
        if if_body.type == "if-scope" :
            # ok -- so the curly braces are there.  Check that they are in the right place.
            if end_paren.line_number != if_body.line_number and not self.intervening_macro_tokens( end_paren, if_body ) :
                # print('ok! we need to move the "{" up to the ")" line')
                if len( if_body.children ) == 1 :
                    # then we're looking at an empty statement
                    # first put the } on the same line as the {
                    # and then move the {} like upwards and after the )
                    lcb = if_body
                    rcb = if_body.children[0]
                    if lcb.line_number != rcb.line_number :
                        orig_rcb_line_number = rcb.line_number
                        self.excise_token_and_move_upwards_and_after( rcb, lcb )
                        # print("toks left on rcb line:", ", ".join( [x.spelling for x in self.line_tokens[orig_rcb_line_number] ]))
                        rcb.start = lcb.one_past_end; rcb.one_past_end = rcb.start+1
                        self.delete_empty_lines_between( lcb.line_number+1, orig_rcb_line_number)
                    # now move the lcb and rcb up behind the end_paren
                    last_visible_ind = self.find_last_visible_token_on_line( rcb.line_number )
                    comments_new_parent = tok.parent
                    prev_child = tok
                    orig_lcb_line_number = lcb.line_number
                    self.move_line_upwards_and_after( lcb, end_paren, last_visible_ind, comments_new_parent, prev_child )
                    self.delete_empty_lines_between( end_paren.line_number+1, orig_lcb_line_number)
                    # print("line start stop", ", ".join( [ "%s %d %d" % ( x.spelling, x.start, x.one_past_end ) for x in self.line_tokens[end_paren.line_number] ]))
                else :
                    old_line = if_body.line_number
                    self.excise_token_and_move_upwards_and_after( if_body, end_paren )
                    if len( self.line_tokens[ old_line] ) == 0 :
                        self.delete_empty_line( old_line )
            # print("if_body.children", ", ".join( [x.spelling for x in if_body.children ]))
            right_curly_brace = if_body.children[-1]
            # print("right_curly_brace", right_curly_brace.spelling, right_curly_brace.line_number)
            assert( right_curly_brace.spelling == "}" )
            # the right curly brace should be on its own line unless its on the same line as the left curly brace
            if if_body.line_number != right_curly_brace.line_number and len( if_body.children ) != 1 :
                last_before_rcb = self.last_descendent( if_body.children[-2] )
                if last_before_rcb.line_number == right_curly_brace.line_number :
                    # print('ok! move rcb to its own line')
                    self.excise_right_curly_brace_and_move_to_next_line_on_own( right_curly_brace )
            elif len( if_body.children ) == 1 :
                # make sure the {}s are on the same line together
                if right_curly_brace.line_number != if_body.line_number and not self.intervening_macro_tokens( if_body, right_curly_brace ) :
                    self.put_empty_curly_braces_on_line_together( if_body, right_curly_brace )
        else :
            # ok -- just a single unscoped statement
            last_is_desc = self.last_descendent( if_body )
            last_ic_desc = self.last_descendent( if_cond )
            if last_is_desc.line_number != last_ic_desc.line_number and not self.intervening_macro_tokens(tok,if_body) :
                #print('we need to add a pair of curly braces')
                lcb_tok = self.add_left_curly_brace_after( end_paren )
                lcb_tok.type = "if-scope"
                self.adjust_token_tree( lcb_tok, tok.children[ tok.children.index(if_cond)+1 ], if_body )
                if if_body.line_number == lcb_tok.line_number :
                    #print('move these tokens down to their own line')
                    self.move_tokens_after_to_their_own_next_line( lcb_tok )
                if self.visible_tokens_on_line_after( last_is_desc ) :
                    print("Possible bug identified on line", last_is_desc.orig_line_number+1, "of file", self.filename, "3")
                    self.move_tokens_after_to_their_own_next_line( last_is_desc )
                self.add_right_curly_brace_on_own_line_after( last_is_desc, lcb_tok )
            else :
                # make sure that there are no other tokens on this line, except an "else" or a scope-ending "}"
                if self.line_tokens[ last_is_desc.line_number ][ -1 ] is not last_is_desc :
                    if self.visible_tokens_on_line_after( last_is_desc ) :
                        # print('move these tokens onto their own line')
                        next_visible = self.next_visible_token_on_line_after( last_is_desc )
                        if not self.token_is_scope_ender( next_visible ) and next_visible.type != "else" :
                            print("Possible bug identified on line", last_is_desc.orig_line_number+1, "of file", self.filename, "4")
                            self.move_tokens_after_to_their_own_next_line( last_is_desc )
        self.adjust_lines_for_children( tok ) # recurse


    def adjust_for( self, tok ) :
        # ok -- let's look at the descendents; they should be
        # the for-declaration descendent and then either a statement of some form
        # or a for-scope

        last_desc = self.last_descendent( tok )
        for_dec_tok = None
        for_body = None

        for child in tok.children :
            if child.type == "for-declaration" :
                for_dec_tok = child
            elif child.is_visible :
                for_body = child

        end_paren = for_dec_tok.children[-1]
        assert( end_paren.spelling == ")" )

        if for_body.type == "for-scope" :
            # ok -- so the curly braces are there.  Check that they are in the right place.
            if end_paren.line_number != for_body.line_number :
                # print('ok! we need to move the "{" up to the ")" line')
                # ask -- are there other tokens on the "{" line?
                if for_body.children[0].spelling == "}" :
                    pass; # handled below
                elif len( self.line_tokens[ for_body.line_number ] ) != 1 :
                    # ok, the line is not going to be empty after we remove the "{"
                    self.excise_token_and_move_upwards_and_after( for_body, end_paren )
                else :
                    # print('ok, excise the "{" and delete the line it was on')
                    old_line = for_body.line_number
                    self.excise_token_and_move_upwards_and_after( for_body, end_paren )
                    self.delete_empty_line( old_line )

            right_curly_brace = for_body.children[-1]
            assert( right_curly_brace.spelling == "}" )
            # the right curly brace should be on its own line unless its the sole child
            # or if it's on the same line as the left curly brace
            if len(for_body.children) != 1 :
                last_before_rcb = self.last_descendent( for_body.children[-2] )
                if last_before_rcb.line_number == right_curly_brace.line_number and right_curly_brace.line_number != for_body.line_number :
                    # print('ok! move rcb to its own line')
                    self.excise_right_curly_brace_and_move_to_next_line_on_own( right_curly_brace )
            else :
                # make sure the {} are on the same line together
                rcb = for_body.children[0]
                assert( rcb.spelling == "}" )
                if rcb.line_number != for_body.line_number and not self.intervening_macro_tokens( for_body, rcb ) :
                    self.put_empty_curly_braces_on_line_together( for_body, rcb )
        else :
            # no curly braces
            last_fs_desc = self.last_descendent( for_body )
            last_fd_desc = self.last_descendent( for_dec_tok )
            if last_fs_desc.line_number != last_fd_desc.line_number and not self.intervening_macro_tokens(tok,for_body) :
                # print('we need to add a pair of curly braces')
                lcb_tok = self.add_left_curly_brace_after( end_paren )
                lcb_tok.type = "for-scope"
                self.adjust_token_tree( lcb_tok, tok.children[ tok.children.index( for_dec_tok ) + 1 ], for_body )
                if for_body.line_number == lcb_tok.line_number :
                    # print('move these tokens down to their own line')
                    self.move_tokens_after_to_their_own_next_line( lcb_tok )
                if self.visible_tokens_on_line_after( last_fs_desc ) :
                    print("Possible bug identified on line", last_fs_desc.orig_line_number+1, "of file", self.filename, "5")
                    self.move_tokens_after_to_their_own_next_line( last_fs_desc )
                self.add_right_curly_brace_on_own_line_after( last_fs_desc, lcb_tok )
            else :
                # make sure that there are no other tokens on this line
                if self.line_tokens[ last_fs_desc.line_number ][ -1 ] is not last_fs_desc :
                    if self.visible_tokens_on_line_after( last_fs_desc ) :
                        # print('move these tokens onto their own line')
                        print("Possible bug identified on line", last_fs_desc.orig_line_number+1, "of file", self.filename, "6")
                        self.move_tokens_after_to_their_own_next_line( last_fs_desc )
        self.adjust_lines_for_children( tok ) # recurse

    def update_all_tokens_list_from_line_tokens( self ) :
        count = 0
        for line in self.line_tokens :
            for tok in line :
                self.all_tokens[ count ] = tok
                count += 1

    # go through all the tokens in the tokens list;
    # through this function, the line_tokens array will get updated
    # and tokens moved around.  Token indices will not be correct, but
    # the all_tokens array will be updated.
    # Token line numbers will be correct
    def adjust_lines_as_necessary( self ) :
        self.adjust_lines_for_children( self.toplevel_token )
        self.update_all_tokens_list_from_line_tokens()
        self.renumber_tokens()

    def adjust_lines_for_token( self, token ) :
        if not token.is_visible :
            return
        elif token.type == "for" :
            self.adjust_for( token )
        elif token.type == "if" :
            self.adjust_if( token )
        elif token.type == "else" :
            self.adjust_else( token )
        else :
            self.adjust_lines_for_children( token )

    def adjust_lines_for_children( self, token ) :
        # iterate across all children, adjust them
        if len( token.children ) == 0 : return
        child_iter = token.children[ 0 ]
        ind = 0
        while True :
            self.adjust_lines_for_token( child_iter )

            # find the next token to adjust
            if ind < len(token.children) and child_iter is token.children[ind] :
                new_ind = ind+1
            else :
                new_ind = token.children.index( child_iter )
                new_ind += 1

            if new_ind == len( token.children ) :
                break

            # o/w increment
            ind = new_ind
            child_iter = token.children[ ind ]


    def beautify_code( self ) :
        # first pass: add lines, remove lines, or move code between lines and renumber
        self.adjust_lines_as_necessary()

        # second pass: change spacing between code elements
        self.correct_spacing()

        # final pass: change indentation

        self.new_lines = []
        self.line_indentations = [0] * len(self.line_tokens)
        for line_number, line in enumerate(self.line_tokens) :
            self.determine_indentation_level( line_number )
            new_line = self.line_from_line_tokens( line_number )
            self.new_lines.append( new_line )

    #################### CODE TO TEST IF TWO FILES ARE EQUIVALENT ###################

    def two_toks( self, other, i_this, i_other ) :
        return self.all_tokens[i_this], other.all_tokens[i_other]

    def two_next_visible( self, other, i_this, i_other ) :
        i_this = self.find_next_visible_token( i_this )
        i_other = other.find_next_visible_token( i_other )
        return i_this, i_other

    def equiv_if_both_empty_or_if_neither_empty_and_equiv( self, other, i_this, i_other, stack ) :
        if debug_equiv: self.enter_equiv( other, "equiv_if_both_empty_or_if_neither_empty_and_equiv", i_this, i_other, stack )
        tok_this, tok_other = self.two_toks( other, i_this, i_other )
        if tok_this.type == "empty-statement" or tok_other == "empty-statement" :
            if tok_this.equivalent( tok_other ) :
                return True, i_this+1, i_other+1
            elif tok_this.type == "empty-statement" and tok_other.type in self.scope_types :
                return self.equiv_to_scoped_statement( other, i_this, i_other, stack, self.empty_statement_equiv )
            elif tok_this.type in self.scope_types :
                retval = other.equiv_to_scoped_statement( self, i_other, i_this, stack, other.empty_statement_equiv )
                return retval[0], retval[2], retval[1]
            else :
                #print("tok_this", tok_this.spelling, tok_this.line_number+1, tok_this.type, \)
                #    "tok_other", tok_other.spelling, tok_other.line_number+1, tok_other.type
                return False, i_this, i_other
        else :
            #print("tok_this:", tok_this.spelling,"tok_other:", tok_other.spelling)
            return self.statement_equiv( other, i_this, i_other, stack )

    def equiv_to_scoped_statement( self, other, i_this, i_other, stack, func ) :
        assert( other.all_tokens[i_other].type in other.scope_types )
        # this might be equivalent if tok_other contains only a single statement
        # print(func)
        still_good, i_this, i_other = func( other, i_this, other.find_next_visible_token(i_other+1), stack )
        if not still_good : return still_good, i_this, i_other
        i_other = other.find_next_visible_token( i_other )
        tok_other = other.all_tokens[i_other]
        if tok_other.parents_last_child() :
            return True, i_this, i_other+1
        else :
            return False, i_this, i_other


    def for_equiv( self, other, i_this, i_other, stack ) :
        if debug_equiv: self.enter_equiv( other, "for_equiv", i_this, i_other, stack )
        tok_this, tok_other = self.two_toks( other, i_this, i_other )
        assert( tok_this.type == "for" )
        if tok_other.type == "empty-statement" : return self.for_equiv( other, i_this, other.find_next_visible_token( i_other+1 ), stack )
        if tok_other.type in self.scope_types :
            return self.equiv_to_scoped_statement( other, i_this, i_other, stack, self.for_equiv )
        if tok_other.type != "for" :
            return False, i_this, i_other

        stack.append( "for_equiv" )
        i_this += 1
        i_other += 1
        still_good, i_this, i_other = self.paren_block_equiv( other, i_this, i_other, "for", "for-declaration", stack )
        if not still_good :
            stack.pop("paren block not equiv")
            return still_good, i_this, i_other
        retval = self.equiv_if_both_empty_or_if_neither_empty_and_equiv( other, i_this, i_other, stack )
        stack.pop()
        return retval

    def if_equiv( self, other, i_this, i_other, stack ) :
        if debug_equiv: self.enter_equiv( other, "if_equiv", i_this, i_other, stack )
        tok_this, tok_other = self.two_toks( other, i_this, i_other )
        assert( tok_this.type == "if" )
        if tok_other.type == "empty-statement" : return self.if_equiv( other, i_this, other.find_next_visible_token( i_other+1 ), stack )
        if tok_other.type in self.scope_types :
            return self.equiv_to_scoped_statement( other, i_this, i_other, stack, self.if_equiv )
        if tok_other.type != "if" : return False, i_this, i_other

        stack.append( "if_equiv" )
        still_good, i_this, i_other = self.paren_block_equiv( other, i_this, i_other, "if", "if-condition", stack )
        if not still_good : stack.pop(); return still_good, i_this, i_other
        i_this, i_other = self.two_next_visible( other, i_this, i_other )
        still_good, i_this, i_other = self.equiv_if_both_empty_or_if_neither_empty_and_equiv( other, i_this, i_other, stack )
        if not still_good : stack.pop("statements not equiv"); return still_good, i_this, i_other

        # ok -- let's go fishing for an else
        i_this2, i_other2 = self.two_next_visible( other, i_this, i_other )
        tok_this, tok_other = self.two_toks( other, i_this2, i_other2 )
        if tok_this.type == "else" or tok_other.type == "else" :
            if tok_this.type != tok_other.type : stack.pop("type mismatch"); return False, i_this2, i_other2
            i_this, i_other = self.two_next_visible( other, i_this2+1, i_other2+1 )
            retval = self.equiv_if_both_empty_or_if_neither_empty_and_equiv( other, i_this, i_other, stack )
            stack.pop("following else")
            return retval
        else :
            stack.pop("else not found " + tok_this.spelling + " " + tok_other.spelling )
            return True, i_this, i_other

    def paren_block_equiv( self, other, i_this, i_other, end_context, left_paren_type, stack ) :
        if debug_equiv: self.enter_equiv( other, "paren_block_equiv", i_this, i_other, stack )
        stack.append( "paren_block_equiv" )
        found_open_paren = False
        while i_this < len(self.all_tokens) and i_other < len(other.all_tokens) :
            tok_this = self.all_tokens[i_this]
            tok_other = other.all_tokens[i_other]
            if tok_this.context() == end_context and tok_this.type != left_paren_type and found_open_paren :
                stack.pop()
                if debug_equiv: print(" "*len(stack) + "leaving paren_block_equiv", i_this, i_other)
                return True, i_this, i_other
            if tok_this.type == left_paren_type :
                found_open_paren = True
            if not tok_this.equivalent( tok_other ) :
                stack.pop()
                return False, i_this, i_other
            i_this, i_other = self.two_next_visible( other, i_this+1, i_other+1 )

        print("Ran out of tokens in paren_block_equiv, end_context =", end_context,"left_paren_type =", left_paren_type)
        stack.pop()
        return False, i_this, i_other

    def do_while_equiv( self, other, i_this, i_other, stack ) :
        if debug_equiv: self.enter_equiv( other, "do_while_equiv", i_this, i_other, stack )
        tok_this, tok_other = self.two_toks( other, i_this, i_other )
        assert( tok_this.type == "do-while" )
        if tok_other.type == "empty-statement" : return self.do_while_equiv( other, i_this, other.find_next_visible_token( i_other+1 ), stack )
        if tok_other.type in self.scope_types :
            return self.equiv_to_scoped_statement( other, i_this, i_other, stack, self.do_while_equiv )
        if tok_other.type != "do-while" : return False, i_this, i_other

        stack.append( "do_while_equiv" )
        i_this, i_other = self.two_next_visible( other, i_this+1, i_other+1 )

        still_good, i_this, i_other = self.equiv_if_both_empty_or_if_neither_empty_and_equiv( other, i_this, i_other, stack )
        if not still_good : stack.pop(); return still_good, i_this, i_other


        still_good, i_this, i_other = self.paren_block_equiv( other, i_this, i_other, "do-while", "do-while-condition", stack )
        if not still_good : stack.pop(); return still_good, i_this, i_other

        i_this, i_other = self.two_next_visible( other, i_this, i_other )
        tok_this, tok_other = self.two_toks( other, i_this, i_other )

        # print("leaving do_while_equiv", i_this, i_other,)
        if tok_this.spelling == ";" and tok_other.spelling == ";" :
            # print("returning true")
            stack.pop()
            return True, i_this+1, i_other+1
        else :
            # print("returning false")
            stack.pop()
            return False, i_this, i_other

    def while_equiv( self, other, i_this, i_other, stack ) :
        if debug_equiv: self.enter_equiv( other, "while_equiv", i_this, i_other, stack )
        tok_this, tok_other = self.two_toks( other, i_this, i_other )
        assert( tok_this.type == "while" )
        if tok_other.type == "empty-statement" : return self.while_equiv( other, i_this, other.find_next_visible_token( i_other+1 ), stack )

        if tok_other.type in self.scope_types :
            return self.equiv_to_scoped_statement( other, i_this, i_other, stack, self.while_equiv )
        if tok_other.type != "while" : return False, i_this, i_other
        stack.append( "while_equiv" )
        i_this, i_other = self.two_next_visible( other, i_this+1, i_other+1 )
        still_good, i_this, i_other = self.paren_block_equiv( other, i_this, i_other, "while", "while-condition", stack )
        if not still_good : stack.pop(); return still_good, i_this, i_other

        retval = self.equiv_if_both_empty_or_if_neither_empty_and_equiv( other, i_this, i_other, stack )
        stack.pop()
        return retval

    def namespace_equiv( self, other, i_this, i_other, stack ) :
        if debug_equiv: self.enter_equiv( other, "namespace_equiv", i_this, i_other, stack )
        tok_this, tok_other = self.two_toks( other, i_this, i_other )
        if tok_other.type == "empty-statement" : return self.namespace_equiv( other, i_this, other.find_next_visible_token( i_other+1 ), stack )
        if not tok_this.equivalent( tok_other ) : return False, i_this, i_other
        i_this, i_other = self.two_next_visible( other, i_this+1, i_other+1 )
        while i_this < len(self.all_tokens) and i_other < len(self.all_tokens) :
            tok_this2, tok_other2 = self.two_toks( other, i_this, i_other )
            if tok_this2.type == "namespace-scope" or tok_other2 == "namespace-scope" :
                if tok_this2.type != tok_other2.type :
                    return False, i_this, i_other
                else :
                    stack.append( "namespace_equiv" )
                    retval = self.statement_equiv( other, i_this, i_other, stack )
                    stack.pop()
                    return retval
            if tok_this2.parent is not tok_this or tok_other2.parent is not tok_other :
                if debug_equiv: print(" "*len(stack) + "leaving snamespace_equiv", tok_this2.spelling, tok_this2.parent.spelling, tok_other2.spelling, tok_other2.parent.spelling, i_this, i_other,)
                if tok_this2.parent is not tok_this and tok_other2.parent is not tok_other :
                    if debug_equiv: print("returning true")
                    return True, i_this, i_other
                else :
                    if debug_equiv: print("returning false")
                    return False, i_this, i_other
            if not tok_this2.equivalent( tok_other2 ) :
                return False, i_this, i_other
            i_this, i_other = self.two_next_visible( other, i_this+1, i_other+1 )

        print("Ran out of tokens in namespace_equiv")
        return False, i_this, i_other

    def class_equiv( self, other, i_this, i_other, stack ) :
        if debug_equiv: self.enter_equiv( other, "class_equiv", i_this, i_other, stack )
        tok_this, tok_other = self.two_toks( other, i_this, i_other )
        if tok_other.type == "empty-statement" : return self.class_equiv( other, i_this, other.find_next_visible_token( i_other+1 ), stack )
        if not tok_this.equivalent( tok_other ) : return False, i_this, i_other

        stack.append( "class_equiv" )
        i_this, i_other = self.two_next_visible( other, i_this+1, i_other+1 )
        while i_this < len(self.all_tokens) and i_other < len(other.all_tokens) :
            tok_this, tok_other = self.two_toks( other, i_this, i_other )
            if tok_this.is_inside_string or tok_other.is_inside_string :
                pass
            elif tok_this.type == "class-scope" or tok_other.type == "class-scope" :
                if not tok_this.equivalent( tok_other ) :
                    stack.pop()
                    return False, i_this, i_other
                else :
                    still_good, i_this, i_other = self.statement_equiv( other, i_this, i_other, stack )
                    if not still_good : stack.pop(); return still_good, i_this, i_other
                    continue
            elif tok_this.spelling == ";" or tok_other.spelling == ";" :
                stack.pop()
                if not tok_this.equivalent( tok_other ) :
                    return False, i_this, i_other
                else :
                    return True, i_this+1, i_other+1

            if not tok_this.equivalent( tok_other ) :
                stack.pop()
                return False, i_this, i_other
            i_this, i_other = self.two_next_visible( other, i_this+1, i_other+1 )

        stack.pop()
        return i_this == len(self.all_tokens) and i_other == len(other.all_tokens), i_this, i_other

    def class_privacy_equiv( self, other, i_this, i_other, stack ) :
        if debug_equiv: self.enter_equiv( other, "class_privacy_equiv", i_this, i_other, stack )
        tok_this, tok_other = self.two_toks( other, i_this, i_other )
        if tok_other.type == "empty-statement" : return self.class_privacy_equiv( other, i_this, other.find_next_visible_token(i_other+1), stack)
        if not tok_this.equivalent( tok_other ) : return False, i_this, i_other

        stack.append( "class_privacy_equiv" )
        i_this, i_other = self.two_next_visible(other,i_this, i_other)
        while i_this < len(self.all_tokens) and i_other < len(other.all_tokens) :
            tok_this, tok_other = self.two_toks(other,i_this,i_other)
            if tok_this.spelling == ":" or tok_other.spelling == ":" :
                stack.pop()
                if tok_this.equivalent(tok_other) :
                    return True, i_this+1, i_other+1
                else :
                    return False, i_this, i_other
            i_this, i_other = self.two_next_visible(other,i_this+1,i_other+1)

        stack.pop()
        return i_this == len(self.all_tokens) and i_other == len(other.all_tokens), i_this, i_other

    def union_equiv( self, other, i_this, i_other, stack ) :
        if debug_equiv: self.enter_equiv( other, "union_equiv", i_this, i_other, stack )
        tok_this, tok_other = self.two_toks( other, i_this, i_other )
        if tok_other.type == "empty-statement" : return self.union_equiv( other, i_this, other.find_next_visible_token( i_other+1 ), stack )
        if not tok_this.equivalent( tok_other ) : return False, i_this, i_other

        stack.append( "union_equiv" )
        i_this, i_other = self.two_next_visible( other, i_this+1, i_other+1 )
        while i_this < len(self.all_tokens) and i_other < len(other.all_tokens) :
            tok_this, tok_other = self.two_toks( other, i_this, i_other )
            if tok_this.type == "union-scope" or tok_other.type == "union-scope" :
                if tok_this.type != tok_other.type :
                    stack.pop()
                    return False, i_this, i_other
                else :
                    still_good, i_this, i_other = self.statement_equiv( other, i_this, i_other, stack )
                    if not still_good : stack.pop(); return still_good, i_this, i_other
                    i_this, i_other = self.two_next_visible( other, i_this, i_other )
            if tok_this.spelling == ";" or tok_other.spelling == ";" :
                stack.pop()
                if tok_this.spelling != tok_other.spelling :
                    return False, i_this, i_other
                else :
                    return True, i_this+1, i_other+1
            if not tok_this.equivalent( tok_other ) :
                stack.pop()
                return False, i_this, i_other
            i_this, i_other = self.two_next_visible( other, i_this+1, i_other+1 )

        # print("Ran out of tokens in union_equiv")
        stack.pop()
        return i_this == len(self.all_tokens) and i_other == len(other.all_tokens), i_this, i_other

    def switch_equiv( self, other, i_this, i_other, stack ) :
        if debug_equiv: self.enter_equiv( other, "switch_equiv", i_this, i_other, stack )
        tok_this, tok_other = self.two_toks( other, i_this, i_other )
        if tok_other.type == "empty-statement" : return self.switch_equiv( other, i_this, other.find_next_visible_token( i_other+1 ), stack )
        if tok_other.type in self.scope_types :
            return self.equiv_to_scoped_statement( other, i_this, i_other, stack, self.switch_equiv )
        if not tok_this.equivalent( tok_other ) : return False, i_this, i_other

        stack.append( "switch_equiv" )
        i_this, i_other = self.two_next_visible( other, i_this+1, i_other+1 )
        while i_this < len(self.all_tokens) and i_other < len(other.all_tokens) :
            tok_this, tok_other = self.two_toks( other, i_this, i_other )
            if tok_this.type == "switch-scope" or tok_other.type == "switch-scope" :
                if not tok_this.equivalent( tok_other ) :
                    stack.pop()
                    return False, i_this, i_other
                else :
                    retval = self.switch_scope_equiv( other, i_this, i_other, stack )
                    stack.pop()
                    return retval
            if not tok_this.equivalent( tok_other ) :
                stack.pop()
                return False, i_this, i_other
            i_this, i_other = self.two_next_visible( other, i_this+1, i_other+1 )

        print("Ran out of tokens in switch_equiv")
        stack.pop()
        return False, i_this, i_other

    def switch_scope_equiv( self, other, i_this, i_other, stack ) :
        if debug_equiv: self.enter_equiv( other, "switch_scope_equiv", i_this, i_other, stack )
        tok_this, tok_other = self.two_toks( other, i_this, i_other )
        if not tok_this.equivalent( tok_other ) : return False, i_this, i_other

        stack.append( "switch_scope_equiv" )
        i_this, i_other = self.two_next_visible( other, i_this+1, i_other+1 )
        while i_this < len( self.all_tokens ) and i_other < len( other.all_tokens ) :
            tok_this, tok_other = self.two_toks( other, i_this, i_other )
            if tok_this.type == "case" or tok_other.type == "case" :
                if tok_this.type != tok_other.type : stack.pop( "cases not equal" ); return False, i_this, i_other
                still_good, i_this, i_other = self.case_equiv( other, i_this, i_other, stack )
                if not still_good : stack.pop("case not equiv"); return still_good, i_this, i_other
                continue
            if ( tok_this.parents_last_child() and tok_this.regular_rcb() ) or ( tok_other.parents_last_child() and tok_other.regular_rcb() ) :
                stack.pop("last child" + tok_this.spelling + " " + tok_other.spelling )
                if tok_this.equivalent( tok_other ) : return True, i_this+1, i_other+1
                else : return False, i_this, i_other

            if not tok_this.equivalent( tok_other ) : stack.pop(); return False, i_this, i_other
            i_this, i_other = self.two_next_visible( other, i_this+1, i_other+1 )

        print("Ran out of tokens in switch_scope_equiv")
        stack.pop()
        return False, i_this, i_other

    def case_equiv( self, other, i_this, i_other, stack ) :
        if debug_equiv: self.enter_equiv( other, "case_equiv", i_this, i_other, stack )
        tok_this, tok_other = self.two_toks( other, i_this, i_other )
        if tok_other.type in self.scope_types :
            return self.equiv_to_scoped_statement( other, i_this, i_other, stack, self.case_equiv )
        if not tok_this.equivalent( tok_other ) : return False, i_this, i_other

        stack.append( "case_equiv" )
        i_this, i_other = self.two_next_visible( other, i_this+1, i_other+1 )
        while i_this < len(self.all_tokens) and i_other < len(other.all_tokens) :
            tok_this, tok_other = self.two_toks( other, i_this, i_other )
            if tok_this.context() == "case" or tok_other.context() == "case" :
                if not tok_this.equivalent( tok_other ) : stack.pop(); return False, i_this, i_other
                i_this, i_other = self.two_next_visible( other, i_this+1, i_other+1 )
            elif tok_this.type == "empty-statement" or tok_other.type == "empty-statement" :
                if tok_this.type == "empty-statement" : i_this += 1
                if tok_other.type == "empty-statement" : i_other += 1
                i_this, i_other = self.two_next_visible( other, i_this, i_other )
                continue
            elif tok_this.type == "case" or tok_other.type == "case" :
                stack.pop()
                return True, i_this, i_other
            elif (tok_this.context() == "switch-scope" and tok_this.parents_last_child()) or \
                    (tok_other.context() == "switch-scope" and tok_other.parents_last_child()) :
                stack.pop()
                if tok_this.spelling != tok_other.spelling : return False, i_this, i_other
                return True, i_this, i_other
            else :
                still_good, i_this, i_other = self.statement_equiv( other, i_this, i_other, stack )
                if not still_good : stack.pop(); return still_good, i_this, i_other
                i_this, i_other = self.two_next_visible( other, i_this, i_other )

        print("Ran out of tokens in case_equiv")
        stack.pop()
        return False, i_this, i_other

    def scope_equiv( self, other, i_this, i_other, stack ) :
        if debug_equiv: self.enter_equiv( other, "scope_equiv", i_this, i_other, stack )
        tok_this, tok_other = self.two_toks( other, i_this, i_other )
        if tok_other.type == "empty-statement" : return self.scope_equiv( other, i_this, other.find_next_visible_token( i_other+1 ), stack )

        if tok_other.type not in self.scope_types :
            # this might be equivalent if tok_this contains only a single statement
            still_good, i_this, i_other = self.statement_equiv( other, self.find_next_visible_token( i_this+1), i_other, stack )
            if not still_good : return still_good, i_this, i_other
            i_this = self.find_next_visible_token( i_this )
            tok_this = self.all_tokens[i_this]
            if tok_this.parents_last_child() :
                return True, i_this+1, i_other
            else :
                return False, i_this, i_other

        if not tok_this.equivalent( tok_other ) : return False, i_this, i_other
        stack.append( "scope_equiv" )
        i_this, i_other = self.two_next_visible( other, i_this+1, i_other+1 )
        while i_this < len(self.all_tokens) and i_other < len(other.all_tokens) :
            tok_this, tok_other = self.two_toks( other, i_this, i_other )
            if debug_equiv: print(" "*len(stack)+"scope equiv while loop", tok_this.spelling, tok_other.spelling)
            if tok_this.type == "empty-statement" or tok_other.type == "empty-statement" :
                if tok_this.type == "empty-statement" : i_this += 1
                if tok_other.type == "empty-statement" : i_other += 1
            elif tok_this.parents_last_child() or tok_other.parents_last_child() :
                # is this the last child?  If so, it's the rcb
                if debug_equiv: print(" "*len(stack) + "leaving scope_equiv", i_this, i_other)
                stack.pop()
                if tok_this.spelling != tok_other.spelling : return False, i_this, i_other
                return True, i_this+1, i_other+1
            else :
                still_good, i_this, i_other = self.statement_equiv( other, i_this, i_other, stack )
                if debug_equiv: print(" "*len(stack) + "scope equiv returned from self.statement_equiv", i_this, i_other)
                if not still_good : stack.pop(); return still_good, i_this, i_other
            #print("scope equiv while loop", i_this, i_other)
            i_this, i_other = self.two_next_visible( other, i_this, i_other )

        # print("Ran out of tokens in scope_equiv" -- this isn't an error condition!)
        stack.pop()
        return i_this == len(self.all_tokens) and i_other == len(other.all_tokens), i_this, i_other

    def empty_statement_equiv( self, other, i_this, i_other, stack ) :
        if debug_equiv: self.enter_equiv( other, "empty_statement_equiv", i_this, i_other, stack )
        tok_this, tok_other = self.two_toks( other, i_this, i_other )
        assert( tok_this.type == "empty-statement" )
        if tok_other.type == "empty-statement" : return True, i_this+1, i_other+1
        elif tok_other.type in self.scope_types :
            return self.equiv_to_scoped_statement( other, i_this, i_other, stack, self.empty_statement_equiv )
        else :
            return self.statement_equiv( other, self.find_next_visible_token(i_this+1), i_other, stack )

    def simple_statement_equiv( self, other, i_this, i_other, stack ) :
        if debug_equiv: self.enter_equiv( other, "simple_statement_equiv", i_this, i_other, stack )
        tok_this, tok_other = self.two_toks( other, i_this, i_other )
        if tok_other.type == "empty-statement" : return self.simple_statement_equiv( other, i_this, other.find_next_visible_token(i_other+1), stack )

        if tok_other.type in self.scope_types :
            return self.equiv_to_scoped_statement( other, i_this, i_other, stack, self.simple_statement_equiv )
        if not tok_this.equivalent( tok_other ) : return False, i_this, i_other
        i_this, i_other = self.two_next_visible( other, i_this+1, i_other+1 )

        while i_this < len( self.all_tokens ) and i_other < len( other.all_tokens ) :
            tok_this2, tok_other2 = self.two_toks( other, i_this, i_other )
            if debug_equiv: print(" "*len(stack)+"simple statement equiv while loop", tok_this2.spelling, tok_other2.spelling)
            if tok_this2.parent is not tok_this or tok_other2.parent is not tok_other :
                if debug_equiv: print(" "*len(stack) + "leaving simple_statement_equiv", tok_this2.spelling, tok_this2.parent.spelling, tok_other2.spelling, tok_other2.parent.spelling, i_this, i_other,)
                if tok_this2.parent is not tok_this and tok_other2.parent is not tok_other :
                    if debug_equiv: print("returning true")
                    return True, i_this, i_other
                else :
                    if debug_equiv: print("returning false")
                    return False, i_this, i_other
            if tok_this2.type in self.scope_types :
                still_good, i_this, i_other = self.equiv_if_both_empty_or_if_neither_empty_and_equiv( other, i_this, i_other, stack )
                if not still_good:
                    print
                    print("not still good from scope within simple_statement_equiv")
                    return still_good, i_this, i_other
                continue # do not increment i_this and i_other
            if not tok_this2.equivalent( tok_other2 ) :
                if debug_equiv: print(" " *len(stack) + "toks not equivalent:'" + tok_this2.spelling + "' vs '" + tok_other2.spelling + "'")
                return False, i_this, i_other
            i_this, i_other = self.two_next_visible( other, i_this+1, i_other+1 )
        #print("simple_statement_equiv should not have reached here", i_this, i_other)
        return True, i_this, i_other

    def substatement_equiv( self, other, i_this, i_other, stack ) :
        if debug_equiv: self.enter_equiv( other, "substatement_equiv", i_this, i_other, stack )
        tok_this, tok_other = self.two_toks( other, i_this, i_other )

        #if tok_other.type == "empty-statement" : return self.simple_statement_equiv( other, i_this, other.find_next_visible_token(i_other+1), stack )

        #if tok_other.type in self.scope_types :
        #    return self.equiv_to_scoped_statement( other, i_this, i_other, stack, self.simple_statement_equiv )

        if not tok_this.equivalent( tok_other ) : return False, i_this, i_other
        i_this, i_other = self.two_next_visible( other, i_this+1, i_other+1 )

        while i_this < len( self.all_tokens ) and i_other < len( other.all_tokens ) :
            tok_this2, tok_other2 = self.two_toks( other, i_this, i_other )
            if tok_this2.parent is not tok_this or tok_other2.parent is not tok_other :
                #print("Finished substatement_equiv:", tok_ths2.spelling, tok_this2.parent.spelling, tok_other2.spelling,  tok_other2.parent.spelling)
                if debug_equiv: print(" "*len(stack) + "leaving substatement_equiv", tok_this2.spelling, tok_this2.parent.spelling, tok_other2.spelling, tok_other2.parent.spelling, i_this, i_other,)
                if tok_this2.parent is not tok_this and tok_other2.parent is not tok_other :
                    if debug_equiv: print("returning true")
                    return True, i_this, i_other
                else :
                    if debug_equiv: print("returning false")
                    return False, i_this, i_other
            if tok_this2.type in self.scope_types :
                still_good, i_this, i_other = self.equiv_if_both_empty_or_if_neither_empty_and_equiv( other, i_this, i_other, stack )
                #still_good, i_this, i_other = self.equiv_to_scoped_statement( other, i_this, i_other, stack, self.simple_statement_equiv )
                if not still_good:
                    return still_good, i_this, i_other
            elif not tok_this2.equivalent( tok_other2 ) :
                if debug_equiv: print(" " *len(stack) + "toks not equivalent:'" + tok_this2.spelling + "' vs '" + tok_other2.spelling + "'")
                return False, i_this, i_other
            i_this, i_other = self.two_next_visible( other, i_this+1, i_other+1 )
        #print("simple_statement_equiv should not have reached here", i_this, i_other)
        return True, i_this, i_other

    def template_equiv( self, other, i_this, i_other, stack ) :
        if debug_equiv: self.enter_equiv( other, "template_equiv", i_this, i_other, stack )
        tok_this, tok_other = self.two_toks( other, i_this, i_other )
        if tok_other.type == "empty-statement" : return self.function_equiv( other, i_this, other.find_next_visible_token(i_other+1), stack)
        if not tok_this.equivalent( tok_other ) : return False, i_this, i_other

        i_this, i_other = self.two_next_visible( other, i_this+1, i_other+1 )

        while i_this < len(self.all_tokens) and i_other < len(other.all_tokens) :
            if debug_equiv: print("template_equiv while loop:", i_this, i_other)
            tok_this, tok_other = self.two_toks( other, i_this, i_other )
            if tok_this.spelling == ">" or tok_other.spelling == ">" :
                if tok_this.spelling == ">" and tok_this.parents_last_child() :
                    if tok_this.equivalent(tok_other) and tok_other.parents_last_child() :
                        return True, i_this+1, i_other+1
                    return False, i_this, i_other

            if not tok_this.equivalent( tok_other ) : return False, i_this, i_other
            i_this, i_other = self.two_next_visible( other, i_this+1, i_other+1 )

        return i_this == len(self.all_tokens) and i_other == len(other.all_tokens), i_this, i_other

    def function_equiv( self, other, i_this, i_other, stack ) :
        if debug_equiv: self.enter_equiv( other, "function_equiv", i_this, i_other, stack )
        tok_this, tok_other = self.two_toks( other, i_this, i_other )
        if tok_other.type == "empty-statement" : return self.function_equiv( other, i_this, other.find_next_visible_token(i_other+1), stack )
        if not tok_this.equivalent( tok_other ) : return False, i_this, i_other

        stack.append( "function_equiv" )
        i_this, i_other = self.two_next_visible( other, i_this+1, i_other+1 )

        while i_this < len(self.all_tokens) and i_other < len(other.all_tokens) :
            # if debug_equiv: print("function_equiv while loop:", i_this, i_other)
            tok_this, tok_other = self.two_toks( other, i_this, i_other )
            if tok_this.type == "function-scope" or tok_other.spelling == "function-scope" :
                retval = self.statement_equiv( other, i_this, i_other, stack )
                stack.pop()
                return retval
            elif ( not tok_this.is_inside_string and tok_this.spelling == ";" ) or \
                    ( not tok_other.is_inside_string and tok_other.spelling == ";" ) :
                stack.pop()
                if tok_this.spelling == tok_other.spelling : return True, i_this+1, i_other+1
                else : return False, i_this, i_other
            else :
                if not tok_this.equivalent( tok_other ) : stack.pop(); return False, i_this, i_other
            i_this, i_other = self.two_next_visible( other, i_this+1, i_other+1 )

        stack.pop()
        return i_this == len(self.all_tokens) and i_other == len(other.all_tokens), i_this, i_other

    def try_equiv( self, other, i_this, i_other, stack ) :
        if debug_equiv: self.enter_equiv( other, "try_equiv", i_this, i_other, stack )
        tok_this, tok_other = self.two_toks( other, i_this, i_other )
        if tok_other.type == "empty-statement" : return self.try_equiv( other, i_this, other.find_next_visible_token(i_other+1), stack )
        if not tok_this.equivalent( tok_other ) : return False, i_this, i_other

        # descend into the try block
        stack.append( "try_equiv" )
        i_this, i_other = self.two_next_visible( other, i_this+1, i_other+1 )
        still_good, i_this, i_other = self.statement_equiv( other, i_this, i_other, stack )
        if not still_good : stack.pop(); return still_good, i_this, i_other

        # descend into catches
        i_this, i_other = self.two_next_visible( other, i_this, i_other )
        last_i_this, last_i_other = i_this, i_other
        while i_this < len(self.all_tokens) and i_other < len( other.all_tokens ) :
            tok_this, tok_other = self.two_toks( other, i_this, i_other )
            if tok_this.type != "catch" and tok_other.type != "catch" :
                stack.pop()
                return True, i_this, i_other
            if not tok_this.equivalent( tok_other ) : stack.pop(); return False, i_this, i_other

            # catch arg
            i_this, i_other = self.two_next_visible( other, i_this+1, i_other+1 )
            still_good, i_this, i_other = self.paren_block_equiv( other, i_this, i_other, "catch", "catch-arg", stack )
            if not still_good : stack.pop(); return still_good, i_this, i_other

            # catch scope
            i_this, i_other = self.two_next_visible( other, i_this, i_other )
            still_good, i_this, i_other = self.statement_equiv( other, i_this, i_other, stack )
            if not still_good : stack.pop(); return still_good, i_this, i_other

    def enum_equiv( self, other, i_this, i_other, stack ) :
        if debug_equiv: self.enter_equiv( other, "enum_equiv", i_this, i_other, stack )
        tok_this, tok_other = self.two_toks( other, i_this, i_other )
        if tok_other.type == "empty-statement" : return self.enum_equiv( other, i_this, other.find_next_visible_token(i_other+1), stack )
        if not tok_this.equivalent( tok_other ) : return False, i_this, i_other

        i_this, i_other = self.two_next_visible( other, i_this+1, i_other+1 )
        while i_this < len(self.all_tokens) and i_other < len(other.all_tokens) :
            tok_this, tok_other = self.two_toks( other, i_this, i_other )
            if tok_this.context() != tok_other.context() : return False, i_this, i_other
            if tok_this.context() != "enum" and tok_this.context() != "enum-scope" : return True, i_this, i_other
            if not tok_this.equivalent( tok_other ) : return False, i_this, i_other
            i_this, i_other = self.two_next_visible( other, i_this+1, i_other+1 )

        return i_this == len(self.all_tokens) and i_other == len(other.all_tokens), i_this, i_other


    def statement_equiv( self, other, i_this, i_other, stack ) :
        if debug_equiv: self.enter_equiv( other, "statement_equiv", i_this, i_other, stack )
        stack.append( "statement_equiv" )
        tok_this = self.all_tokens[ i_this ]
        if tok_this.type == "for" :
            retval = self.for_equiv( other, i_this, i_other, stack )
        elif tok_this.type == "if" :
            retval = self.if_equiv( other, i_this, i_other, stack )
        elif tok_this.type == "do-while" :
            retval = self.do_while_equiv( other, i_this, i_other, stack )
        elif tok_this.type == "while" :
            retval = self.while_equiv( other, i_this, i_other, stack )
        elif tok_this.type == "namespace" :
            retval = self.namespace_equiv( other, i_this, i_other, stack )
        elif tok_this.type == "class" :
            retval = self.class_equiv( other, i_this, i_other, stack )
        elif tok_this.type == "class-privacy" :
            retval = self.class_privacy_equiv( other, i_this, i_other, stack )
        elif tok_this.type == "union" :
            retval = self.union_equiv( other, i_this, i_other, stack )
        elif tok_this.type == "switch" :
            retval = self.switch_equiv( other, i_this, i_other, stack )
        elif tok_this.type == "case" :
            retval = self.case_equiv( other, i_this, i_other, stack )
        elif tok_this.type in self.scope_types :
            retval = self.scope_equiv( other, i_this, i_other, stack )
        elif tok_this.type == "empty-statement" :
            retval = self.empty_statement_equiv( other, i_this, i_other, stack )
        elif tok_this.type == "statement" :
            retval = self.simple_statement_equiv( other, i_this, i_other, stack )
        elif tok_this.type == "substatement" :
            retval = self.substatement_equiv( other, i_this, i_other, stack )
        elif tok_this.type == "template" :
            retval = self.template_equiv( other, i_this, i_other, stack )
        elif tok_this.type == "function" :
            retval = self.function_equiv( other, i_this, i_other, stack )
        elif tok_this.type == "try" :
            retval = self.try_equiv( other, i_this, i_other, stack )
        elif tok_this.type == "enum" :
            retval = self.enum_equiv( other, i_this, i_other, stack )
        else :
            print("We should not have reached here!", tok_this.type, tok_this.line_number+1, tok_this.spelling)
            retval = (False, i_this, i_other)
        stack.pop()
        return retval

    def enter_equiv( self, other, fname, i_this, i_other, stack ) :
        print("%sentering" % (" "*len(stack)), fname, i_this, i_other,)
        print("%s %d " % ( self.all_tokens[i_this].spelling, self.all_tokens[i_this].line_number+1 ) if i_this < len(self.all_tokens) else "None",)
        print("%s %d " % ( other.all_tokens[i_other].spelling, other.all_tokens[i_other].line_number+1 ) if i_other < len(other.all_tokens) else "None")

    def equivalent( self, other ) :
        stack = SmartStack()
        stack.debug = debug_equiv
        i_this = 0
        i_other = 0
        i_this, i_other = self.two_next_visible( other, i_this, i_other )
        while i_this < len( self.all_tokens ) and i_other < len( other.all_tokens ) :
            still_equiv, i_this, i_other = self.statement_equiv( other, i_this, i_other, stack )
            if not still_equiv :
                return False, i_this, i_other
            i_this, i_other = self.two_next_visible( other, i_this, i_other )

        if i_this < len( self.all_tokens ) :
            # make sure there are no visible tokens remaining for self
            i_this = self.find_next_visible_token( i_this )
            if i_this < len( self.all_tokens ) :
                print("unprocessed tokens in self for", self.filename, ":", i_this+1, "of", len(self.all_tokens),"processed")
                print("(other;", i_other, "of", len(other.all_tokens), "processed)")
                #print(", ".join( [ x.spelling for x in self.all_tokens[ i_this: ] ]))
                return False, i_this, i_other
        elif i_other < len( other.all_tokens ) :
            # make sure there are no visible tokens remaining for other
            i_other = self.find_next_visible_token( i_other )
            if i_other < len( other.all_tokens ) :
                print("unprocessed tokens in other")
                return False, i_this, i_other
        return True, i_this, i_other

    #def equivalent( self, other ) :
        # recurse through the tree

def canonical_set_of_macros_to_test() :
    return [ [],
             [ "USE_MPI" ],
             [ "USE_MPI", "SERIALIZATION" ],
             [ "MULTI_THREADED" ],
             [ "USE_OPENMP" ],
             [ "BOINC_GRAPHICS" ],
             [ "WINDOWS" ],
             [ "WIN32" ],
             [ "_WIN32" ],
             [ "MAC" ],
             [ "WIN_PYROSETTA" ] ]



def beautify_file( filename, overwrite, opts = None ) :

    if opts and opts.macro_sets :
        macro_sets = opts.macro_sets
    else :
        macro_sets = canonical_set_of_macros_to_test()


    last_beaut=None
    orig_beaut=None

    orig_lines = open( filename ).readlines()
    if ( any( [ line.endswith('\r\n') for line in orig_lines ] ) ):
       print("WARNING: file %s contains DOS-style line endings. This may cause the beautifier to fail."%filename)
    for macro_set in macro_sets :

        if orig_beaut :
            # print("Macro Inquiries:", ", ".join( orig_beaut.macros_inquired_about ))
            found_any = False
            for macro in macro_set :
                if macro in orig_beaut.macros_inquired_about :
                    found_any = True
                    break
            if not found_any : continue

        beaut = Beautifier()
        beaut.filename = filename
        if opts :
            if opts.pound_if_setting == "take_if" or opts.pound_if_setting == "take_else" :
                #print("Setting beautifier pound_if_setting:", opts.pound_if_setting)
                beaut.pound_if_setting = opts.pound_if_setting

        #print("beautifying", filename, "with macros:" + ", ".join( [ x for x in macro_set ] ))
        for macro in macro_set :
            beaut.defined_macros.append( macro )

        if not last_beaut :
            lines = orig_lines
        else :
            lines = last_beaut.new_lines
        for line in lines :
            beaut.tokenize_line( line )

        beaut.minimally_parse_file()
        beaut.beautify_code()

        if not orig_beaut : orig_beaut = beaut
        last_beaut = beaut

        if debug :
            print("\n\nFINISHED macro set = [", ", ".join( macro_set ), "]")
            for line in beaut.new_lines :
                print(line,)

    if overwrite or debug:
        # make sure that the beautified code is identical to the original code
        all_good = True
        for macro_set in macro_sets :

            if macro_set :
                # print("Macro Inquiries:", ", ".join( orig_beaut.macros_inquired_about ))
                found_any = False
                for macro in macro_set :
                    if macro in orig_beaut.macros_inquired_about :
                        found_any = True
                        break
                if not found_any : continue

            b1 = Beautifier(); b1.filename = filename
            b2 = Beautifier(); b2.filename = filename
            if opts and ( opts.pound_if_setting == "take_if" or opts.pound_if_setting == "take_else" ) :
                b1.pound_if_setting = opts.pound_if_setting
                b2.pound_if_setting = opts.pound_if_setting

            for macro in macro_set :
                b1.defined_macros.append( macro )
                b2.defined_macros.append( macro )
            for line in orig_lines :
                b1.tokenize_line( line )
            for line in beaut.new_lines :
                b2.tokenize_line( line )
            b1.minimally_parse_file()
            b2.minimally_parse_file()
            equiv, t1, t2 = b1.equivalent( b2 )
            if not equiv :
                all_good = False
                print("MACRO set: " + ", ".join( macro_set ) + " did not produce the same tree in the original and beautified code for", filename)

        if overwrite and all_good :
            open( filename, "w" ).writelines( beaut.new_lines )
            #pass
        elif overwrite and not all_good :
            print("Did not beautify", filename, "because tree differed in the presence of some macros")
        elif not overwrite and not all_good :
            print("Beautified file does not match structure of original!")

        if debug and not overwrite :
            open( filename +".beaut", "w" ).writelines( beaut.new_lines )

    elif not overwrite :
        open( filename +".beaut", "w" ).writelines( beaut.new_lines )

class SmartStack :
    def __init__( self ) :
        self.debug = False
        self.stack = []
    def __len__( self ) :
        return len( self.stack )
    def append( self, str ) :
        if self.debug : print(" "*len(self.stack)+"appending", str)
        self.stack.append( str )
    def pop( self, msg = "" ) :
        str = self.stack.pop()
        if self.debug : print(" "*len(self.stack) + "popped", str, msg)



class BeautifierOpts :
    def __init__( self ) :
        self.pound_if_setting = ""
        self.macro_sets = None

if __name__ == "__main__" :
    with blargs.Parser(locals()) as p :
        p.str("filename").required()
        p.flag("overwrite")
        p.str("pound_if_setting").default("take_if") # can be take_if or take_else
        p.multiword( "macros" ).cast( lambda x : x.split() )

    opts = BeautifierOpts()
    if pound_if_setting : opts.pound_if_setting = pound_if_setting
    if macros : opts.macro_sets = [ macros ]

    # print("pound if setting:", opts.pound_if_setting)
    beautify_file( filename, overwrite, opts )
