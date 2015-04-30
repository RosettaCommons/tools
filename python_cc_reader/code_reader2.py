import sys

debug = False

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
                "ctor-initializer-list",
                "block-comment",
                "statement",
                "substatement",
                "template",
                "template-arg-list" ]

class Token :
    def __init__( self ) :
        self.spelling = ""
        self.start = 0
        self.one_past_end = 0
        self.line_number = 0
        self.parent = None
        self.children = []
        self.is_reserve = False
        self.is_commented = False
        self.is_visible = True
        self.invisible_by_macro = False
        self.is_string = False
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

class AdvancedCodeReader :
    def __init__( self ) :
        self.line = ""
        self.all_lines = []
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
        self.reserve_words = set( [ "using", "namespace", "class", "for", "while", "do", \
            "repeat", "public", "private", "protected", "template", "typedef", "typename", \
            "operator", "inline", "explicit", "static", "mutable", "virtual", "friend", \
            "unsigned", "struct", "union" ] )
        self.whitespace = set([" ", "\t", "\n"])
        self.dividers = set([";",":",",","(",")","{","}","=","[","]","<",">","&","|"])
        self.privacy_types = set([ "public", "protected", "private"])
        self.binary_math_symbols = set(["*","-","+","/"])
        self.in_comment_block = False
        self.in_string = False

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
            if self.line[i] in self.whitespace :
                if start >= 0 :
                    self.take_token_for_range( start, i )
                start = -1
                i+=1
                continue
            if self.line[ i ] in self.dividers :
                if start >= 0 :
                    self.take_token_for_range( start, i )
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
                    i+=2
                    continue
            # make block-comment ending symbols their own tokens
            if self.line[ i ] == "*" :
                if i+1 < len( self.line ) and self.line[i+1] == "/" :
                    if start >= 0 :
                        self.take_token_for_range( start, i )
                    self.take_token_for_range( i, i+2 )
                    i+=2
                    continue
            # make string opening or closing symbols their own tokens
            if self.line[ i ] == '"' :
                if i > 0 and ( self.line[i-1] != "\\" ) :
                    if start >= 0 :
                        self.take_token_for_range( start, i )
                    self.take_token_for_range( i,i+1)
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
        self.this_line_tokens.append( tok )

    def process_strings_and_comments( self ) :
        # first check to see if we're in a comment block
        # and then check to see if the first token begins with a #, if it does, then
        # the whole line is invisible, and we have to decide if the next lines
        # are visible or not
        process_as_macro_line = False
        if not self.in_comment_block and len(self.this_line_tokens) > 0 and self.this_line_tokens[0].spelling[0] == "#" :
            process_as_macro_line = True

        for i, tok in enumerate( self.this_line_tokens ) :
            if not self.macro_definitions_say_line_visible() :
                tok.is_visible = False
            if self.in_comment_block :
                tok.is_visible = False
                tok.is_commented = True
                if tok.spelling == "*/" and self.macro_definitions_say_line_visible() :
                    self.in_comment_block = False
            elif self.in_string :
                tok.is_string = True
                if tok.spelling == '"' :
                    self.in_string = False
            elif len(tok.spelling) == 2 and ( tok.spelling == "//" or ( tok.spelling == "/*" and self.macro_definitions_say_line_visible() ) ) :
                # comment begin
                tok.is_visible = False
                tok.is_comment = True
                if tok.spelling == "//" :
                    #the rest of the line is a comment
                    for j in xrange(i,len(self.this_line_tokens)) :
                        self.this_line_tokens[j].is_visible = False
                        self.this_line_tokens[j].is_commented = True
                    return
                if tok.spelling == "/*" :
                    #we're in a comment block
                    self.in_comment_block = True
            elif tok == '"' :
                tok.is_string = True
                self.in_string = True

        if process_as_macro_line : self.process_macro_line()

    def process_macro_line( self ) :
        # this line has been shown to start with "#" so the entire line is being treated as a macro definition
        for tok in self.this_line_tokens :
            tok.is_visible = False
        tok0 = self.this_line_tokens[ 0 ]
        if tok0.spelling == "#include" :
            return
        elif tok0.spelling == "#define" :
            self.defined_macros.add( self.this_line_tokens[1] )
        elif tok0.spelling == "#ifdef" :
            if len(self.this_line_tokens) > 1 :
                name = self.this_line_tokens[1]
                if name in self.defined_macros :
                    # say "true" if the name is in the set of defined macros, and
                    # the parent macro visibility was also "true"
                    self.nested_ifdefs.append( (name, "ifdef", self.nested_ifdefs[-1][2] ) )
                else :
                    self.nested_ifdefs.append( (name, "ifdef", False ))
        elif tok0.spelling == "#ifndef" :
            if len(self.this_line_tokens) > 1 :
                name = self.this_line_tokens[1]
                if name not in self.defined_macros :
                    self.nested_ifdefs.append( (name, "ifndef", self.nested_ifdefs[-1][2] ) )
                else :
                    self.nested_ifdefs.append( (name, "ifndef", False ))
        elif tok0.spelling == "#if" :
            # conservative -- don't try to parse or logically evaluate this line
            # just say "yeah, it's probably true" and if you hit an else later, then
            # instead of negating it, say the else is also true.
            # this could cause problems
            name = " ".join( [ x.spelling for x in self.this_line_tokens[1:] if not x.is_comment ] )
            self.nested_ifdefs.append( (name, "if", self.nested_ifdefs[-1][2] ))
        elif tok0.spelling == "#else" :
            name, iftype, truthval = self.nested_ifdefs.pop()
            if iftype == "else" :
                # we have hit a problem
                print "ERROR: Encountered '#else' inappropriately.  Ifdef stack: "
                for n2, it2, tv2 in reversed( self.nested_ifdefs ) :
                    print "   ", n2, it2, tv2
                assert( iftype != "else" )
            if self.nested_ifdefs[-1][2] :
                if iftype == "ifdef" or iftype == "ifndef" :
                    self.nested_ifdefs.append( (name, "else", not truthval) )
                else :
                    assert( iftype == "if" )
                    self.nested_ifdefs.append( (name, "else", True) )
            else :
                self.netsted_ifdefs.append( (name, "else", "False" ))
        elif tok0.spelling == "#endif" :
            assert( self.nested_ifdefs[-1][1] in set( ["if", "ifdef", "ifndef", "else" ]) )
            self.nested_ifdefs.pop()
        else :
            print "WARNING: Unhandled macro:", tok0.spelling

    def find_next_visible_token( self, i, stack=None ) :
        while i < len(self.all_tokens) :
            if self.all_tokens[i].is_visible : return i
            if stack != None :
                self.set_depth_and_context( i, stack )
            i += 1
        return i
    def print_entry( self, fname, i, stack ) :
        print (" "*len(stack)), fname, i, self.all_tokens[i].spelling, "line number", self.all_tokens[i].line_number+1

    def renumber_tokens( self ) :
        for i,tok in enumerate(self.all_tokens) :
            tok.index = i

    def minimally_parse_file( self ) :
        # at this point, the file has been split into a large set of tokens
        # and we're going to go through the tokens and, to a limited extent,
        # interpret the tokens into a structure
        self.renumber_tokens()
        i = 0
        toplevel_token = Token()
        toplevel_token.type = "namespace-scope"
        stack = [ toplevel_token ] # you start at namespace scope
        while i < len( self.all_tokens ) :
            i = self.process_statement( i, stack )
        if debug : self.print_depth_stack( stack )

    def process_statement( self, i, stack ) :
        if debug : self.print_entry("process_statement",i,stack)
        i = self.find_next_visible_token(i,stack)
        if i == len( self.all_tokens ) : return i
        i_spelling = self.all_tokens[i].spelling
        if i_spelling == "for" :
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
            return self.process_union_decl(i,stack)
        elif i_spelling == "switch" :
            return self.process_switch(i,stack)
        elif i_spelling == "case" or i_spelling == "default" :
            return self.process_case(i,stack)
        elif i_spelling == "{" :
            return self.process_scope(i,stack)
        elif i_spelling == ";" :
            self.set_depth_and_context(i,stack)
            return i+1;
        #elif i_spelling == "}" :
        #    return self.process_scope_end(i,stack)
        elif i_spelling == "using" or i_spelling == "typedef" :
            return self.process_simple_statement(i,stack)
        elif i_spelling == "public" or i_spelling == "protected" or i_spelling == "private" :
            return self.process_privacy_declaration(i,stack)
        elif i_spelling == "template" :
            # just ignore templates, treat them like statements that don't end in ;'s
            i = self.process_template(i,stack)
            return self.process_statement(i,stack)
        elif stack[-1].context() == "namespace-scope" or stack[-1].context() == "class-scope" :
            return self.process_function_preamble_or_variable(i,stack)
        else :
            return self.process_simple_statement(i,stack)

    def set_parent( self, i, stack, token_type = "substatement" ) :
        self.all_tokens[i].type   = token_type
        self.all_tokens[i].parent = stack[-1]
        stack[-1].children.append( self.all_tokens[i] )

    def handle_read_all_tokens_error( self, function_name, stack ) :
        print "ERROR: Ran out of tokens in function", function_name
        self.print_depth_stack( stack )
        sys.exit(1)

    def print_depth_stack( self, stack ) :
        print "Depth stack: "
        for elem in reversed( stack ) :
            print "  ", elem[0], ("\n" if elem[1] == None else "%d %s" % ( elem[1].line_number+1, self.all_lines[ elem[1].line_number ])),

    def read_to_end_paren( self, i, stack ) :
        # the starting "(" should have already been read
        paren_depth = 1
        while i < len( self.all_tokens ) :
            if debug : print (" "*len(stack)), "read to end paren", i, self.all_tokens[i].spelling, self.all_tokens[i].line_number+1
            self.set_parent(i,stack)
            if self.all_tokens[i].spelling == ")" :
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
        template_params = "not-begun"
        while i < len(self.all_tokens) :
            if self.all_tokens[i].spelling == "<" :
                assert( template_params == "not-begun" )
                self.set_parent(i,stack,"template-arg-list")
                template_params = "started"
                continue
            elif self.all_tokens[i].spelling == ">" :
                self.set_parent(i,stack)
                stack.pop() #pop template-arg-list
                stack.pop() #pop template
                return i+1
            self.set_parent(i,stack)
            i+=1
        self.handle_read_all_tokens_error( "process_template", stack )
                
    def process_simple_statement(self,i,stack) :
        if debug : self.print_entry("process_simple_statement",i,stack)

        self.set_parent(i,stack,"statement")
        stack.append(self.all_tokens[i])
        i+=1;
        # read all the tokes from i until the next semicolon
        last = i
        while i < len(self.all_tokens) :
            if self.all_tokens[i].is_visible and self.all_tokens[i].spelling == ";" :
                # print "statement:", " ".join( [ x.spelling for x in self.all_tokens[i:i+1] ] )
                self.set_parent(i,stack)
                stack.pop()
                return i+1
            self.set_parent(i,stack)
            i+=1
        self.handle_read_all_tokens_error( "process_simple_statement", stack )

    def current_classname( self, i, stack ) :
        assert( stack[-1].context() == "class-scope" )
        assert( len(stack) > 2 )
        # let's get the parent "class" token
        class_token = stack[-2][1]
        assert( class_token )
        # and then the next visible token should be the class name
        # class names are never namespace qualified are they? That'd make that annoying
        i = self.find_next_visible_token( class_token.index+1 )
        return self.all_tokens[i].spelling


    def process_function_preamble_or_variable(self,i,stack) :
        # print "process_function_preamble_or_variable", i, depth
        # treat function declarations and variable declarations interchangably.
        # so we should gobble up the tokens that make this a function
        # and see if we have a constructor, in which case, we should also
        # be watchful for an intializer list
        if debug : self.print_entry("process_function_preamble_or_variable",i,stack)
        classname = ""
        if stack[-1][0] == "class-scope" :
            classname = self.current_classname( i, stack )
        arglist = "not-begun"
        is_ctor = False
        found_init_list = False

        i_initial = i # save this in case we aren't actually entering a function; treat it like a statement in that case
        self.set_parent(i,stack,"function")
        stack.append( self.all_tokens[i] )
        i+=1
        while i < len(self.all_tokens) :
            # print (" "*len(stack)), "process function preamble or variable", i, self.all_tokens[i].spelling, self.all_tokens[i].line_number+1
            if not self.all_tokens[i].is_visible :
                # don't do anything special; just don't try and interpret this token as having meaning
                pass
            elif self.all_tokens[i].spelling == ":" :
                # print "found : -- ", classname, is_ctor, arglist, found_init_list
                if is_ctor and arglist == "ended" and not found_init_list :
                    # we're now reading the initializer list
                    # stack.pop() # remove the function-decl-argument-list
                    self.set_parent( i,stack,"ctor-initializer-list")
                    stack.append(self.all_tokens[i])
                    found_init_list = True
                    continue
                elif arglist == "not-begun" :
                    # this might be constructor, we don't know
                    # this might also be some namespace qualifier on a return type
                    if i+1 < len(self.all_tokens) and self.all_tokens[i+1].spelling == ":" :
                        classname = self.all_tokens[i-1].spelling
            elif self.all_tokens[i].spelling == "(" and arglist == "not-begun" :
                funcname = self.all_tokens[i-1].spelling
                if funcname != "operator" :
                    if funcname == classname :
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
                return i

            elif self.all_tokens[i].spelling == "{" :
                if ( found_init_list ) :
                    stack.pop() # remove ctor-initializer-list
                # ok! now descend into the function block
                i = self.process_statement( i, stack )
                stack.pop()
                return i

            self.set_parent(i,stack)
            i+=1
        self.handle_read_all_tokens_error( "process_function_preamble_or_variable", stack )

    def process_for(self,i,stack) :
        if debug : self.print_entry("process_for",i,stack)

        assert( self.all_tokens[i].spelling == "for" )
        self.set_parent(i,stack,"for")
        stack.append(self.all_tokens[i] )
        i+=1
        arglist = "not-yet-begun"
        while i < len(self.all_tokens) :
            if not self.all_tokens[i].is_visible :
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
            if not self.all_tokens[i].is_visible :
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
            if not self.all_tokens[i].is_visible :
                pass
            elif self.all_tokens[i].spelling == "(" :
                self.set_parent(i,stack,"do-while-condition")
                stack.append(self.all_tokens[i])
                i+=1
                i = self.read_to_end_paren( i, stack )
                stack.pop() # pop do-while-condition
                stack.pop() # pop do-while
                return i
            self.set_parent(i,stack)
            i+=1

        self.handle_read_all_tokens_error( "process_do_while", stack )

    def process_while(self,i,stack) :
        assert( self.all_tokens[i].spelling == "while" )
        if debug : self.print_entry("process_while",i,stack)
        self.set_parent(i,stack,"while")
        stack.append(self.all_tokens[i])
        i+=1
        while i < len(self.all_tokens) :
            if not self.all_tokens[i].is_visible :
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
                return i
            self.set_parent(i,stack)
            i+=1
        # we should not have arrived here
        self.handle_read_all_tokens_error( "process_namespace", stack )

    def process_class_decl(self,i,stack) :
        assert( self.all_tokens[i].spelling == "class" )
        if debug : self.print_entry("process_class_decl",i,stack)

        seen_inheritance_colon = False

        self.set_parent(i,stack,"class")
        stack.append(self.all_tokens[i])
        i+=1

        while i < len( self.all_tokens ) :
            if not self.all_tokens[i].is_visible :
                pass
            elif self.all_tokens[i].spelling == "{" :
                if seen_inheritance_colon :
                    stack.pop() # remove class-inheritance-list
                i = self.process_statement( i, stack )
                i = self.find_next_visible_token( i, stack )
                assert( self.all_tokens[i].spelling == ";" )
                self.set_parent( i, stack )
                stack.pop()
                return i+1
            elif self.all_tokens[i].spelling == ":" :
                if not seen_inheritance_colon :
                    seen_inheritance_colon = True
                    self.set_parent(i,stack,"class-inheritance-list")
                    stack.append(self.all_tokens[i])
                    continue
            elif self.all_tokens[i].spelling == ";" :
                # this is just a forward declaration
                self.set_parent(i,stack)
                stack.pop()
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
        while i < len( self.all_tokens ) :
            if not self.all_tokens[i].is_visible :
                pass
            elif self.all_tokens[i].spelling == "{" :
                i = self.process_statement(i,stack)
                # read semi-colon
                i = self.find_next_visible_token(i,stack)
                assert( self.all_tokens[i].spelling == ";" )
                self.set_parent(i,stack)
                stack.pop() # remove union
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
            if not self.all_tokens[i].is_visible :
                pass
            if self.all_tokens[i].spelling == "(" :
                self.set_parent(i,stack,"switch-expression")
                stack.append(self.all_tokens[i])
                i+=1
                i = self.read_to_end_paren( i, stack )
                stack.pop() # remove switch-expression
                i = self.process_statement( i, stack )
                stack.pop() # remove switch
                return i
            self.set_parent(i,stack)
            i+=1
        self.handle_read_all_tokens_error( "process_class_decl", stack )
    def process_case(self,i,stack) :
        assert( self.all_tokens[i] == "case" or self.all_tokens[i] == "default" )
        if debug : self.print_entry("process_case",i,stack)

        self.set_parent(i,stack,"case")
        stack.append(self.all_tokens[i])
        i+=1
        while i < len(self.all_tokens) :
            if not self.all_tokens[i].is_visible :
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

    def process_scope(self,i,stack) :
        if debug : self.print_entry("process_scope",i,stack)

        scopename = stack[-1][0]+"-scope"
        if scopename not in token_types :
            scopename = "scope"
        self.set_parent(i,stack,scopename)
        stack.append(self.all_tokens[i])
        i+=1
        while i < len(self.all_tokens) :
            i = self.find_next_visible_token(i,stack)
            if self.all_tokens[i].spelling == "}" :
                # print (" "*len(stack)), "process_scope", i, "found end }", self.all_tokens[i].line_number+1
                #ok, done processing at this depth
                self.set_parent(i,stack)
                stack.pop() # remove "scope" or "*-scope"
                return i+1
            else :
                # print (" "*len(stack)), "process_scope encountered non-{, non-} at", i, self.all_tokens[i].spelling, self.all_tokens[i].line_number+1
                i = self.process_statement(i,stack)
        self.handle_read_all_tokens_error( "process_scope", stack )

    # def process_scope_end(self,i,stack) :
    #     # print "process_scope_end", i, depth
    #     self.set_parent(i,depth-1,stack.pop())
    #     return i+1, depth-1
    def process_privacy_declaration(self,i,stack) :
        if debug : self.print_entry("process_privacy_declaration",i,stack)
        self.set_parent(i,stack,"class-privacy")
        stack.append(self.all_tokens[i])
        i+=1
        while i < len(self.all_tokens) :
            if not self.all_tokens[i].is_visible :
                pass
            elif self.all_tokens[i].spelling == ":" :
                self.set_parent(i,stack)
                stack.pop()
                return i+1
            self.set_parent(i,stack)
            i+=1
        self.handle_read_all_tokens_error( "process_privacy_declaration", stack )


    def whole_line_invisible( self, line_number ) :
        for tok in self.line_tokens[ line_number ] :
            if tok.is_visible : return False
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
        # print "first token: ", first_token.spelling, "" if not first_token.parent else first_token.parent.spelling

        if first_token.context() == "namespace" :
            # declaring new namespace -- use the indentation level of previous line
            if line_number != 0 :
                self.line_indentations[line_number] = self.line_indentations[line_number-1]
        elif first_token.context() == "namespace-scope" :
            # namespaces scopes don't indent
            if line_number != 0 :
                self.line_indentations[line_number] = self.line_indentations[line_number-1]
        elif first_token.context() == "scope" :
            if first_token.spelling == "}" :
                self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number]
            else :
                # indent one from the parent's indentation
                self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number] + 1
        elif first_token.context() == "for" :
            # this is an odd case -- let's indent twice
            self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number] + 2
        elif first_token.context() == "for-declaration" :
            # this is where where the three statements in a for loop take more than one line
            # indent twice
            self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number] + 2
        elif first_token.context() == "for-scope" :
            if first_token.spelling == "}" :
                self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number]
            else :
                self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number]+1
        elif first_token.context() == "do-while" :
            # this is an odd case -- between the "do" and the "{"
            self.line_indentations[line_number] = self.line_indentions[first_token.parent.line_number]
        elif first_token.context() == "do-while-condition" :
            # the do-while condition has occupied more than one line, I guess.  Indent twice like for and if
            self.line_indentaitons[line_number] = self.line_indentations[first_token.parent.line_number]+2
        elif first_token.context() == "do-while-scope" :
            if first_token.spelling == "}" :
                self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number]
            else :
                self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number]+1
        elif first_token.context() == "if" :
            # this can only happen for the "(" at the beginning of the if-condition or if you have a comment
            # after the if-condition but before the {
            self.line_indentaitons[line_number] = self.line_indentations[first_token.parent.line_number]+2
        elif first_token.context() == "if-condition" :
            # indent twice if the if-condition has gone on so long that it wraps to a second line
            self.line_indentaitons[line_number] = self.line_indentations[first_token.parent.line_number]+2
        elif first_token.context() == "if-scope" :
            if first_token.spelling == "}" :
                self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number]
            else :
                self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number]+1
        elif first_token.context() == "else" :
            # indent to the same level as the if
            self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number]
        elif first_token.context() == "else-scope" :
            if first_token.spelling == "}" :
                self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number]
            else :
                self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number]+1
        elif first_token.context() == "while" :
            # donno, this shouldn't really happen
            self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number]+2
        elif first_token.context() == "while-condition" :
            self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number]+2
        elif first_token.context() == "while-scope" :
            if first_token.spelling == "}" :
                self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number]
            else :
                self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number]+1
        elif first_token.context() == "switch" :
            self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number]+2
        elif first_token.context() == "switch-expression" :
            self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number]+1
        elif first_token.context() == "switch-scope" :
            # don't indent for case: and default: statements, only for case-block and case-block-scopes.
            self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number]
            #if first_token.spelling == "}" :
            #    self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number]
            #else :
            #    self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number]+1
        elif first_token.context() == "case" :
            # dunno -- is this a case condition that went for more than one line?
            self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number]
        elif first_token.context() == "case-block" :
            self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number]+1
            pass
        elif first_token.context() == "case-block-scope" :
            if first_token.spelling == "}" :
                self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number]
            else :
                self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number]+1
        elif first_token.context() == "class" :
            # kind of an odd duck
            # indent once if we're in the process of declaring a class unless we're a "{"
            if self.line_tokens[line_number][0].spelling == "{" :
                self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number]
            else :
                self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number]+1
        elif first_token.context() == "class-inheritance-list" :
            self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number]+1
            pass
        elif first_token.context() == "class-scope" :
            # print "class scope ", line_number, " and parent on ", first_token.parent.line_number
            if first_token.spelling == "}" or first_token.type == "class-privacy" :
                self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number]
            else :
                self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number]+1
        elif first_token.context() == "class-privacy" :
            # don't indent privacy declarations
            self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number]
        elif first_token.context() == "union" :
            pass
        elif first_token.context() == "union-scope" :
            self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number]+1
        elif first_token.context() == "function" :
            # we're in the middle of a function declaration -- e.g. "inline \n void \n etc"
            # or we're declaring a variable -- it's hard to tell
            self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number]
        elif first_token.context() == "function-decl-argument-list" :
            self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number]+1
        elif first_token.context() == "function-scope" :
            if first_token.spelling == "}" :
                self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number]
            else :
                self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number]+1
        elif first_token.context() == "ctor-initializer-list" :
            self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number]+1
        elif first_token.context() == "statement" :
            # this statement has run on to a second line, indent once
            self.line_indentations[line_number] = self.line_indentations[first_token.parent.line_number]+1

    def line_from_line_tokens( self, line_number ) :
        toks = self.line_tokens[line_number]
        spellings = ["\t" * self.line_indentations[line_number] ]
        for i, tok in enumerate( toks ) :
            spellings.append( tok.spelling )
            if i+1 != len( toks ) :
                spellings.append( " " * ( toks[i+1].start - tok.one_past_end ) ) #preserve spaces
        spellings.append("\n")
        return "".join( spellings )

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
            print "Failed to find tok:", tok.spelling, "on line", tok.line_number
            sys.exit(1)
        while i < len(line_toks) :
            line_toks[i].start = line_toks[i].start+1
            line_toks[i].one_past_end = line_toks[i].one_past_end+1
            i+=1

    def enforce_space_between( self, tok1, tok2 ) :
        # print "tok1", tok1.spelling, tok1.line_number, tok1.start, tok1.one_past_end
        # print "tok2", tok2.spelling, tok2.line_number, tok2.start, tok2.one_past_end
        if tok1.line_number == tok2.line_number and tok1.one_past_end == tok2.start :
            self.insert_space_after( tok1 )

    def correct_spacing( self ) :
        for i,tok in enumerate(self.all_tokens) :
            if not tok.is_visible : continue

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
                    self.enforce_space_between( tok, nexttok )
                elif tok.spelling == "else" :
                    # make sure there's a space between the if's "}" and the next "{"
                    assert( prevtok )
                    if prevtok.spelling == "}" :
                        self.enforce_space_between( prevtok, tok )
                    if nexttok.spelling == "{" :
                        self.enforce_space_between( tok, nexttok )
            elif tok.context() == "if-condition" :
                if tok.spelling == ")" and nexttok.context() == "if" :
                    self.enforce_space_between( tok, nexttok )
            elif tok.context() == "else" :
                if tok.spelling == "{" :
                    self.enforce_space_between( tok, nexttok )
            elif tok.context() == "else-scope" :
                if tok.spelling == "}" :
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
                elif tok.spelling == ")" and nexttok.context() == "for" :
                    self.enforce_space_between( prevtok, tok ) # for (..++i_) 
            elif tok.context() == "while" :
                if tok.spelling == "(" :
                    self.enforce_space_between( prevtok, tok ) # while_(
                    self.enforce_space_between( tok, nexttok ) # while (_condition
            elif tok.context() == "while-condition" :
                if tok.spelling == ")" and nexttok.context() == "while" :
                    self.enforce_space_between( prevtok, tok ) # while ( condition_)
                    self.enforce_space_between( tok, nexttok ) # while ( .. )_{
            elif tok.context() == "case" :
                if tok.spelling == ":" and nexttok.context() == "case-block" :
                    self.enforce_space_between( prevtok, tok )
                    self.enforce_space_between( tok, nexttok )
            elif tok.context() == "function" :
                if tok.spelling == "(" and nexttok.spelling != ")" :
                    self.enforce_space_between( tok, nexttok )
            elif tok.context() == "function-decl-argument-list" :
                if tok.spelling == ")" and nexttok.context() == "function" and prevtok.spelling != "(" :
                    self.enforce_space_between( prevtok, tok )
            elif tok.context() == "namespace" :
                if tok.spelling == "{" :
                    self.enforce_space_between( prevtok, tok )

    def last_descendent( self, token ) :
        # return the last descendent of a given token
        # does not require that the index of the token be up-to-date
        if not token.children :
            return token
        else :
            return self.last_descendent( token )

    def excise_left_curly_brace_and_move_after( self, lcb_tok, dest_tok ) :
        # dest_tok == the token after which the left-curly-brace (lcb) is to be placed.
        pass

    def delete_empty_line( self, old_line ) :
        # iterate across all the tokens and decrement the line_number for any that are after
        # this old line number; also delete one line from the line_toks table
        for tok in self.all_tokens) :
            if tok.line_number > old_line : tok.line_number -= 1
        self.line_tokens = self.line_tokens[:old_line] + self.line_tokens[old_line+1:]

    def insert_line_with_single_token( self, tok ) :
        # the token should know what line it belongs on
        self.line_tokens.insert( tok.line_number, [ tok ] )
        for i in xrange( tok.line_number+1, len( self.line_tokens ) ) :
            for tok in self.line_tokens[i] :
                tok.line_number += 1

    def excise_right_curly_brace_and_move_to_next_line( self, rcb ) :
        parent = rcb.parent
        if parent.children[-1] is not rcb :
            # case where there are comments at the end of the line after rcb
            parent.children.remove( rcb )
            parent.children.append( rcb )
        rcb_ind = self.line_tokens[ rcb.line_number ].index( rcb )
        if rcb_ind + 1 != len( self.line_tokens[ rcb.line_number ] ) :
            # ok, we have comments or something at the end of this line
            # if its nothing but comments, then leave them on this line
            # otherwise, they should be on their own line after the }

            any_visible = False
            for i in xrange( rcb_ind+1, len( self.line_tokens[ rcb.line_number ] ) ) :
                if self.line_tokens[ rcb.line_number ][ i ].is_visible :
                    any_visible = True
                    break;
            if any_visible :
                # ok! then we need to move these statements to their own line
                # that's on the to-do list.

    def add_left_curly_brace_after( self, tok_before ) :
        pass

    def add_right_curly_brace_on_own_line_after( self, lcb_tok, tok_before ) :
        rcb_tok = Token()
        rcb_tok.spelling = "}"
        rcb_tok.start = 0 # don't worry about the indentation
        rcb_tok.one_past_end = 1
        rcb_tok.line_number = tok_before.line_number + 1
        self.all_tokens.insert( self.all_tokens.index( tok_before )+1, rcb_tok )
        self.insert_line_with_single_token( rcb_tok )
        rcb_tok.set_parent( lcb_tok )

    def adjust_token_tree( self, new_parent, child ) :
        # adds "new_parent" as the last child of child's former parent
        # and sets new_parent up as the parent of child.
        old_parent = child.parent
        old_parent.children.remove( child )
        old_parent.append( new_parent )
        child.set_parent( new_parent )

    def adjust_for( self, i ) :
        # ok -- let's look at the descendents; they should be
        # the for-declaration descendent and then either a statement of some form
        # or a for-scope
        tok = self.all_tokens[i]
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
                # ok! we need to move the "{" up to the ")" line
                # ask -- are there other tokens on the "{" line?
                if for_body.children[0].spelling == "}" :
                    # ok, keep the "{" and "}" pair on the same line together
                    # to do!
                    assert( False )
                elif for_body.children[0].line_number == for_body.line_number :
                    # ok, the line is not going to be empty after we remove the "{"
                    self.excise_left_curly_brace_and_move_after( for_body, end_paren )
                else :
                    # ok, excise the "{" and delete the line it was on
                    old_line = for_body.line_number
                    self.excise_left_curly_brace_and_move_after( for_body, end_paren )
                    self.delete_empty_line( old_line )

            right_curly_brace = for_body.children[-1]
            assert( right_curly_brace.spelling == "}" )
            # the right curly brace should be on its own line unless its the sole child
            if len(for_body.children) != 1 :
                last_before_rcb = self.last_descendent( for_body.children[-2] )
                if last_before_rcb.line_number == right_curly_brace.line_number :
                    # ok! move rcb to its own line
                    self.excise_right_curly_brace_and_move_to_next_line_on_own( right_curly_brace )
        else :
            last_fs_desc = self.last_descendent( for_body )
            if last_fs_desc.line_number != tak.line_number :
                # we need to add a pair of curly braces
                lcb_tok = self.add_left_curly_brace_after( end_paren )
                lcb_tok.type = "for-scope"
                self.adjust_token_tree( lcb_tok, for_body )
                self.add_right_curly_brace_on_own_line_after( last_fs_desc, lcb_tok )

    # go through all the tokens in the all_tokens list;
    # through this function, the line_tokens array will get updated
    # and tokens moved around.  Token indices will not be correct.
    # Token line numbers will be correct
    def adjust_lines_as_necessary( self ) :
        return
        self.new_lines = []
        self.new_line_tokens = []
        i = 0
        while i < len(self.all_tokens) :
            if not self.all_tokens[i].is_visible : i+=1; continue;
            if self.all_tokens[i].type == "for" :
                i = self.adjust_for( i )
                
            

        self.all_lines = self.new_lines
        self.line_tokens = self.new_line_tokens
        self.renumber_tokens()

    def beautify_code( self ) :
        # first pass: add lines, remove lines, or move code between lines and renumber
        self.adjust_lines_as_necessary()

        # second pass: change spacing between code elements
        self.correct_spacing()

        # final pass: change indentation
        
        self.new_lines = []
        self.line_indentations = [0] * len(self.all_lines)
        for line_number, line in enumerate(self.all_lines) :
            self.determine_indentation_level( line_number )
            new_line = self.line_from_line_tokens( line_number )
            self.new_lines.append( new_line )
