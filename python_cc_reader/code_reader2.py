import sys

contexts = [ "namespace",
             "namespace-scope"
             "scope",
             "for",
             "for-declaration",
             "for-scope"
             "do-while",
             "do-while-condition",
             "do-while-scope",
             "if",
             "if-condition",
             "if-scope",
             "else",
             "else-scope",
             "while",
             "while-condition"
             "while-scope"
             "switch",
             "switch-expression",
             "switch-scope"
             "case",
             "case-scope",
             "default-case",
             "default-case-scope",
             "class",
             "class-inheritance-list",
             "class-scope"
             "class-privacy",
             "struct",
             "struct-inheritance-list",
             "struct-scope"
             "union",
             "union-scope"
             "function",
             "function-decl-argument-list",
             "function-scope",
             "ctor-initializer-list",
             "block-comment",
             "invisible-by-macro" ]

class Token :
    def __init__( self ) :
        self.spelling = ""
        self.start = 0
        self.one_past_end = 0
        self.line_number = 0
        self.filename = ""
        self.is_reserve = False
        self.is_commented = False
        self.is_visible = True
        self.is_string = False
        self.is_macro = False
        self.context = ""
        self.depth = 0

class AdvancedCodeReader :
    def __init__( self ) :
        self.line = ""
        self.all_lines = []
        self.line_number = 0
        self.this_line_tokens = []
        self.all_tokens = [] # tokens for the entire file
        self.line_tokens = [] # the tokens for each line
        self.namespace_stack = []
        self.indentation_depth = 0
        self.context_stack = []
        self.nested_ifdefs  = [ ( "None", "None", True ) ]
        self.defined_macros = []
        self.reserve_words = set( [ "using", "namespace", "class", "for", "while", "do", \
            "repeat", "public", "private", "protected", "template", "typedef", "typename", \
            "operator", "inline", "explicit", "static", "mutable", "virtual", "friend", \
            "unsigned", "struct", "union" ] )
        self.whitespace = set([" ", "\t", "\n"])
        self.dividers = set([";",":",",","(",")","{","}","=","[","]","<",">","&","|"])
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

    def find_next_visible_token( self, i, depth=None, context=None ) :
        while i < len(self.all_tokens) :
            if self.all_tokens[i].is_visible : return i
            if depth != None :
                self.set_depth_and_context( i, depth, context )
            i += 1
        return i

    def minimally_parse_file( self ) :
        # at this point, the file has been split into a large set of tokens
        # and we're going to go through the tokens and, to a limited extent,
        # interpret the tokens into a structure
        i = 0
        stack = [ ("namespace-scope",None) ] # you start at namespace scope
        while i < len( self.all_tokens ) :
            # print "while top", i, depth
            i = self.process_statement( i, stack )
        self.print_depth_stack( stack )

    def process_statement( self, i, stack ) :
        print "process_statement", i, self.all_tokens[i].spelling
        i = self.find_next_visible_token(i,len(stack),stack[-1][0])
        if i > len( self.all_tokens ) : return i
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
            self.set_depth_and_context(i,len(stack),stack)
            return i+1;
        #elif i_spelling == "}" :
        #    return self.process_scope_end(i,stack)
        elif i_spelling == "using" or i_spelling == "typedef" :
            return self.process_simple_statement(i,stack)
        elif i_spelling == "public" or i_spelling == "protected" or i_spelling == "private" :
            return self.process_privacy_declaration(i,stack)
        elif stack[-1][0] == "namespace-scope" or stack[-1][0] == "class-scope" or stack[-1][0] == "struct-scope" :
            return self.process_function_preamble_or_variable(i,stack)
        else :
            return self.process_simple_statement(i,stack)

    def set_depth_and_context( self, i, depth, context ) :
        self.all_tokens[i].depth = depth
        self.all_tokens[i].context = context

    def handle_read_all_tokens_error( self, function_name, stack ) :
        print "ERROR: Ran out of tokens in function", function_name
        self.print_depth_stack( stack )
        sys.exit(1)

    def print_depth_stack( self, stack ) :
        print "Depth stack: "
        for elem in reversed( stack ) :
            print "  ", elem[0], ("\n" if elem[1] == None else "%d %s" % ( elem[1].line_number, self.all_lines[ elem[1].line_number - 1 ])),

    def read_to_end_paren( self, i, stack ) :
        paren_depth = 0
        while i < len( self.all_tokens ) :
            self.set_depth_and_context(i,len(stack),stack[-1][0])
            if self.all_tokens[i].spelling == ")" :
                paren_depth -= 1
                if paren_depth == 0 :
                    return i
            elif self.all_tokens[i].spelling == "(" :
                paren_depth += 1
            i += 1
        self.handle_read_all_tokens_error("process_for", stack)

    def process_simple_statement(self,i,stack) :
        print "process_simple_statement", i, self.all_tokens[i].spelling
        # print "process simple statement", i, len(stack)
        # read all the tokes from i until the next semicolon
        last = i
        while last < len(self.all_tokens) :
            if self.all_tokens[last].is_visible and self.all_tokens[last].spelling == ";" :
                # print "statement:", " ".join( [ x.spelling for x in self.all_tokens[i:last+1] ] )
                break
            last += 1
        if last == len(self.all_tokens) :
            self.handle_read_all_tokens_error( "process_statement", stack );
        for j in xrange(i,last+1) :
            self.set_depth_and_context(j,len(stack),stack[-1][0])
        return last+1

    def process_function_preamble_or_variable(self,i,stack) :
        # print "process_function_preamble_or_variable", i, depth
        # treat function declarations and variable declarations interchangably.
        # so we should gobble up the tokens that make this a function
        # and see if we have a constructor, in which case, we should also
        # be watchful for an intializer list
        print "process_function_preamble_or_variable", i, self.all_tokens[i].spelling
        classname = ""
        arglist = "not-begun"
        is_ctor = False
        found_init_list = False
        stack.append( "function" )
        while i < len(self.all_tokens) :
            if not self.all_tokens[i].is_visible :
                # don't do anything special; just don't try and interpret this token as having meaning
                pass
            elif self.all_tokens[i].spelling == ":" :
                if is_ctor and arglist == "ended" :
                    # we're now reading the initializer list
                    stack.append(( "ctor-initializer-list", self.all_tokens[i] ))
                    found_init_list = True
                elif arglist != "started" :
                    assert( arglist == "not-begun" )
                    # this might be constructor, we don't know
                    # this might also be some namespace qualifier on a return type
                    if i+1 < len(self.all_tokens) and self.all_tokens[i+1].spelling == ":" :
                        classname = self.all_tokens[i-1].spelling
            elif self.all_tokens[i].spelling == "(" :
                funcname = self.all_tokens[i-1].spelling
                if funcname != "operator" :
                    if funcname == classname :
                        is_ctor = True
                    arglist = "started"
                    stack.append(( "function-decl-argument-list", self.all_tokens[i] ))
                    i = self.read_to_end_paren( i, stack )
                    arglist = "ended"
                    stack.pop() # remove decl-argument-list from stack
                    i += 1
                    continue
            elif self.all_tokens[i].spelling == ";" :
                # ok, there's no function body here, just a function name declaration
                # so we're not going to increase the stack depth
                # in fact, this may not have been a function at all -- it might have
                # been a variale!
                assert( not found_init_list )
                self.set_depth_and_context(i,len(stack),context)
                stack.pop(); # get rid of "function" on stack
                return i+1

            elif self.all_tokens[i].spelling == "{" :
                if ( found_init_list ) :
                    stack.pop()
                # ok! now descend into the function block
                i = self.process_statement( i, stack )
                stack.pop()
                return i

            self.set_depth_and_context(i,len(stack),stack[-1][0])
            i+=1

    def process_for(self,i,stack) :
        print "process_for", i, self.all_tokens[i].spelling
        arglist = "not-yet-begun"
        stack.append(( "for",self.all_tokens[i] ))
        while i < len(self.all_tokens) :
            if not self.all_tokens[i].is_visible :
                pass
            elif self.all_tokens[i].spelling == "(" :
                if arglist == "not-yet-begun" :
                    arglist = "started"
                    stack.append(( "for-declaration",self.all_tokens[i] ))
                    i = self.read_to_end_paren( i, stack )
                    # now remove for-declaration from the stack and process the next statment
                    stack.pop()
                    i = self.process_statement( i+1, stack )
                    stack.pop()
                    return i
            self.set_depth_and_context(i,len(stack),stack[-1][0])
            i += 1
        # whoa! we shouldn't have reached here!
        self.handle_read_all_tokens_error( "process_for", stack )

    def process_if(self,i,stack) :
        print "process_if", i, self.all_tokens[i].spelling
        stack.append(("if",self.all_tokens[i]))
        while i < len(self.all_tokens) :
            if not self.all_tokens[i].is_visible :
                pass
            elif self.all_tokens[i].spelling == "(" :
                stack.append(( "if-condition",self.all_tokens[i] ))
                i = self.read_to_end_paren(i,stack)
                stack.pop() #remove if-condition
                i = self.process_statement( i+1, stack )
                stack.pop() #remove if
                j = self.find_next_visible_token( i )
                if self.all_tokens[j].spelling == "else" :
                    stack.append(( "else",self.all_tokens[i] ))
                    i = self.find_next_visible_token(i,len(stack),"else" )
                    self.set_depth_and_context( i, len(stack), "else" )
                    i = self.process_statement(i+1,stack)
                    stack.pop() #remove else
                return i
            self.set_depth_and_context(i,len(stack),stack[-1][0])
            i += 1

        # whoa! We did not expect to get here!
        self.handle_read_all_tokens_error( "process_if", stack )

    def process_do_while( self, i, stack ) :
        print "process_do_while", i, self.all_tokens[i].spelling
        stack.append(( "do-while",self.all_tokens[i] ))
        self.set_depth_and_context(i,len(stack),stack[-1][0])
        i+=1
        i = self.process_statement(i,stack)
        while i < len(self.all_tokens ) :
            if not self.all_tokens[i].is_visible :
                pass
            elif self.all_tokens[i].spelling == "(" :
                stack.append(( "do-while-condition",self.all_tokens[i] ))
                i = self.read_to_end_paren( i, stack )
                stack.pop()
                return i
            self.set_depth_and_context(i,len(stack),stack[-1][0])
            i+=1

        self.handle_read_all_tokens_error( "process_do_while", stack )

    def process_while(self,i,stack) :
        print "process_while", i, self.all_tokens[i].spelling
        stack.append(( "while",self.all_tokens[i] ))
        while i < len(self.all_tokens) :
            if not self.all_tokens[i].is_visible :
                pass
            elif self.all_tokens[i].spelling == "(" :
                stack.append(( "while-condition",self.all_tokens[i] ))
                i = self.read_to_end_paren( i, stack )
                # ok -- now process the body of the while loop
                self.set_depth_and_context( i, len(stack),stack[-1][0] )
                stack.pop() # remove while-condition
                i = self.process_statement( i, stack )
                stack.pop() # remove while
                return i
            self.set_depth_and_context( i, len(stack), stack[-1][0] )
            i+=1
        # we should not have arrived here
        self.handle_read_all_tokens_error( "process_while", stack )

    def process_namespace(self,i,stack) :
        print "process_namespace", i, self.all_tokens[i].spelling
        stack.append(( "namespace",self.all_tokens[i] ))
        while i < len(self.all_tokens) :
            i = self.find_next_visible_token( i, len(stack), stack[-1][0] )
            if self.all_tokens[i].spelling == "{" :
                i = self.process_statement( i, stack )
                stack.pop()
                return i
            self.set_depth_and_context(i,len(stack),stack[-1][0])
            i+=1
        # we should not have arrived here
        self.handle_read_all_tokens_error( "process_namespace", stack )

    def process_class_decl(self,i,stack) :
        print "process_class_decl", i, self.all_tokens[i].spelling
        seen_inheritance_colon = False
        stack.append(( "class",self.all_tokens[i] ))
        while i < len( self.all_tokens ) :
            if not self.all_tokens[i].is_visible :
                pass
            elif self.all_tokens[i].spelling == "{" :
                if seen_inheritance_colon :
                    stack.pop() # remove class-inheritance-list
                i = self.process_statement( i, stack )
                stack.pop()
                return i
            elif self.all_tokens[i].spelling == ":" :
                if not seen_inheritance_colon :
                    seen_inheritance_colon = True
                    stack.append(( "class-inheritance-list",self.all_tokens[i] ))
            self.set_depth_and_context( i, len(stack), stack[-1][0] )
            i += 1
        # we should not have arrived here
        self.handle_read_all_tokens_error( "process_class_decl", stack )

    def process_union(self,i,stack) :
        print "STUB process_union"
    def process_switch(self,i,stack) :
        print "process_switch", i, self.all_tokens[i].spelling
        stack.append(( "switch",self.all_tokens[i] ))
        while i < len(self.all_tokens) :
            if not self.all_tokens[i].is_visible :
                pass
            if self.all_tokens[i].spelling == "(" :
                stack.append(( "switch-expression",self.all_tokens[i] ))
                i = self.read_to_end_paren( i, stack )
                stack.pop()
                i = self.process_statement( i, stack )
                stack.pop()
                return i
            self.set_depth_and_context(i,len(stack),stack[-1][0])
            i+=1
        self.handle_read_all_tokens_error( "process_class_decl", stack )
    def process_case(self,i,stack) :
        print "process_case", i, self.all_tokens[i].spelling
        stack.append(( "case",self.all_tokens[i] ))
        while i < len(self.all_tokens) :
            if not self.all_tokens[i].is_visible :
                pass
            elif self.all_tokens[i].spelling == ":"  :
                if self.all_tokens[i+1].spelling != ":" :
                  self.set_depth_and_context(i,len(stack),stack[-1][0])
                  while i < len( self.all_tokens ) :
                      i = self.find_next_visible_token(i,len(stack),stack[-1][0])
                      i_spelling = self.all_tokens[i].spelling
                      if i_spelling == "case" or i_spelling == "default" or i_spelling == "}" :
                          stack.pop()
                          return i
                      else :
                          i = self.process_statement(i,stack)
                else :
                    self.set_depth_and_context( i, len(stack), stack[-1][0] )
                    i+=1
            self.set_depth_and_context(i,len(stack),stack[-1][0])
            i+=1 
        # we should not have gotten here
        self.handle_read_all_tokens_error( "process_case_or_default_case", stack )

    def process_scope(self,i,stack) :
        print "process_scope", i, self.all_tokens[i].spelling
        self.set_depth_and_context(i,len(stack),stack[-1][0])
        i+=1
        stack.append(( stack[-1][0]+"-scope",self.all_tokens[i] ))
        if stack[-1][0] not in contexts :
            stack.pop()
            stack.append(( "scope", self.all_tokens[i] ))
        while i < len(self.all_tokens) :
            i = self.find_next_visible_token(i,len(stack),stack[-1][0])
            if self.all_tokens[i].spelling == "}" :
                #ok, done processing at this depth
                stack.pop()
                self.set_depth_and_context(i,len(stack),stack[-1][0])                
                return i+1
            else :
                i = self.process_statement(i,stack)
        self.handle_read_all_tokens_error( "process_scope", stack )

    # def process_scope_end(self,i,stack) :
    #     # print "process_scope_end", i, depth
    #     self.set_depth_and_context(i,depth-1,stack.pop())
    #     return i+1, depth-1
    def process_privacy_declaration(self,i,stack) :
        print "process_privacy_declaration", i, self.all_tokens[i].spelling
        stack.append(( "class-privacy",self.all_tokens[i] ))
        while i < len(self.all_tokens) :
            if not self.all_tokens[i].is_visible :
                pass
            elif self.all_tokens[i].spelling == ":" :
                self.set_depth_and_context(i,len(stack),stack[-1][0])
                stack.pop()
                return i+1
            self.set_depth_and_context(i,len(stack),stack[-1][0])
            i+=1
        self.handle_read_all_tokens_error( "process_privacy_declaration", stack )
