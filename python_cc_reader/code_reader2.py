import sys


contexts = [ "namespace",
             "solitary-scope",
             "for",
             "for-declaration",
             "do-while",
             "do-while-declaration",
             "do-while-condition",
             "if",
             "if-condition",
             "else",
             "else-if",
             "while",
             "while-condition"
             "switch",
             "switch-variable",
             "case",
             "class",
             "class-inheritance-list",
             "struct",
             "struct-inheritance-list",
             "union",
             "function",
             "function-decl-argument-list",
             "function-call-argument-list",
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
        self.dividers = set([";",":",",","(",")","{","}","=","[","]"])
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
        for tok in self.this_line_tokens : tok.visible = False
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

    def find_next_visible_token( self, i ) :
        while i < len(self.all_tokens) :
            if self.all_tokens[i].is_visible : return i
            i += 1
        return i

    def minimally_parse_file( self ) :
        # at this point, the file has been split into a large set of tokens
        # and we're going to go through the tokens and, to a limited extent,
        # interpret the tokens into a structure
        i = 0
        depth = 1
        stack = ["namespace"] # you start at namespace scope
        while i < len( self.all_tokens ) :
            print "while top", i, depth
            i = self.find_next_visible_token(i)
            if i > len( self.all_tokens ) : break
            i_spelling = self.all_tokens[i].spelling
            if i_spelling == "for" :
                i,depth = self.process_for(i,depth,stack)
            elif i_spelling == "if" :
                i,depth = self.process_if(i,depth,stack)
            elif i_spelling == "else" :
                i,depth = self.process_else(i,depth,stack)
            elif i_spelling == "while" :
                i,depth = self.process_while(i,depth,stack)
            elif i_spelling == "namespace" :
                i,depth = self.process_namespace(i,depth,stack)
            elif i_spelling == "class" or i_spelling == "struct" :
                i,depth = self.process_class_decl(i,depth,stack)
            elif i_spelling == "union" :
                i,depth = self.process_union_decl(i,depth,stack)
            elif i_spelling == "switch" :
                i,depth = self.process_switch(i,depth,stack)
            elif i_spelling == "case" :
                i,depth = self.process_case(i,depth,stack)
            elif i_spelling == "{" :
                i,depth = self.process_solitary_scope(i,depth,stack)
            elif i_spelling == "}" :
                i,depth = self.process_scope_end(i,depth,stack)
            elif i_spelling == "using" or i_spelling == "typedef" :
                i,depth = self.process_statement(i,depth,stack)
            elif i_spelling == "public" or i_spelling == "protected" or i_spelling == "private" :
                i,depth = self.process_privacy_declaration(i,depth,stack)
            elif stack[-1] == "namespace" or stack[-1] == "class" or stack[-1] == "struct" or stack[-1] == "union" :
                i,depth = self.process_function_preamble_or_variable(i,depth,stack)
            else :
                i,depth = self.process_statement(i,depth,stack)

    def print_depth_stack( self, stack ) :
        print "Depth stack: "
        for elem in reverse( stack ) :
            print "  ", elem


    def process_statement(self,i,depth,stack) :
        # print "process statement", i, depth
        # read all the tokes from i until the next semicolon
        last = i
        while last < len(self.all_tokens) :
            if self.all_tokens[last].is_visible and self.all_tokens[last].spelling == ";" :
                # print "statement:", " ".join( [ x.spelling for x in self.all_tokens[i:last+1] ] )
                break
            last += 1
        if last == len(self.all_tokens) :
            print "Unexpectedly ran out of tokens!"
            self.print_depth_stack( stack )
            sys.exit(1)
        for j in xrange(i,last+1) :
            self.all_tokens[j].context = stack[-1]
            self.all_tokens[j].depth = depth
        return last+1,depth

    def process_function_preamble_or_variable(self,i,depth,stack) :
        # print "process_function_preamble_or_variable", i, depth
        # treat function declarations and variable declarations interchangably.
        # so we should gobble up the tokens that make this a function
        # and see if we have a constructor, in which case, we should also
        # be watchful for an intializer list
        classname = ""
        arglist = "not-begun"
        is_ctor = False
        context = stack[-1]
        paren_depth = 0
        found_init_list = False
        while i < len(self.all_tokens) :
            if not self.all_tokens[i].is_visible :
                # don't do anything special; just don't try and interpret this token as having meaning
                pass
            elif self.all_tokens[i].spelling == ":" :
                if arglist == "ended" :
                    # we're now reading the initializer list
                    context = "ctor-initializer-list"
                    depth += 1
                    found_init_list = True
                elif arglist != "started" :
                    assert( arglist == "not-begun" )
                    # this might be constructor, we don't know
                    # this might also be some namespace qualifier on a return type
                    if i+1 < len(self.all_tokens) and self.all_tokens[i+1].spelling == ":" :
                        classname = self.all_tokens[i-1].spelling
            elif self.all_tokens[i].spelling == "(" :
                if arglist == "started" :
                   paren_depth += 1
                else :
                    funcname = self.all_tokens[i-1].spelling
                    if funcname != "operator" :
                        if funcname == classname :
                            is_ctor = True
                        arglist = "started"
                        context = "function-decl-argument-list"
            elif self.all_tokens[i].spelling == ")" :
                if arglist == "started" :
                    paren_depth -= 1
                    if paren_depth == 0 :
                        arglist = "ended"
                        self.all_tokens[i].context = context
                        self.all_tokens[i].depth = depth
                        context = "function"
                        i += 1
                        continue
            elif self.all_tokens[i].spelling == ";" :
                # ok, there's no function body here, just a function name declaration
                # so we're not going to increase the stack depth
                # in fact, this may not have been a function at all -- it might have
                # been a variale!
                assert( not found_init_list )
                self.all_tokens[i].context = context
                self.all_tokens[i].depth = depth
                return i+1, depth
            elif self.all_tokens[i].spelling == "{" :
                if ( found_init_list ) :
                    depth -= 1
                # ok! we're about to add a function onto the stack
                self.all_tokens[i].context = "function"
                self.all_tokens[i].depth = depth
                stack.append( "function" )
                return i+1, depth+1

            self.all_tokens[i].context = context
            self.all_tokens[i].depth = depth
            i+=1            
        # whoa! this function wasn't supposed to reach here
        print "ERROR: reach end of tokens while in process_function_preamble_or_variable"
        self.print_depth_stack( stack )
    def process_for(self,i,depth,stack) :
        print "STUB process_for"
    def process_if(self,i,depth,stack) :
        print "STUB process_if"
    def process_else(self,i,depth,stack) :
        print "STUB process_else"
    def process_while(self,i,depth,stack) :
        print "STUB process_while"
    def process_namespace(self,i,depth,stack) :
        print "STUB process_namespace"
    def process_class_decl(self,i,depth,stack) :
        print "STUB process_class_decl"
    def process_union(self,i,depth,stack) :
        print "STUB process_union"
    def process_switch(self,i,depth,stack) :
        print "STUB process_switch"
    def process_case(self,i,depth,stack) :
        print "STUB process_case"
    def process_solitary_scope(self,i,depth,stack) :
        print "STUB process_solitary_scope"
    def process_scope_end(self,i,depth,stack) :
        # print "process_scope_end", i, depth
        self.all_tokens[i].context = stack.pop()
        return i+1, depth-1
    def process_privacy_declaration(self,i,depth,stack) :
        print "STUB process_privacy_declaration"

    
