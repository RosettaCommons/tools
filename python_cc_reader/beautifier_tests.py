import code_reader2
import blargs
import sys

all_funcs = []

def simple_function1() :
    return [ "void add_one( int & input )\n",
             "{\n",
             "input += 1;\n",
             "}\n" ]
all_funcs.append( simple_function1 )

def two_functions() :
    return [ "void add_one( int & input )\n",
             "{\n",
             "input += 1;\n",
             "}\n",
             "\n",
             "void add_two( int & input, std::vector< int > const & some_vector ) {\n",
             "input += 2;",
             "some_vector[0] += 2;\n",
             "}\n" ];
             
all_funcs.append( two_functions )


def function_w_for() :
    return [ "void to_ten()\n",
             "{\n",
             "for(int i = 0; i < 10; ++i){\n",
             "std::cout << i << std::endl;\n",
             "}\n",
             "}\n" ];
all_funcs.append( function_w_for )

def function_w_while() :
    return [ "void to_ten()\n",
             "{\n",
             "int i = 0;",
             "while( i < 10 ) {\n",
             "std::cout << i << std::endl;\n",
             "++i;\n",
             "}\n",
             "}\n" ];
all_funcs.append( function_w_while )

def namespace_and_function() :
    return [ "namespace testing {\n",
             "void add_one( int & input )\n",
             "{\n",
             "input += 1;\n",
             "}\n",
             "}\n" ]
all_funcs.append( namespace_and_function )

def class_dec_w_inheritance() :
    return [ "class MyClass : public MyOtherClass {\n",
             "MyClass();\n",
             "~MyClass();\n",
             "void testing() { ++i; };\n",
             "};" ]
all_funcs.append( class_dec_w_inheritance )

def class_dec_w_ctor() :
    return [ "class MyClass {\n",
             "MyClass() : my_int_( 5 ) {}\n",
             "int my_int_;\n",
             "};\n" ]
all_funcs.append( class_dec_w_ctor )

def class_dec_w_privacy() :
    return [ "class MyClass {\n",
             "private :\n",
             "MyClass();\n",
             "MyClass( int );\n",
             "};\n" ]
all_funcs.append( class_dec_w_privacy )

def test_code_reader_on_lines( lines ) :
    cr = code_reader2.AdvancedCodeReader()
    for line in lines :
        print line,
        cr.tokenize_line( line )
    print
    
    cr.minimally_parse_file()
    for tok in cr.all_tokens :
        assert( tok.parent != None )
        assert( tok in tok.parent.children )

    # for i,line_toks in enumerate( cr.line_tokens ):
    #     print lines[i],
    #     for tok in line_toks :
    #         print "(sp=%s,v=%s,cxt=%s,d=%d)" % ( tok.spelling, tok.is_visible, tok.context, tok.depth ),
    #     print
    # print len(cr.all_tokens)

    cr.beautify_code()
    print "Beautified"
    for line in cr.new_lines :
        print line,
    print

class BeautifierTest :
    def __init__( self, li, lf, name ) :
        self.lines_initial = li
        self.lines_final = lf
        self.name = name

def depth( token ) :
    if token.parent :
        return 1 + depth( token.parent )
    else :
        return 0

def print_token( token ) :
    print "%s (%s %d) %s %d %s" % ( " " * depth(token), token.parent.type, token.parent.line_number, token.type, token.line_number, token.spelling )

def print_all_tokens( cr ) :
    for token in cr.all_tokens :
        print_token( token )

def test_code_reader( lines_initial, lines_final ) :
    cr = code_reader2.AdvancedCodeReader()
    for line in lines_initial :
        # print line,
        cr.tokenize_line( line )
    # for i,tok_line in enumerate(cr.line_tokens) :
    #     print i, [x.spelling for x in tok_line ]
    cr.minimally_parse_file()
    cr.beautify_code()

    # print the tokens
    # for tok in cr.all_tokens :
    #    print tok.spelling, tok.line_number, tok.type, tok.parent.type

    good = True

    # make sure line_tokens and all_tokens agree
    token_counter = -1
    for line_number in xrange(len( cr.line_tokens )) :
        for tok in cr.line_tokens[line_number] :
            token_counter += 1
            if cr.all_tokens[ token_counter ] is not tok :
                good = False
                print "all_tokens and line_tokens discrepancy", self.all_tokens[ token_counter ].spelling, "vs", tok.spelling

    for i, line in enumerate( cr.new_lines ) :
        if i >= len(lines_final) :
            print "Generated too many output lines", line,
            good = False
        elif line != lines_final[i] :
            print "Generated line",i+1, "of:", line[:-1], "did not match", lines_final[i][:-1]
            good = False
    if len(cr.new_lines) < len(lines_final) :
        good = False
        for i in xrange(len(cr.new_lines),len(lines_final)) :
            print "Failed to generate expected output line", i+1, ":", lines_final[i],

    if not good :
        print_all_tokens( cr )

        print "Input:"
        for line in lines_initial :
            print line,
        print "Expected:"
        for line in lines_final :
            print line,
        print "Generated:"
        for line in cr.new_lines :
            print line,
        print "for failed test."
    else :
        # now, make sure that all the tokens in the tree after beautification
        # represent the same tree that'd be created by parsing the beautified
        # code!
        cr2 = code_reader2.AdvancedCodeReader()
        for line in lines_final :
            cr2.tokenize_line( line )
        cr2.minimally_parse_file()
        if len(cr2.all_tokens) != len(cr.all_tokens) :
            good = False
            print "Reparsing the final lines produces a different number of tokens!"
            print "Beautified:", len(cr.all_tokens), "Parsed Expected Output:", len(cr2.all_tokens)
        else :
            for i, tok in enumerate(cr.all_tokens) :
                tok2 = cr2.all_tokens[i]
                if not tok.equivalent( tok2 ) :
                    tok2
                    print "Tree mismatch:"
                    print "  ", tok.spelling, "vs", tok2.spelling
                    print "  ", tok.line_number, "vs", tok2.line_number
                    print "  ", tok.type, "vs", tok2.type
                    print "  ", tok.context(), "vs", tok2.context()
                    good = False


    cr = code_reader2.AdvancedCodeReader()
    for line in lines_initial :
        cr.tokenize_line( line )
    cr.minimally_parse_file()

    cr2 = code_reader2.AdvancedCodeReader()
    for line in lines_final :
        cr2.tokenize_line( line )
    cr2.minimally_parse_file()

    good2, i_cr, i_cr2 = cr.equivalent( cr2 )
    if not good2 :
        print "Input lines for test were not found equivalent"
        print "They differ at tokens: "
        cr_tok = cr.all_tokens[i_cr]
        cr2_tok = cr2.all_tokens[i_cr2]
        print cr_tok.type, cr_tok.line_number, cr_tok.spelling
        print cr2_tok.type, cr2_tok.line_number, cr2_tok.spelling
        print "lines initial:"
        for line in lines_initial :
            print line,
        print "lines final:"
        for line in lines_final :
            print line,

    return good and good2

def replace_leading_spaces_w_tabs( line ) :
    # print "line", line
    newline = []
    found_non_space = False
    i = 0
    while i+4 < len(line) :
        # print i, line[i:(i+4)]
        if line[i:(i+4)] == "    " :
            newline.append( "\t" )
            i+=4
        else :
            found_non_space = True
            break

    # now look at the characters we didn't examine; are any of them non-whitespace?
    for i in xrange(i,len(line)) :
        if line[i] != " " and line[i] != "\n" and line[i] != "\t" :
            found_non_space = True

    if found_non_space :
        newline.append( line.strip() )
        newline.append( "\n" )
        return "".join(newline)
    else :
        return "\n"

def trim_empty_blank_lines( lines ) :
    i = len( lines ) - 1
    while i > 0 :
        if lines[i] != "\n" :
            return lines[:(i+1)]
        i -= 1
    return []
        

def read_test_cases( fname ) :
    lines = open( fname ).readlines()
    tests = []
    init_lines = []
    final_lines = []
    for line in lines :
        if len(line) > 3 :
            if line[0:3] == "---" :
                init_lines = trim_empty_blank_lines( init_lines )
                final_lines = trim_empty_blank_lines( final_lines )
                # print init_lines
                # print final_lines
                test = BeautifierTest( init_lines, final_lines, line[3:].strip() )
                init_lines = []; final_lines = []
                tests.append( test )
            elif line[0] == "#" :
                continue
            else :
                cols = line.split("|")
                assert( len(cols) == 2 )
                init_line = replace_leading_spaces_w_tabs( cols[0] )
                init_lines.append( init_line )
                if len(cols[1]) > 0 :
                    final_line = replace_leading_spaces_w_tabs( cols[1][1:] )
                    final_lines.append( final_line )
    assert( not init_lines ) # test cases should always end with ---
    assert( not final_lines ) # test cases should always end with ---
    return tests
                    

if __name__ == "__main__" :
    with blargs.Parser(locals()) as p :
        p.str( "filename" )
        p.flag( "old_tests" )
        p.int("just_one")

    if filename :
        lines = open( filename ).readlines()
        test_code_reader_on_lines( lines )
    elif old_tests :
        for func in all_funcs :
            lines = func()
            test_code_reader_on_lines( lines )
    else :
        tests = read_test_cases( "beautifier_test_cases.txt" )
        count_pass = 0
        for i,test in enumerate( tests ) :
            if just_one is not None and i+1 != just_one : continue
            ok = test_code_reader( test.lines_initial, test.lines_final )
            if not ok :
                print "Failed test", test.name
                print
            else :
                count_pass += 1
        print "Passed", count_pass, "of", len(tests),"tests."
        if count_pass != len(tests) :
            sys.exit(1)
