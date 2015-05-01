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

def test_code_reader( lines_initial, lines_final ) :
    cr = code_reader2.AdvancedCodeReader()
    for line in lines_initial :
        cr.tokenize_line( line )
    cr.minimally_parse_file()
    cr.beautify_code()
    good = True
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
                    print "  ", tok.line, "vs", tok2.line
                    print "  ", tok.type, "vs", tok2.type
                    print "  ", tok.context(), "vs", tok2.context()
                    good = False

    return good

def replace_leading_spaces_w_tabs( line ) :
    print "line", line
    newline = []
    found_non_space = False
    if len(line) < 4 :
        "print shorter than 4"
        for i in xrange(len(line)) :
            if line[i] != " " and line[i] != "\n" and line[i] != "\t" :
                found_non_space = True
    else :
        i = 0
        while i+4 < len(line) :
            print i, line[i:(i+4)]
            if line[i:(i+4)] == "    " :
                newline.append( "\t" )
                i+=4
            else :
                found_non_space = True
                break

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
                print init_lines
                print final_lines
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
        for test in tests :
            ok = test_code_reader( test.lines_initial, test.lines_final )
            if not ok :
                print "Failed test", test.name
                print
            else :
                count_pass += 1
        print "Passed", count_pass, "of", len(tests),"tests."
        if count_pass != len(tests) :
            sys.exit(1)
