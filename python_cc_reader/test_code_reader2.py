import code_reader2
import blargs

all_funcs = []

def simple_function1() :
    return [ "void add_one( int & input )\n",
             "{\n",
             "\tinput += 1;\n",
             "}\n" ]
all_funcs.append( simple_function1 )

def two_functions() :
    return [ "void add_one( int & input )\n",
             "{\n",
             "\tinput += 1;\n",
             "}\n",
             "\n",
             "void add_two( int & input, std::vector< int > const & some_vector ) {\n",
             "\tinput += 2;",
             "\tsome_vector[0] += 2;\n",
             "}\n" ];
             
all_funcs.append( two_functions )


def function_w_for() :
    return [ "void to_ten()\n",
             "{\n",
             "\tfor ( int i = 0; i < 10; ++i ) {\n",
             "\t\tstd::cout << i << std::endl;\n",
             "\t}\n",
             "}\n" ];
all_funcs.append( function_w_for )

def function_w_while() :
    return [ "void to_ten()\n",
             "{\n",
             "\tint i = 0;",
             "\twhile( i < 10 ) {\n",
             "\t\tstd::cout << i << std::endl;\n",
             "\t\t++i;\n",
             "\t}\n",
             "}\n" ];
all_funcs.append( function_w_while )

def namespace_and_function() :
    return [ "namespace testing {\n",
             "void add_one( int & input )\n",
             "{\n",
             "\tinput += 1;\n",
             "}\n",
             "}\n" ]
all_funcs.append( namespace_and_function )

def class_dec_w_inheritance() :
    return [ "class MyClass : public MyOtherClass {\n",
             "\tMyClass();\n",
             "\t~MyClass();\n",
             "\tvoid testing() { ++i; };\n",
             "};" ]
all_funcs.append( class_dec_w_inheritance )

def test_code_reader_on_lines( lines ) :
    cr = code_reader2.AdvancedCodeReader()
    for line in lines :
        print line,
        cr.tokenize_line( line )
        
    cr.minimally_parse_file()
    for i,line_toks in enumerate( cr.line_tokens ):
        print lines[i],
        for tok in line_toks :
            print "(sp=%s,v=%s,cxt=%s,d=%d)" % ( tok.spelling, tok.is_visible, tok.context, tok.depth ),
        print
    print len(cr.all_tokens)

if __name__ == "__main__" :
    with blargs.Parser(locals()) as p :
        p.str( "filename" )

    if filename :
        lines = open( filename ).readlines()
        test_code_reader_on_lines( lines )
    else :
        for func in all_funcs :
            lines = func()
            test_code_reader_on_lines( lines )
