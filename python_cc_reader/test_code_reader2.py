import code_reader2
import blargs

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
