from python_cc_reader.beauty import beautifier
from python_cc_reader.external.blargs import blargs
import sys


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
    print("%s (%s %d) %s %d %s" % ( " " * depth(token), token.parent.type, token.parent.line_number, token.type, token.line_number, token.spelling ))

def print_all_tokens( beaut ) :
    for token in beaut.all_tokens :
        print_token( token )

def test_code_reader(test_name, lines_initial, lines_final):
    beaut = beautifier.Beautifier()
    beaut.filename = test_name
    for line in lines_initial :
        # print line,
        beaut.tokenize_line( line )
    # for i,tok_line in enumerate(beaut.line_tokens) :
    #     print i, [x.spelling for x in tok_line ]
    beaut.minimally_parse_file()
    beaut.beautify_code()

    # print the tokens
    # for tok in beaut.all_tokens :
    #    print tok.spelling, tok.line_number, tok.type, tok.parent.type

    good = True

    # make sure line_tokens and all_tokens agree
    token_counter = -1
    for line_number in range(len( beaut.line_tokens )) :
        for tok in beaut.line_tokens[line_number] :
            token_counter += 1
            if beaut.all_tokens[ token_counter ] is not tok :
                good = False
                print("all_tokens and line_tokens discrepancy", self.all_tokens[ token_counter ].spelling, "vs", tok.spelling)

    for i, line in enumerate( beaut.new_lines ) :
        if i >= len(lines_final) :
            print("Generated too many output lines", line, end=' ')
            good = False
        elif line != lines_final[i] :
            print("Generated line",i+1, "of:", line[:-1], "did not match", lines_final[i][:-1])
            good = False
    if len(beaut.new_lines) < len(lines_final) :
        good = False
        for i in range(len(beaut.new_lines),len(lines_final)) :
            print("Failed to generate expected output line", i+1, ":", lines_final[i], end=' ')

    if not good :
        print_all_tokens( beaut )

        print("Input:")
        for line in lines_initial :
            print(line, end=' ')
        print("Expected:")
        for line in lines_final :
            print(line, end=' ')
        print("Generated:")
        for line in beaut.new_lines :
            print(line, end=' ')
        print("for failed test.")
    else :
        # now, make sure that all the tokens in the tree after beautification
        # represent the same tree that'd be created by parsing the beautified
        # code!
        beaut2 = beautifier.Beautifier()
        for line in lines_final :
            beaut2.tokenize_line( line )
        beaut2.minimally_parse_file()
        if len(beaut2.all_tokens) != len(beaut.all_tokens) :
            good = False
            print("Reparsing the final lines produces a different number of tokens!")
            print("Beautified:", len(beaut.all_tokens), "Parsed Expected Output:", len(beaut2.all_tokens))
        else :
            for i, tok in enumerate(beaut.all_tokens) :
                tok2 = beaut2.all_tokens[i]
                if not tok.equivalent( tok2 ) :
                    tok2
                    print("Tree mismatch:")
                    print("  ", tok.spelling, "vs", tok2.spelling)
                    print("  ", tok.line_number, "vs", tok2.line_number)
                    print("  ", tok.type, "vs", tok2.type)
                    print("  ", tok.context(), "vs", tok2.context())
                    good = False


    beaut = beautifier.Beautifier()
    for line in lines_initial :
        beaut.tokenize_line( line )
    beaut.minimally_parse_file()

    beaut2 = beautifier.Beautifier()
    for line in lines_final :
        beaut2.tokenize_line( line )
    beaut2.minimally_parse_file()

    good2, i_beaut, i_beaut2 = beaut.equivalent( beaut2 )
    if not good2 :
        print("Input lines for test were not found equivalent")
        print("They differ at tokens: ")
        beaut_tok = beaut.all_tokens[i_beaut]
        beaut2_tok = beaut2.all_tokens[i_beaut2]
        print(beaut_tok.type, beaut_tok.line_number, beaut_tok.spelling)
        print(beaut2_tok.type, beaut2_tok.line_number, beaut2_tok.spelling)
        print("lines initial:")
        for line in lines_initial :
            print(line, end=' ')
        print("lines final:")
        for line in lines_final :
            print(line, end=' ')

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
    for i in range(i,len(line)) :
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
            elif line[0] == "%" :
                # matlab style comments; we need to allow c++ style comments and c++ style macros
                # so // and # are both out as possible comment delimiters
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
            ok = test_code_reader(test.name, test.lines_initial, test.lines_final)
            if not ok :
                print("Failed test", test.name)
                print()
            else :
                count_pass += 1
        print("Passed", count_pass, "of", len(tests),"tests.")
        if count_pass != len(tests) :
            sys.exit(1)
