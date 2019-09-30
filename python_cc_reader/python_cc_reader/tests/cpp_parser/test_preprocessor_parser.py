from python_cc_reader.cpp_parser.preprocessor_parser import *

def test_tokenize_1():
    scanner = PreProcScanner()
    toks = scanner.tokenize_line("defined(MULTITHREADED)")
    assert len(toks) == 4
    assert toks[0].type == PreProcTokenTypes.DEFINED
    assert toks[1].type == PreProcTokenTypes.LPAREN
    assert toks[2].type == PreProcTokenTypes.VAR
    assert toks[2].spelling == "MULTITHREADED"
    assert toks[3].type == PreProcTokenTypes.RPAREN

def test_tokenize_2():
    scanner = PreProcScanner()
    toks = scanner.tokenize_line("defined(WIN32) || defined(__CYGWIN__)")
    assert len(toks) == 9
    assert toks[0].type == PreProcTokenTypes.DEFINED
    assert toks[1].type == PreProcTokenTypes.LPAREN
    assert toks[2].type == PreProcTokenTypes.VAR
    assert toks[3].type == PreProcTokenTypes.RPAREN
    assert toks[4].type == PreProcTokenTypes.OR
    assert toks[5].type == PreProcTokenTypes.DEFINED
    assert toks[6].type == PreProcTokenTypes.LPAREN
    assert toks[7].type == PreProcTokenTypes.VAR
    assert toks[8].type == PreProcTokenTypes.RPAREN
    

def test_preproc_parser_1():
    parser = PreProcParser()
    ast = parser.tree_from_line("defined(WIN32) || defined(__CYGWIN__)")
    assert ast
    assert ast.type == PreProcTokenTypes.OR
    left_def = ast.left
    right_def = ast.right
    assert left_def.type == PreProcTokenTypes.DEFINED
    assert right_def.type == PreProcTokenTypes.DEFINED
    left_var = left_def.sub_expr
    right_var = right_def.sub_expr
    assert left_var.type == PreProcTokenTypes.VAR
    assert right_var.type == PreProcTokenTypes.VAR
    assert left_var.varname == "WIN32"
    assert right_var.varname == "__CYGWIN__"
