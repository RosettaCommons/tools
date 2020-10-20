from python_cc_reader.cpp_parser.code_reader import *

def test_if_preprocessor_parsing():
    cr = CodeReader()
    cr.push_new_file("unit test")
    cr.examine_line("#define TEST123")
    cr.examine_line("#define TEST456")
    cr.examine_line("#if defined(TEST123) && defined(TEST456)")
    cr.examine_line("#define ABCD")
    cr.examine_line("#else")
    cr.examine_line("#define EFGH")
    cr.examine_line("#endif")

    assert "ABCD" in cr.defined_macros
    assert "EFGH" not in cr.defined_macros

