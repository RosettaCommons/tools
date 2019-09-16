from ..external.blargs import blargs
from . import beautifier
import sys

def parse_file(lines):
    beaut = beautifier.Beautifier()
    for line in lines:
        beaut.tokenize_line(line)
    beaut.minimally_parse_file()
    return beaut


def compare_whether_two_files_are_equivalent(file1, file2):
    b1 = parse_file(open(file1).readlines())
    b2 = parse_file(open(file2).readlines())
    good, i1, i2 = b1.equivalent(b2)
    if not good:
        print("The files", file1, "and", file2, "are not equivalent.")
        print("They differ at tokens:")
        if i1 < len(b1.all_tokens):
            tok1 = b1.all_tokens[i1]
            print(
                "type:",
                tok1.type,
                "on line",
                tok1.line_number,
                "with spelling",
                tok1.spelling,
            )
        else:
            print("File1 token is out of range")
        if i2 < len(b2.all_tokens):
            tok2 = b2.all_tokens[i2]
            print(
                "type:",
                tok2.type,
                "on line",
                tok2.line_number,
                "with spelling",
                tok2.spelling,
            )
        else:
            print("File2 token is out of range")
        sys.exit(1)
