
import re
from ..cpp_parser import code_reader
from .add_headers import write_file
from ..cpp_parser.code_utilities import load_source_tree


def remove_duplicate_headers_from_filelines(filename, filelines):
    included_files = {}
    include_line = re.compile("\#include")
    quote_include = re.compile('\#include\s*"')
    angle_include = re.compile("\#include\s*<")

    cr = code_reader.CodeReader()
    cr.push_new_file(filename)
    newlines = []

    for line in filelines:
        cr.examine_line(line)

        if cr.line_is_visible() and include_line.match(line):
            if quote_include.match(line):
                include = line.split('"')[1]
            elif angle_include.match(line):
                include = line.split("<")[1].split(">")[0]
            else:
                print(line)
                sys.exit(1)

            if include not in included_files:
                included_files[include] = 0
            else:
                line = "// Auto-header: duplicate removed " + line
            newlines.append(line)
        else:
            newlines.append(line)

    return newlines


def remove_duplicate_headers_from_file(filename):
    newlines = remove_duplicate_headers_from_filelines(
        filename, open(filename).readlines()
    )
    write_file(filename, newlines)


if __name__ == "__main__":
    compilable_files, all_includes, file_contents = load_source_tree()
    for cc in compilable_files:
        newlines = remove_duplicate_headers_from_filelines(cc, file_contents[cc])
        write_file(cc, newlines)
