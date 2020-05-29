
import re
from .add_headers import write_file
from ..cpp_parser import code_reader

# if preserve_regular_header = True, then
# only remove auto-headers from the file
def remove_header_from_filelines(
    fname, filelines, header, preserve_regular_header=False
):
    is_header_angle = re.compile("#include <" + header + ">")
    is_header_quote = re.compile('#include "' + header + '"')
    auto_header_block = re.compile("//Auto Headers")
    ifdef_line = re.compile("#ifdef")
    win32_line = re.compile("#ifdef WIN32")
    endif_line = re.compile("#endif")

    found_auto_block = False
    removed_sought_header = False
    n_nested_win32_ifdefs = 0

    cr = code_reader.CodeReader()
    cr.push_new_file(fname)

    newlines = []
    for line in filelines:
        cr.examine_line(line)

        if removed_sought_header:
            newlines.append(line)
            continue

        if found_auto_block and cr.line_is_visible():
            if not is_header_angle.match(line) and not is_header_quote.match(line):
                newlines.append(line)
            else:
                removed_sought_header = True
        else:
            if auto_header_block.match(line):
                newlines.append(line)
                found_auto_block = True
            elif not cr.line_is_visible():
                newlines.append(line)
            elif is_header_angle.match(line) or is_header_quote.match(line):
                if not preserve_regular_header:
                    line = "// AUTO-REMOVED " + line
                # else, we have a regular header and we've been asked to preserve it
                # so we should leave this line intact
                removed_sought_header = True
                newlines.append(line)
            else:
                newlines.append(line)
    return newlines


def remove_header_from_file(filename, header):
    filelines = open(filename, "r").readlines()
    newlines = remove_header_from_filelines(filename, filelines, header)
    write_file(filename, newlines)


def remove_autoheader_but_not_regular_header(fname, filelines, header):
    return remove_header_from_filelines(fname, filelines, header, True)
