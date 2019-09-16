import re
from ..cpp_parser import code_reader
from .add_headers import write_file


def auto_ns_comment():
    return " // AUTO USING NS"


# add auto-namespaces to a file, and return the list of namespaces that have been added
# if a namespace has been declared at global scope already, then don't add it.
def add_using_namespaces_to_lines(
    filename, filelines, input_namespaces=["std", "ObjexxFCL", "ObjexxFCL::fmt"]
):

    auto_using_start = "//Auto using namespaces"
    auto_using_end = auto_using_start + " end"
    namespaces = list(input_namespaces)  # deep copy

    auto_using_start_line = 0
    auto_using_block = re.compile(auto_using_start)
    auto_using_block_end = re.compile(auto_using_end)
    using_namespace = "using namespace "
    is_using_ns_dec = re.compile(using_namespace)
    is_include = re.compile("#include")

    using_ns_decs = {}
    for ns in input_namespaces:
        using_ns_decs[ns] = re.compile("using namespace " + ns + ";")

    count_line = 0
    last_auto_using_line = 0
    last_include_line = 0

    cr = code_reader.CodeReader()
    cr.push_new_file(filename)

    # look for the beginning of the auto-namespace block, if one already exists
    # look for using namespace declarations at global scope
    for line in filelines:
        cr.examine_line(line)
        count_line += 1
        if cr.line_is_visible() and is_include.match(line):
            last_include_line = count_line
        if auto_using_start_line > 1 and last_auto_using_line == 0:
            if is_using_ns_dec.search(cr.commentless_line):
                ns = cr.commentless_line.split(using_namespace)[1].split(";")[0]
                if ns not in namespaces:
                    namespaces.append(ns)

        if auto_using_block.match(line):
            if auto_using_start_line != 0:
                if auto_using_block_end.match(line):
                    last_auto_using_line = count_line
            else:
                auto_using_start_line = count_line

        # look for globally-scoped using namespace declarations.
        # last_auto_using_line is >= to auto_using_start_line iff we're outside the auto-using block
        if (
            cr.line_is_visible()
            and cr.scope_level == 0
            and last_auto_using_line >= auto_using_start_line
        ):
            if is_using_ns_dec.search(cr.commentless_line):
                for ns in list(using_ns_decs.keys()):
                    if using_ns_decs[ns].search(cr.commentless_line):
                        # this using namespace already exists at global scope; don't bother adding it
                        if ns in namespaces:
                            namespaces.remove(ns)

    if last_include_line == 0:
        last_include_line = 1  # just put the using declarations at the top!

    if auto_using_start_line != 0 and last_auto_using_line == 0:
        print("Error: While reading", filename)
        print(
            "Error: Found auto-using-namespace start line, but failed to find end line"
        )
        exit

    count_line = 0
    newlines = []

    last_line = ""
    for line in filelines:
        count_line += 1
        if (
            auto_using_start_line != 0
            and count_line >= auto_using_start_line
            and count_line <= last_auto_using_line
        ):
            continue
        write_autoheaders = False
        if auto_using_start_line == 0:
            if count_line == last_include_line + 1:
                write_autoheaders = True
        else:
            if count_line == last_auto_using_line + 1:
                write_autoheaders = True
        if write_autoheaders:
            if last_line != "\n":
                newlines.append("\n")  # pad with one blank line
            newlines.append(auto_using_start + "\n")
            for ns in namespaces:
                ns_decl = ""
                for subns in ns.split("::"):
                    ns_decl += "namespace " + subns + " { "
                for subns in ns.split("::"):
                    ns_decl += "} "
                ns_decl += "using namespace " + ns + ";" + auto_ns_comment() + "\n"
                newlines.append(ns_decl)
            newlines.append(auto_using_end + "\n")
        newlines.append(line)
        last_line = line

    return namespaces, newlines


def add_using_namespaces(
    filename, input_namespaces=["std", "ObjexxFCL", "ObjexxFCL::fmt"]
):
    filelines = open(filename, "r").readlines()
    namespaces, newlines = add_using_namespaces_to_lines(
        filename, filelines, input_namespaces
    )
    write_file(filename, newlines)
    return namespaces


def remove_using_namespace_from_lines(filename, filelines, namespace):
    auto_using_start = "//Auto using namespaces"
    auto_using_end = auto_using_start + " end"

    auto_using_block = re.compile(auto_using_start)
    auto_using_block_end = re.compile(auto_using_end)
    line_sought = re.compile("using namespace " + namespace + ";" + auto_ns_comment())
    inside_auto_using_block = False
    found_sought_namespace = False

    newlines = []

    for line in filelines:
        if not found_sought_namespace:
            if inside_auto_using_block:
                if line_sought.search(line):
                    found_sought_namespace = True
                    continue
                if auto_using_block_end.match(line):
                    print(
                        "Error: While reading",
                        filename,
                        "failed to find",
                        namespace,
                        "in auto-namespace block",
                    )
                    exit
            if auto_using_block.match(line):
                inside_auto_using_block = True
        newlines.append(line)

    if not found_sought_namespace:
        print(
            "Error: While reading",
            filename,
            "failed to find sought namespace or auto-header block!",
        )
        exit

    return newlines


def remove_using_namespace(filename, namespace):
    filelines = open(filename, "r").readlines()
    newlines = remove_using_namespace_from_lines(filename, filelines, namespace)
    write_file(filename, newlines)


def cleanup_auto_namespace_block_lines(filelines):
    auto_using_start = "//Auto using namespaces"
    auto_using_end = auto_using_start + " end"

    auto_using_block = re.compile(auto_using_start)
    auto_using_block_end = re.compile(auto_using_end)

    inside_auto_using_block = False
    newlines = []
    for line in filelines:
        if inside_auto_using_block:
            if not auto_using_block_end.match(line):
                return filelines  # non-emtpy block, do not overwrite
            else:
                inside_auto_using_block = False
                continue
        if auto_using_block.match(line):
            inside_auto_using_block = True
            if newlines[-1] == "\n":
                newlines.pop()  # delete the spacer line, if there is one.
            # print "now inside auto using block"
            continue
        newlines.append(line)
    return newlines


def cleanup_auto_namespace_block(filename):
    filelines = open(filename, "r").readlines()
    newlines = cleanup_auto_namespace_block_lines(filelines)
    write_file(filename, newlines)
