import re
from ..cpp_parser import code_reader


def write_file(filename, filelines):
    file = open(filename, "w")
    file.writelines(filelines)
    file.close()


def group_headers(headers):

    # order should be: core, protocols, devel, apps,
    # utility, numeric, ObjexxFCL, STL, everything else
    lib_names = [
        "platform",
        "core",
        "protocols",
        "devel",
        "apps",
        "utility",
        "numeric",
        "ObjexxFCL",
    ]
    other_names = ["stl", "other"]
    header_sets = {}
    for name in lib_names:
        header_sets[name] = {}
    for name in other_names:
        header_sets[name] = []
    for h in headers:
        if len(h.split("/")) > 1:
            libname = h.split("/")[0]
            if libname in lib_names:
                nsname = h[0 : h.rfind("/")]
                if nsname not in header_sets[libname]:
                    header_sets[libname][nsname] = []
                header_sets[libname][nsname].append(h)
            else:
                header_sets["other"].append(h)
        else:
            header_sets["stl"].append(h)
    for libname in lib_names:
        for nsname in header_sets[libname]:
            header_sets[libname][nsname].sort()
    for name in other_names:
        header_sets[name].sort()
    return lib_names, other_names, header_sets


def group_and_sort_headers(headers):
    lib_names, other_names, header_sets = group_headers(headers)
    all_headers = []
    for libname in lib_names:
        namespaces = list(header_sets[libname].keys())
        namespaces.sort()
        for ns in namespaces:
            # print libname, ns
            all_headers.extend(header_sets[libname][ns])
    for name in other_names:
        all_headers.extend(header_sets[name])
    return all_headers


def add_autoheaders_to_filelines(fname, filelines, headers):
    if len(headers) == 0:
        return filelines
    auto_header_start = "//Auto Headers"
    sorted_headers = group_and_sort_headers(headers)

    auto_header_start_line = 0
    auto_header_block = re.compile(auto_header_start)
    is_include = re.compile("#include")
    quote_include = re.compile('#include\s*"')
    angle_include = re.compile("#include\s*<")
    pound_define = re.compile("#define")
    count_line = 0
    last_include_line = 0
    pound_define_line = 0

    cr = code_reader.CodeReader()
    cr.push_new_file(fname)
    for line in filelines:
        cr.examine_line(line)
        count_line += 1
        if cr.line_is_visible() and cr.scope_level == 0 and is_include.match(line):
            last_include_line = count_line
            if auto_header_start_line > 1:
                if quote_include.match(line):
                    include = line.split('"')[1]
                elif angle_include.match(line):
                    include = line.split("<")[1].split(">")[0]
                sorted_headers.append(include)
        if auto_header_block.match(line):
            auto_header_start_line = count_line
            if count_line > last_include_line:
                last_include_line = count_line  # treat this line as the last include
        if (
            cr.line_is_visible()
            and pound_define.match(line)
            and cr.scope_level == 0
            and pound_define_line == 0
        ):
            pound_define_line = count_line

    # print last_include_line

    if auto_header_start_line > 0:
        sorted_headers = group_and_sort_headers(sorted_headers)
        last_header = ""
        temp = []
        for header in sorted_headers:
            if header != last_header:
                temp.append(header)
            last_header = header
        sorted_headers = temp

    if last_include_line == 0:
        if pound_define_line == 0:
            last_include_line = 1  # just put the headers at the top!
        else:
            last_include_line = pound_define_line  # put the new headers after a #define

    count_line = 0
    newlines = []

    for line in filelines:
        count_line += 1
        if (
            auto_header_start_line != 0
            and count_line >= auto_header_start_line
            and count_line <= last_include_line
        ):
            # assert( is_include( line ).match() ) # don't skip anything that's not a header!
            continue
        if count_line == last_include_line + 1:
            if auto_header_start_line == 0:
                newlines.append("\n")
            newlines.append(auto_header_start + "\n")
            for h in sorted_headers:
                newlines.append("#include <" + h + ">\n")
            if auto_header_start_line == 0:
                newlines.append("\n")
        newlines.append(line)

    return newlines


def add_autoheaders_to_file(filename, headers):
    if len(headers) == 0:
        return

    filelines = open(filename, "r").readlines()
    newlines = add_autoheaders_to_filelines(filename, filelines, headers)
    write_file(filename, newlines)


def add_autoheaders_to_file_list(list_of_files, headers):
    file = open(list_of_files, "r")
    file_list = []
    for line in file:
        file_list.append(line.rstrip())

    for path in file_list:
        add_headers.add_autoheaders_to_file(path, headers)

def count_autoheaders(fname, filelines):
    auto_header_start = "//Auto Headers"
    is_include = re.compile("#include")

    cr = code_reader.CodeReader()
    cr.push_new_file(fname)
    include_count = 0
    auto_headers_begun = False
    for line in filelines:
        cr.examine_line(line)
        if not auto_headers_begun:
            if line.startswith(auto_header_start):
                auto_headers_begun = True
        
        elif cr.line_is_visible() and cr.scope_level == 0 and is_include.match(line):
            include_count += 1
    return include_count

def remove_autoheader_scar(fname, filelines):
    if count_autoheaders(fname, filelines) > 0:
        return
    auto_header_start = "//Auto Headers"

    cr = code_reader.CodeReader()
    cr.push_new_file(fname)

    auto_headers_begun = False
    blank_line_after = False
    
    newlines = []
    for line in filelines:
        cr.examine_line(line)
        if not auto_headers_begun:
            if line.startswith(auto_header_start):
                auto_headers_begun = True
            else:
                newlines.append(line)
        else:
            if not blank_line_after:
                if line == "\n":
                    blank_line_after = True
                    continue
            newlines.append(line)
    with open(fname, "w") as fid:
        fid.writelines(newlines)

        
