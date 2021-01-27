
# from test_compile import find_includes_for_file

from . import code_reader
import re, sys, os


def regex_subset(filelist, regex):
    subset = []
    for file in filelist:
        if regex.match(file):
            subset.append(file)
    return subset


def regex_exclude(filelist, regex):
    subset = []
    for file in filelist:
        if not regex.match(file):
            subset.append(file)
    return subset


def cc_subset(filelist):
    is_ccfile = re.compile("\S*\.cc$")
    return regex_subset(filelist, is_ccfile)


def fwdhh_regex():
    return re.compile("\S*\.fwd\.hh$")


def fwd_hh_subset(filelist):
    is_fwdhh = fwdhh_regex()
    return regex_subset(filelist, is_fwdhh)


def non_fwd_hh_subset(filelist):
    is_fwdhh = fwdhh_regex()
    all_hh = re.compile("\S*\.hh$")
    allhhs = regex_subset(filelist, all_hh)
    return regex_exclude(allhhs, is_fwdhh)


def libraries_with_ccfiles_to_examine():
    """The list of libraries in rosetta which contain .cc files that are built by scons.

   This list may also be used to determine the subset of .hh files whose headers
   likely contain extra #includes that need to be removed.  Headers contained
   in other libraries should not be modified."""
    return [
        "basic",
        "core.1",
        "core.2",
        "core.3",
        "core.4",
        "core.5",
        "protocols",
        "devel",
        "apps",
        "pilot_apps",
    ]


def libraries_with_hhfiles_to_examine():
    return libraries_with_ccfiles_to_examine() + ["utility", "numeric", "ObjexxFCL"]


def directories_with_ccfiles_to_examine():
    """The list of directories in rosetta which are turned into libraries that are built by scons
   This list is used to determine if a file should be examined.  Note that
   there is no longer a perfect correspondence between the libraries that are
   produced and the directories that contain them"""
    return ["basic", "core", "protocols", "devel", "apps"]


def directories_with_hhfiles_to_examine():
    return directories_with_ccfiles_to_examine() + ["utility", "numeric", "ObjexxFCL"]


# read a line which begins with #include and return the file name that
# is #included.  This will read both <> includes or "" includes.
def include_for_line(cr, line):
    re_splitter = re.compile("\S+")
    toks = re_splitter.findall(line)
    if len(toks) < 2:
        if len(toks) > 0:
            if toks[0].find("<") != -1:
                return toks[0].split("<")[1].split(">")[0]
        print("Error", cr.curr_file(), cr.curr_line(), line)
        print("#include is missing a filename")
        sys.exit(1)
    if toks[1].find("<") != -1:
        include = toks[1].split("<")[1].split(">")[0]
    else:
        include = toks[1].split('"')[1]
    return include


def find_all_includes(name, filelines):
    inclusions = []
    cr = code_reader.CodeReader()
    cr.push_new_file(name)
    re_pound_include = re.compile("#include")
    for line in filelines:
        cr.examine_line(line)
        if cr.line_is_visible():
            if re_pound_include.match(line):
                include = include_for_line(cr, line)
                if include not in inclusions:
                    inclusions.append(include)
    return inclusions


def find_includes_at_global_scope(name, filelines):
    inclusions = []
    cr = code_reader.CodeReader()
    cr.push_new_file(name)
    re_pound_include = re.compile("#include")
    for line in filelines:
        cr.examine_line(line)
        if cr.line_is_visible() and cr.scope_level == 0:
            if re_pound_include.match(cr.commentless_line):
                # print(cr.line_is_visible(), cr.curr_file(), cr.curr_line(), line, cr.commentless_line,)
                include = include_for_line(cr, cr.commentless_line)
                if include not in inclusions:
                    inclusions.append(include)
    return inclusions


# return a list of all the non-utility cc files listed for
# compilation in the *.src.settings files.
# NOTE: I should be reading pilot_apps.src.settings.all!
def compiled_cc_files():
    to_be_compiled = []
    projects = rosetta_projects()

    for library in projects["src"]:
        settings_file_name = library + ".src.settings"
        libname = library
        if library == "pilot_apps":
            libname = "apps"
            settings_file_name += ".all"
        f = open(settings_file_name).read()
        # print("reading files for library", libname)
        _globals = {}
        _locals = {}
        exec(f, _globals, _locals)
        sources = _locals["sources"]
        for dir in list(sources.keys()):
            for cc_file in sources[dir]:
                if len(cc_file) > 3 and cc_file[-3:] == ".cu":
                    # skip cuda files
                    continue
                if dir == "":
                    name = (
                        cc_file + ".cc"
                    )  # note, following Sam's change to the .src.settings files, no longer append libname to file name
                elif libname == "apps":
                    name = libname + "/" + dir + "/" + cc_file + ".cc"
                else:
                    name = (
                        dir + "/" + cc_file + ".cc"
                    )  # note, following Sam's change to the .src.settings files, no longer append libname to file name
                # print("file to be compiled: ", name)
                to_be_compiled.append(name)
    return to_be_compiled


def compiled_cxxtest_hh_files():
    to_be_compiled = []
    # the unit tests were not split when core was split
    # so instead of dealing with core.1, core.2 etc, just look at
    # the
    for library in directories_with_ccfiles_to_examine():
        settings_file_name = "../test/" + library + ".test.settings"
        libname = library
        if library == "pilot_apps":
            libname = "apps"
        if not os.path.isfile(settings_file_name):
            print("did not find", settings_file_name, ". Skipping")
            continue
        f = open(settings_file_name).read()
        _globals = {}
        _locals = {}
        exec(f, _globals, _locals)
        sources = _locals["sources"]
        for dir in list(sources.keys()):
            for hh_file in sources[dir]:
                name = "../test/" + libname + "/" + dir + "/" + hh_file + ".cxxtest.hh"
                to_be_compiled.append(name)
    return to_be_compiled


# return all the libraries that are to be built
def rosetta_projects():
    f = open("../projects.settings").read()
    _globals = {}
    _locals = {}
    exec(f, _globals, _locals)
    return _locals["projects"]


# return a set of all the top-level directories in a particular .src.settings file
# intended for the sub-divided libraries and for use in the library_levels.py script
def toplevel_subdirs_of_library(libname):
    settings_file_name = libname + ".src.settings"
    f = open(settings_file_name).read()
    _globals = {}
    _locals = {}
    exec(f, _globals, _locals)
    sources = _locals["sources"]
    subdirs = set([])
    for dir in list(sources.keys()):
        subdir_full = dir.partition("/")[2]
        subdir_toplevel = subdir_full.partition("/")[0]
        if subdir_toplevel not in subdirs:
            subdirs.add(subdir_toplevel)
    return subdirs


# return a dictionary containing the #includes of all
# internal source files.  By "internal", I mean non-stl,
# non-external (e.g. boost) source files.  This dictionary
# includes both .cc files and .hh files.  It only examines
# those files which are reached during compilation:
# those .cc files listed in a .src.settings file, and any
# header file that is #included at global scope
def scan_compilable_files():
    all_includes = {}
    unprocessed_filenames = compiled_cc_files()
    for name in unprocessed_filenames:
        toks = name.split("/")
        if len(toks) < 2 or toks[0] not in directories_with_hhfiles_to_examine():
            # print("skipping past", name)
            continue
        # print("reading includes from", name)
        includes = find_includes_at_global_scope(name, open(name).readlines())
        for include in includes:
            if include not in unprocessed_filenames:
                unprocessed_filenames.append(include)
        all_includes[name] = includes
    return all_includes


# opens the scons *.src.settings files, examines the #includes for all the
# .cc files that the settings files instruct scons to build, and generates
# a list of all the header files that are #included at global scope.
# All these header files and .cc files should be independently compilable.
# It then loads the contents of all files #included at any scope.
# This function returns a pair: the list of all files that should
# be independently compilable, and the contents of all files.
#
# ASSUMPTION: if a header file is "injected" into another file (i.e. it's
# #included at non-global scope, e.g. utility/keys/KeyLookupFunctors.hh),
# then that header does not #include any other header.
# If this assumption is violated, then the contents of that second-order
# injuected header will not be placed into the source dictionary
# that this function returns
def load_source_tree():
    compilable_includes = scan_compilable_files()

    file_contents = {}
    compilable_files = list(compilable_includes.keys())
    all_library_files = list(compilable_includes.keys())

    # find the subset of files that are not #included at global scope,
    # but which are in our rosetta libraries
    for file in compilable_files:
        for included in compilable_includes[file]:
            if included not in compilable_files:
                toks = included.split("/")
                if len(toks) > 1 and toks[0] in directories_with_hhfiles_to_examine():
                    all_library_files.append(included)
    for file in all_library_files:
        file_contents[file] = open(file).readlines()
    return compilable_files, compilable_includes, file_contents


def find_library_files():
    compilable_includes = scan_compilable_files()
    all_library_files = list(compilable_includes.keys())
    for file in list(compilable_includes.keys()):
        for included in compilable_includes[file]:
            if included not in compilable_includes:
                toks = included.split("/")
                if len(toks) > 1 and toks[0] in directories_with_hhfiles_to_examine():
                    all_library_files.append(included)
    return all_library_files


def follow_includes_for_file(
    file, file_contents, cr, lines, re_pound_include, already_included
):
    cr.push_new_file(file)
    for line in file_contents[file]:
        cr.examine_line(line)
        if cr.commentless_line == "\n":
            continue
        if cr.line_is_visible():
            if re_pound_include.match(cr.commentless_line):
                include = include_for_line(cr, cr.commentless_line)
                if include in file_contents:
                    if include in already_included:
                        continue
                    already_included.append(include)
                    follow_includes_for_file(
                        include,
                        file_contents,
                        cr,
                        lines,
                        re_pound_include,
                        already_included,
                    )
                    continue
        lines.append(cr.commentless_line)
    cr.pop_last_file()
    return


def expand_includes_for_file(file, file_contents):
    cr = code_reader.CodeReader()
    lines = []
    re_pound_include = re.compile("#include")
    already_included = []
    follow_includes_for_file(
        file, file_contents, cr, lines, re_pound_include, already_included
    )
    return lines


# some aparent dependencies are not truely circular
# we have to remove them.  This returns a list of ordered pairs
# where the first file #includes the second file, but this dependency
# should not be counted.
def known_circular_dependencies():
    circular_dependencies = [
        ("utility/vector1_bool.hh", "utility/vector1.hh"),
        ("utility/vectorL_bool.hh", "utility/vectorL.hh"),
        ("utility/vector0_bool.hh", "utility/vector0.hh"),
        ("utility/io/zipstream.ipp", "utility/io/zipstream.hpp"),
        (
            "core/io/silent/ProteinSilentStuct.tmpl.hh",
            "core/io/silent/ProteinSilentStruct.hh",
        ),
        (
            "core/scoring/etable/BaseEtableEnergy.hh",
            "core/scoring/etable/BaseEtableEnergy.tmpl.hh",
        )
    ]
    return circular_dependencies


def look_for_destructorless_classes(fname, filelines):
    cr = code_reader.CodeReader()
    cr.push_new_file(fname)
    classname = ""
    class_scope_level = 0
    last_scope_level = 0
    dstor_regex = None
    class_dec = re.compile("\s*class ")
    class_has_begun = False
    dstorless_classes = []
    for line in filelines:
        cr.examine_line(line)
        # print(cr.scope_level, classname, cr.commentless_line,)
        if class_dec.match(cr.commentless_line):
            class_scope_level = last_scope_level
            if cr.scope_level > class_scope_level:
                class_has_begun = True
            classname = (
                cr.commentless_line.split("class")[1]
                .split(":")[0]
                .split("{")[0]
                .strip()
            )
            # print(classname)
            dstor_regex = re.compile("~" + classname + "()")
        if dstor_regex:
            if dstor_regex.search(cr.commentless_line):
                dstor_regex = None
        if class_has_begun and cr.scope_level == class_scope_level:
            # print("End class")
            if dstor_regex:
                dstorless_classes.append(classname)
            class_has_begun = False
            classname = ""
            dstor_regex = None
        if not class_has_begun and dstor_regex and class_scope_level < cr.scope_level:
            class_has_begun = True
        last_scope_level = cr.scope_level
    return dstorless_classes


def look_for_classes_with_header_declared_dstors(fname, filelines):
    cr = code_reader.CodeReader()
    cr.push_new_file(fname)
    classname = ""
    class_scope_level = 0
    last_scope_level = 0
    dstor_regex1 = None
    dstor_regex2 = None
    class_dec = re.compile("\s*class ")
    class_has_begun = False
    inheader_dstor_classes = []
    for line in filelines:
        cr.examine_line(line)
        # print(cr.scope_level, classname, cr.commentless_line,)
        if class_dec.match(cr.commentless_line):
            class_scope_level = last_scope_level
            if cr.scope_level > class_scope_level:
                class_has_begun = True
            classname = (
                cr.commentless_line.split("class")[1]
                .split(":")[0]
                .split("{")[0]
                .strip()
            )
            # print(classname)
            dstor_regex1 = re.compile("~" + classname + "\(\)[^;]*$")
            dstor_regex2 = re.compile("~" + classname + "\(\)\s*\{")
        if dstor_regex1:
            if dstor_regex1.search(cr.commentless_line):
                # print("MATCH")
                dstor_regex1 = None
                dstor_regex2 = None
        if dstor_regex2:
            if dstor_regex2.search(cr.commentless_line):
                # print("MATCH")
                dstor_regex1 = None
                dstor_regex2 = None
        if class_has_begun and cr.scope_level == class_scope_level:
            # print("End class")
            if not dstor_regex1:
                inheader_dstor_classes.append(classname)
            class_has_begun = False
            classname = ""
            dstor_regex1 = None
            dstor_regex2 = None
        if not class_has_begun and dstor_regex1 and class_scope_level < cr.scope_level:
            class_has_begun = True
        last_scope_level = cr.scope_level
    return inheader_dstor_classes


def find_class_dec_lines(fname, filelines):
    cr = code_reader.CodeReader()
    cr.push_new_file(fname)
    classname = ""
    class_scope_level = 0
    last_scope_level = 0
    class_dec_lines = []
    class_dec = re.compile("\s*class ")
    class_has_begun = False
    for line in filelines:
        cr.examine_line(line)
        # print(cr.scope_level, classname, cr.commentless_line,)
        if class_dec.match(cr.commentless_line):
            class_scope_level = last_scope_level
            if cr.scope_level > class_scope_level:
                class_has_begun = True
            classname = (
                cr.commentless_line.split("class")[1]
                .split(":")[0]
                .split("{")[0]
                .strip()
            )
            class_dec_lines.append((classname, cr.curr_line()))
            # print(fname, classname, cr.curr_line())
        if class_has_begun and cr.scope_level == class_scope_level:
            # print("End class")
            class_has_begun = False
            classname = ""
        if not class_has_begun and class_scope_level < cr.scope_level:
            class_has_begun = True
        last_scope_level = cr.scope_level
    return class_dec_lines


def look_for_op_containing_classes(fname, filelines):
    cr = code_reader.CodeReader()
    cr.push_new_file(fname)
    classname = ""
    class_scope_level = 0
    last_scope_level = 0
    op_data_regex = re.compile("\s*[\w:]+OP\s+\w+;")
    # op_data_regex = re.compile( "\s*public:" )
    class_dec = re.compile("\s*class ")
    class_has_begun = False
    classes_containing_an_op = []
    class_contains_no_ops = True
    for line in filelines:
        cr.examine_line(line)
        # print(cr.scope_level, classname, class_has_begun, class_scope_level, class_contains_no_ops, cr.commentless_line, )
        if class_dec.match(cr.commentless_line):
            class_scope_level = last_scope_level
            if cr.scope_level > class_scope_level:
                class_has_begun = True
            classname = (
                cr.commentless_line.split("class")[1]
                .split(":")[0]
                .split("{")[0]
                .strip()
            )
            # print(classname)
            class_contains_no_ops = True
        if class_has_begun and class_contains_no_ops:
            if op_data_regex.match(cr.commentless_line):
                # print("MATCH",)
                # print(cr.commentless_line,)
                classes_containing_an_op.append(classname)
                class_contains_no_ops = False
        if class_has_begun and cr.scope_level == class_scope_level:
            # print("End class")
            class_has_begun = False
            # print("class", classname, "ended on line", cr.curr_line())
            classname = ""
        if not class_has_begun and class_scope_level < cr.scope_level:
            class_has_begun = True
        last_scope_level = cr.scope_level
    return classes_containing_an_op


def read_rosetta_lib_class_declarations():
    dirs_to_search = ["protocols", "core", "devel", "apps"]
    regexes = [(dir, re.compile(dir)) for dir in dirs_to_search]

    compilable_files, all_includes, file_contents = load_source_tree()
    header_subset = non_fwd_hh_subset(compilable_files)
    all_classes = []
    for header in header_subset:
        keep = False
        for dir, regex in regexes:
            if regex.match(header):
                keep = True
        if not keep:
            continue

        # print("Examining header:", header)
        flines = open(header).readlines()
        classes = code_reader.read_classes_from_header(header, flines)
        all_classes.extend(classes)
    return all_classes
