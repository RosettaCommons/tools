from ..cpp_parser.code_utilities import *
from ..cpp_parser.code_reader import *

# from inclusion_graph import *

compilable_files, all_includes, file_contents = load_source_tree()

hh_files = non_fwd_hh_subset(compilable_files)

for hh_file in hh_files:
    if hh_file == "devel/simple_options/option.hh":
        continue
    if hh_file.partition("utility")[2]:
        continue
    if hh_file.partition("numeric")[2]:
        continue
    if hh_file.partition("ObjexxFCL")[2]:
        continue
    print("Processing header", hh_file)
    classes = read_classes_from_header(hh_file, file_contents[hh_file])

    for cl in classes:
        has_copyctor = False
        copyctor_inhh = False
        for func in cl.functions:
            if (
                func.name == "Constructor"
                and len(func.parameters)
                and re.search("[\s]*" + cl.name + "[\s]*", func.parameters[0])
            ):
                has_copyctor = True
                if func.definition_present:
                    copyctor_inhh = True
                break
        if not has_copyctor:
            continue
        if not copyctor_inhh:
            ccfile = hh_file[:-2] + "cc"
        else:
            ccfile = hh_file
        if ccfile not in file_contents:
            print(
                "Skipping class",
                cl.name,
                "because",
                ccfile,
                "is not in the source tree",
            )
            continue
        uncopied = find_data_members_unassigned_in_copy_constructor(
            ccfile, file_contents[ccfile], cl
        )
        for un in uncopied:
            print(
                "Data member",
                un,
                "in class",
                cl.name,
                "is unassigned in the copy constructor",
            )
        # print cl.name, hh_file
        # for func in cl.functions :
        #    print "   ", cl.name, func.name
