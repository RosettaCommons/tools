from ..cpp_parser.code_utilities import *
from ..cpp_parser.code_reader import *
import re

OP_regex1 = re.compile("OP")
OP_regex2 = re.compile("owning_ptr")


def class_lacks_a_cc_dstor(cl):
    dstor_funcs = [func for func in cl.functions if func.name == "Destructor"]
    if dstor_funcs:
        # print cl.name,"has dstor; defined in header?", dstor_funcs[0].definition_present
        return dstor_funcs[0].definition_present
    return True


def class_contains_an_OP_datamember(cl):
    for dmemb in cl.data_members:
        if OP_regex1.search(dmemb[0]) or OP_regex2.search(dmemb[0]):
            print("   Found suspected OP", dmemb[0], "in class ", cl.name)
            return True
    return False


def class_needs_a_cc_dstor(cl):
    return class_contains_an_OP_datamember(cl) and class_lacks_a_cc_dstor(cl)


def class_lacks_a_cc_copy_ctor(cl):
    ctor_funcs = [func for func in cl.functions if func.name == "Constructor"]
    for ctor in ctor_funcs:
        if len(ctor.parameters) != 1:
            continue
        if ctor.parameters[0].find(cl.name) != -1:
            return ctor.definition_present
    return True


def class_declares_a_copy_ctor(cl):
    ctor_funcs = [func for func in cl.functions if func.name == "Constructor"]
    for ctor in ctor_funcs:
        if len(ctor.parameters) != 1:
            continue
        if ctor.parameters[0].find(cl.name) != -1:
            return True
    return False


def find_all_classes_that_require_cc_dstors():
    compilable_files, all_includes, file_contents = load_source_tree()

    hh_files = non_fwd_hh_subset(compilable_files)

    classes_needing_ccdstors = []

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
            print(" Processing class", cl.name)
            if class_needs_a_cc_dstor(cl):
                classes_needing_ccdstors.append(cl)
            if class_lacks_a_cc_copy_ctor(cl):
                if class_declares_a_copy_ctor(cl):
                    print(
                        "  class",
                        cl.name,
                        "declares a copy constructor and implements it in its header file",
                    )
                else:
                    print(
                        "  class",
                        cl.name,
                        "lacks a CC copy constructor declaration and implementation",
                    )
            else:
                print("  class", cl.name, "has a CC copy constructor")

    return classes_needing_ccdstors


if __name__ == "__main__":
    classes_needing_ccdstors = find_all_classes_that_require_cc_dstors()
    for cl in classes_needing_ccdstors:
        print(
            "CCDstor needed for class",
            cl.name,
            "declared in file",
            cl.file_declared_within,
        )
