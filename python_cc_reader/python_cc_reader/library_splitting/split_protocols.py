from library_levels import protocols_levels
from optparse import OptionParser
import pprint
import shutil


def get_letter_index(character):
    char_dict = {
        "a": 0,
        "b": 1,
        "c": 2,
        "d": 3,
        "e": 4,
        "f": 5,
        "g": 6,
        "h": 7,
        "i": 8,
        "j": 9,
    }
    return char_dict[character]


def get_index_letter(index):
    index_dict = {
        0: "a",
        1: "b",
        2: "c",
        3: "d",
        4: "e",
        5: "f",
        6: "g",
        7: "h",
        8: "i",
        9: "j",
    }
    return index_dict[index]


def get_header():
    header = """# -*- mode:python;indent-tabs-mode:nil;show-trailing-whitespace:t; -*-
            #
            # (c) Copyright Rosetta Commons Member Institutions.
            # (c) This file is part of the Rosetta software suite and is made available under license.
            # (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
            # (c) For more information, see http://www.rosettacommons.org. Questions about this can be
            # (c) addressed to University of Washington UW CoMotion, email: license@uw.edu.
            
            # Project settings for rosetta sources
            
            """
    return header


def setup_map_subset(source_map, top_level_namespace_set):
    new_source_map = {}
    for key in source_map:
        key_vector = key.split("/")
        if len(key_vector) < 2:
            continue
        if key_vector[0] == "protocols" and key_vector[1] in top_level_namespace_set:
            new_source_map[key] = source_map[key]
    return new_source_map


def setup_map_inverse_subset(source_map, top_level_namespace_set):
    new_source_map = {}
    for key in source_map:
        key_vector = key.split("/")
        if len(key_vector) < 2:
            new_source_map[key] = source_map[key]
            continue
        if (
            key_vector[0] == "protocols"
            and key_vector[1] not in top_level_namespace_set
        ):
            new_source_map[key] = source_map[key]
    return new_source_map


def pprint_object_as_statement(object, name, outfile):
    outfile.write(name + " = ")
    pp.pprint(object)


def get_upstream_subprojects(level_index, library_levels):
    subprojects = [
        "core.5",
        "core.4",
        "core.3",
        "core.2",
        "core.1",
        "basic",
        "numeric",
        "utility",
        "ObjexxFCL",
        "z",
        "sqlite3",
        "cppdb",
    ]
    if level_index == 1:
        return subprojects
    for index in range(level_index):
        sublevel_count = len(library_levels[index])
        index += 1
        if index == level_index:
            break
        if sublevel_count == 1:
            subprojects.insert(0, "protocols." + str(index))
        else:
            for sub_index in range(sublevel_count):
                subprojects.insert(
                    0,
                    "protocols_" + str(get_index_letter(sub_index)) + "." + str(index),
                )

    return subprojects


if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option(
        "--make_src_settings_for_level",
        dest="level_to_make",
        help="Make a src.settings file for the desired protocols level",
        default="",
    )
    (options, args) = parser.parse_args()

    # parsed_level_name = options.level_to_make.split(".")

    library_levels = protocols_levels()

    library_names = [
        (1,),
        (2,),
        (3,),
        (4, "a"),
        (4, "b"),
        (4, "c"),
        (4, "d"),
        (4, "e"),
        (4, "f"),
        (5, "a"),
        (5, "b"),
        (5, "c"),
        (6, "a"),
        (6, "b"),
        (6, "c"),
    ]

    for name in library_names:
        major_library_level = library_levels[name[0] - 1]
        name_string = ""
        if len(name) == 2:
            x, y = name
            # print name[0]
            print(major_library_level)
            minor_library_level = major_library_level[get_letter_index(y) - 1]
            name_string = "_" + y + "." + str(x)
        else:
            print(name)
            minor_library_level = major_library_level[0]
            print(minor_library_level)
            name_string = "." + str(name[0])

        protocols_src_settings_map = {}
        exec(
            compile(
                open("protocols.src.settings", "rb").read(),
                "protocols.src.settings",
                "exec",
            ),
            protocols_src_settings_map,
        )

        source_file_subset = setup_map_subset(
            protocols_src_settings_map["sources"], minor_library_level
        )
        inverse_source_file_subset = setup_map_inverse_subset(
            protocols_src_settings_map["sources"], minor_library_level
        )

        # pp.pprint(source_file_subset)
        outfile = open("protocols" + name_string + ".src.settings", "w")
        pp = pprint.PrettyPrinter(indent=1, stream=outfile)
        outfile.write(get_header())
        pprint_object_as_statement(source_file_subset, "sources", outfile)
        outfile.write("\n")

        pprint_object_as_statement(
            protocols_src_settings_map["include_path"], "include_path", outfile
        )
        pprint_object_as_statement(
            protocols_src_settings_map["library_path"], "library_path", outfile
        )
        pprint_object_as_statement(
            protocols_src_settings_map["libraries"], "libraries", outfile
        )
        pprint_object_as_statement(
            get_upstream_subprojects(int(name[0]), library_levels),
            "subprojects",
            outfile,
        )
        outfile.close()

        shutil.move("protocols.src.settings", "protocols.src.settings.old")
        outfile_orig = open("protocols.src.settings", "w")
        pp = pprint.PrettyPrinter(indent=1, stream=outfile_orig)
        pp = pprint.PrettyPrinter(indent=1, stream=outfile_orig)
        outfile_orig.write(get_header())
        pprint_object_as_statement(inverse_source_file_subset, "sources", outfile_orig)
        outfile_orig.write("\n")
        pprint_object_as_statement(
            protocols_src_settings_map["include_path"], "include_path", outfile_orig
        )
        pprint_object_as_statement(
            protocols_src_settings_map["library_path"], "library_path", outfile_orig
        )
        pprint_object_as_statement(
            protocols_src_settings_map["libraries"], "libraries", outfile_orig
        )
        pprint_object_as_statement(
            protocols_src_settings_map["subprojects"], "subprojects", outfile_orig
        )

        outfile_orig.close()
