# This class will take a ClassDeclaration object and will
# generate C++ that performs (shallow) copies of all the
# data members in both copy-ctors and assignment operators


class CodeWriter:
    def __init__(self):
        pass

    def rosetta_copyright(self):
        lines = []
        lines.append(
            "// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-\n"
        )
        lines.append("// vi: set ts=2 noet:\n")
        lines.append("//\n")
        lines.append("// (c) Copyright Rosetta Commons Member Institutions.\n")
        lines.append(
            "// (c) This file is part of the Rosetta software suite and is made available under license.\n"
        )
        lines.append(
            "// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.\n"
        )
        lines.append(
            "// (c) For more information, see http://www.rosettacommons.org. Questions about this can be\n"
        )
        lines.append(
            "// (c) addressed to University of Washington CoMotion, email: license@uw.edu.\n"
        )
        return lines

    def dstor_for_class(self, cdec):
        cname = cdec.name
        lines = []
        lines.append(cname + "::~" + cname + "() {}\n")
        return lines

    def copy_ctor_for_class(self, cdec):
        cname = cdec.name
        lines = []
        if len(cdec.data_members):
            lines.append(cname + "::" + cname + "( " + cname + " const & src ) :\n")
            for dmemb in cdec.data_members[:-1]:
                lines.append("\t" + dmemb[1] + "( src." + dmemb[1] + " ),\n")
            lines.append(
                "\t"
                + cdec.data_members[-1][1]
                + "( src."
                + cdec.data_members[-1][1]
                + " )\n"
            )
        else:
            lines.append(cname + "::" + cname + "( " + cname + " const & src )\n")
        lines.append("{}\n")
        return lines

    def assignment_operator_for_class(self, cdec):
        cname = cdec.name
        lines = []
        lines.append(cname + " const &\n")
        lines.append(cname + "::operator = ( " + cname + " const & rhs )\n")
        lines.append("{\n")
        lines.append("\t" + "if ( this != & rhs ) {\n")
        for dmemb in cdec.data_members:
            lines.append("\t\t" + dmemb[1] + " = rhs." + dmemb[1] + ";\n")
        lines.append("\t}\n")
        lines.append("\t" + "return *this;\n")
        lines.append("}\n")
        return lines

    def print_dstor_for_class(self, cdec):
        cctor_lines = self.dstor_for_class(cdec)
        for line in cctor_lines:
            print(line, end=" ")

    def print_copy_ctor_for_class(self, cdec):
        cctor_lines = self.copy_ctor_for_class(cdec)
        for line in cctor_lines:
            print(line, end=" ")

    def print_assignment_operator_for_class(self, cdec):
        assignment_op_lines = self.assignment_operator_for_class(cdec)
        for line in assignment_op_lines:
            print(line, end=" ")

    def creator_class_declaration(self, classname, classtype, baseclass):
        lc_classtype = classtype.lowercase()
        lines = []
        creator_name = classname + "Creator"
        lines.append("class " + creator_name + " : public " + baseclass + " {\n")
        lines.append("public:\n")
        lines.append(
            "\tvirtual " + classtype + "OP create_" + lc_classtype + "() const;\n"
        )
        lines.append("\tvirtual std::string keyname() const;\n")
        lines.append("\tstatic std::string " + lc_classtype + "_name();\n")
        lines.append("};\n")
        lines.append("\n")
        return lines

    def mover_creator_header_file_for_class(self, cdec):
        creator_name = cdec.name + "Creator"
        lines = self.rosetta_copyright()
        lines.append("\n")
        scopes = cdec.scope.split("::")[:-1]
        lines.append(
            "///@file " + cdec.scope.replace("::", "/") + creator_name + ".hh\n"
        )
        lines.append(
            "///@brief This class will create instances of Mover "
            + cdec.name
            + " for the MoverFactory\n"
        )
        lines.append(
            "///@author Andrew Leaver-Fay via code_writer.py (aleaverfay@gmail.com)\n"
        )
        lines.append("\n")
        included_string = "INCLUDED"
        for scope in scopes:
            included_string += "_" + scope
        included_string += "_" + cdec.name + "Creator_HH"
        lines.append("#ifndef " + included_string + "\n")
        lines.append("#define " + included_string + "\n")
        lines.append("\n")
        lines.append("#include <protocols/moves/MoverCreator.hh>\n")
        lines.append("\n")
        for scope in scopes:
            lines.append("namespace " + scope + " {\n")
        lines.append("\n")
        classlines = self.creator_class_header(
            cdec.name, "Mover", "protocols::moves::MoverCreator"
        )
        # lines.append("class "+creator_name+" : public protocols::moves::MoverCreator {\n")
        # lines.append("public:\n")
        # lines.append("\tvirtual MoverOP create_mover() const;\n" )
        # lines.append("\tvirtual std::string keyname() const;\n")
        # lines.append("\tstatic std::string mover_name();\n")
        # lines.append("};\n")
        # lines.append("\n")
        lines.extend(classlines)
        for scope in scopes:
            lines.append("}\n")
        lines.append("\n")
        lines.append("#endif\n")
        lines.append("\n")

        return lines
