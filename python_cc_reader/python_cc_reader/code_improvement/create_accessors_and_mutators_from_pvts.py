from . import blargs


def write_accessors(classname, pvts):
    new_hh_lines = []
    new_cc_lines = []
    for line in pvts:
        cols = [x.strip() for x in line.split()]
        varname = cols[-1].replace(";", "")
        accessor_name = varname if varname[-1] != "_" else varname[:-1]
        new_hh_lines.append(
            "\t" + " ".join(cols[:-1]) + " " + accessor_name + "() const;\n"
        )
        new_cc_lines.append(" ".join(cols[:-1]) + "\n")
        new_cc_lines.append(classname + "::" + accessor_name + "() const\n")
        new_cc_lines.append("{\n")
        new_cc_lines.append("\treturn " + varname + ";\n")
        new_cc_lines.append("}\n")
        new_cc_lines.append("\n")
    return new_hh_lines, new_cc_lines


def write_mutators(classname, pvts):
    new_hh_lines = []
    new_cc_lines = []
    for line in pvts:
        cols = [x.strip() for x in line.split()]
        varname = cols[-1].replace(";", "")
        mutator_name = varname if varname[-1] != "_" else varname[:-1]
        var_decl = " ".join(cols[:-1]) + (
            " const & setting" if len(cols) > 2 else " setting"
        )
        new_hh_lines.append("\tvoid " + mutator_name + "( " + var_decl + " );\n")
        new_cc_lines.append("void\n")
        new_cc_lines.append(classname + "::" + mutator_name + "( " + var_decl + " )\n")
        new_cc_lines.append("{\n")
        new_cc_lines.append("\t" + varname + " = setting;\n")
        new_cc_lines.append("}\n")
        new_cc_lines.append("\n")
    return new_hh_lines, new_cc_lines


if __name__ == "__main__":

    with blargs.Parser(locals()) as p:
        p.str("fname").required()
        p.str("classname").required()

    lines = open(fname).readlines()

    acc_new_hh_lines, acc_new_cc_lines = write_accessors(classname, lines)
    mut_new_hh_lines, mut_new_cc_lines = write_mutators(classname, lines)

    print("header:")
    for line in acc_new_hh_lines + mut_new_hh_lines:
        print(line, end=" ")
    print()
    print("ccfile:")
    for line in acc_new_cc_lines + mut_new_cc_lines:
        print(line, end=" ")
    print()
