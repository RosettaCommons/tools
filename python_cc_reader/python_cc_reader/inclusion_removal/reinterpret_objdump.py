# This file will break up the results of a call to objdump
# into sections, re-label all the elements in the section
# and then output each section

import re, sys


def labeled_instructions():
    return ["jne", "jmpq", "jmp", "callq", "je", "jle", "jp"]


def ignore_instructions():
    return ["nop"]


def regexes_for_instructions():
    regexes = {}
    for label in labeled_instructions():
        regexes[label] = re.compile("\s" + label + "\s")
    return regexes


def skip_sections():
    return ["static_initialization_and_destruction", "boost3arg", "GLOBAL__I"]


def relabel_sections(alllines):
    sections_to_skip = skip_sections()
    goodlines = []
    goodsections = {}
    disstart = re.compile("Disassembly")
    re_sect_header = re.compile("^[0-9a-f]* <(.*)>:$")
    re_half_sect_header = re.compile("<.*>:$")
    inside_sect_body = False
    splitter = re.compile("\S+")
    ignore_commands = ignore_instructions()

    goodsection = True
    for line in alllines:
        if disstart.search(line):
            continue
        if inside_sect_body:
            if line == "\n":
                inside_sect_body = False
                continue
            if goodsection:
                columns = line.split("\t")
                # print columns
                if len(columns) < 3:
                    goodsections[currsect].append(line.split("#")[0].strip() + "\n")
                elif columns[2].find("nop") == -1:
                    goodsections[currsect].append(line.split("#")[0].strip() + "\n")
        elif re_sect_header.match(line):
            currsect = line.split("<")[1].split(">")[0]
            inside_sect_body = True
            goodsection = True
            for skip in sections_to_skip:
                if line.find(skip) != -1:
                    goodsection = False
                    break
            if goodsection:
                goodsections[currsect] = []
                goodsections[currsect].append(splitter.findall(line)[1] + "\n")

    sorted_sections = list(goodsections.keys())
    sorted_sections.sort()
    for section in sorted_sections:
        for line in goodsections[section]:
            goodlines.append(line)

    newlines = []
    orig_labels = {}
    count = 0
    regexes = regexes_for_instructions()

    for line in goodlines:
        count += 1
        if re_half_sect_header.match(line):
            currsect = line.split("<")[1].split(">")[0]
            orig_labels[currsect] = {}
            continue
        orig_labels[currsect][line.split(":")[0].strip()] = str(count)

    count = 0
    currsect = "NONE"
    for line in goodlines:
        skip_this_line = False
        if re_half_sect_header.match(
            line
        ):  # no processing needs to be done for this line
            newlines.append(line)
            currsect = line.split("<")[1].split(">")[0]
            continue
        newline = (
            orig_labels[currsect][line.partition(":")[0].strip()] + ":" + line.partition(":")[2]
        )
        for label in regexes:
            if regexes[label].search(newline):
                copy = newline
                newline_split = newline.split(label)
                toks = splitter.findall(newline_split[1])
                if len(toks) < 2:
                    break  # keep this line as is! e.g. 4fe:       ff 95 b0 f8 ff ff       callq  *-0x750(%rbp)
                refsect = line.split("<")[1].split(">")[0].split("+")[0]
                if refsect not in orig_labels:
                    # print "skipping refsect", refsect
                    skip_this_line = True
                    break
                    # print refsect, "not present in keys for orig_labels:",
                    # for key in orig_labels.keys() :
                    #   print key,
                    # sys.exit(1)
                if toks[0].strip() in orig_labels[refsect]:
                    # print newline,
                    newline = (
                        newline_split[0]
                        + label
                        + " " * (7 - len(label))
                        + orig_labels[refsect][toks[0].strip()]
                        + " "
                        + "".join(toks[1:])
                        + "\n"
                    )
                    # print newline,
                else:
                    # print toks[ 0 ], "on line", count, "not found in orig_lables dictionary for reference section", refsect
                    newline = (
                        newline_split[0]
                        + label
                        + " " * (7 - len(label))
                        + "funcend"
                        + " "
                        + "".join(toks[1:])
                        + "\\n"
                    )
                break
        if not skip_this_line:
            newlines.append(newline)
    return newlines


def compare_objdump_lines(lines1, lines2):
    if len(lines1) != len(lines2):
        return False
    offset_found = False
    offsets = set()
    good = True
    count_mismatches = 0  # tolerate 2 mismatches
    for i in range(len(lines1)):
        if lines1[i] == lines2[i]:
            continue
        split1 = lines1[i].split("\t")
        split2 = lines2[i].split("\t")
        if len(split1) < 3 or len(split2) < 3:
            return False
        if split1[2].find("mov") == -1 or split2[2].find("mov") == -1:
            return False
        if split1[2].find("$") == -1 or split2[2].find("$") == -1:
            return False
        addr1 = split1[2].split("$")[1].split(",")[0]
        if addr1 == "":
            return False
        addr2 = split2[2].split("$")[1].split(",")[0]
        if addr2 == "":
            return False
        iad1 = int(addr1, 16)
        iad2 = int(addr2, 16)
        if offset_found and iad1 - iad2 not in offsets:
            # print "Offset mismatch:", offset, iad1 - iad2
            count_mismatches += 1
            if count_mismatches > 2:
                good = (
                    False
                )  # see how many offset mismatches we have, this is a "bad" comparison, but don't stop counting yet
        offset_found = True
        offsets.add(iad1 - iad2)
    if not good:
        print("objdump comparison failed: n_offsets", count_mismatches)
    return good


if __name__ == "__main__":

    newlines = relabel_sections(open(sys.argv[1]).readlines())
    newlines2 = relabel_sections(open(sys.argv[2]).readlines())

    verbose = False
    if len(sys.argv) > 3:
        if sys.argv[3] == "verbose":
            verbose = True

    if compare_objdump_lines(newlines, newlines2):
        # if newlines == newlines2 :
        print("They match")

    else:
        print("They don't match")
        if not verbose:
            sys.exit(1)
        if len(newlines) != len(newlines2):
            print("unequal length", len(newlines), len(newlines2))
            smaller = len(newlines)
            if len(newlines2) < len(newlines):
                smaller = len(newlines2)
            for i in range(smaller):
                print("Compare")
                print(newlines[i], end=" ")
                print(newlines2[i])
        else:
            for i in range(len(newlines)):
                if newlines[i] != newlines2[i]:
                    print(newlines[i], end=" ")
                    print("does not match")
                    print(newlines2[i])
