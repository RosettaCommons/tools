from code_utilities import load_source_tree, directories_with_ccfiles_to_examine
from add_headers import write_file

compilable_files, all_includes, file_contents = load_source_tree()
for fname in file_contents :
    toks = fname.split("/")
    if len(toks) < 0 or toks[0] not in directories_with_ccfiles_to_examine() :
        continue
    print "processing", fname
    flines = file_contents[ fname ]
    newlines = []
    found_auto_headers = False
    auto_header_begin = -1
    count_line = 0
    for line in flines :
        if found_auto_headers :
            if line == "\n" and (count_line == auto_header_begin + 1 or count_line == auto_header_begin + 2 ) :
                count_line += 1
                continue
            else :
                found_auto_headers = False
        elif line == "//Auto Headers\n" :
            found_auto_headers = True
            auto_header_begin = count_line
            count_line += 1
            continue
        newlines.append( line )
        count_line += 1
    write_file( fname, newlines )
