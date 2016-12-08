import blargs
import rewrite_rosetta_script

def append_block( real_blocks, faux_blocks, all_blocks, found_ampersand, begin_line, end_line ) :
    if not found_ampersand :
        real_blocks.append(( begin_line, end_line ))
    else :
        faux_blocks.append( (begin_line, end_line))
    all_blocks.append( (begin_line, end_line ))

def find_xml_blocks_in_lines( lines ) :
    real_blocks = []
    faux_blocks = []
    all_blocks =  []
    
    in_quote_block = False 
    first_line_of_quote_block = False
    in_indentation_block = False
    block_starts_w_lt = False
    found_ampersand = False
    for i,line in enumerate( lines ) :
        #print line, found_ampersand
        if len(line) >= 3 and line[0:3] == "```" :
            if in_indentation_block :
                in_indentaiton_block = False
                if block_starts_w_lt :
                    append_block( real_blocks, faux_blocks, all_blocks, found_ampersand, begin_line, i-1 )
            #print "found ```"
            if in_quote_block :
                if block_starts_w_lt :
                    append_block( real_blocks, faux_blocks, all_blocks, found_ampersand, begin_line, i-1 )
                in_quote_block = False
            else :
                begin_line = i+1
                in_quote_block = True
                first_line_of_quote_block = True
                block_starts_w_lt = False
                found_ampersand = False
        elif not in_quote_block and len(line) >= 3 and line[0:3] == "   " :
            if not in_indentation_block :
                in_indentation_block = True
                block_starts_w_lt = False
                found_ampersand = False
                begin_line = i
                for j in xrange(len(line)) :
                    if line[j] == "\t" or line[j] == " " : continue
                    if line[j] == "<" : block_starts_w_lt = True
                    break

            if not found_ampersand :
                if line.find( "&" ) != -1 :
                    found_ampersand = True

        elif in_indentation_block :
            in_indentation_block = False
            if block_starts_w_lt :
                append_block( real_blocks, faux_blocks, all_blocks, found_ampersand, begin_line, i-1 )
        elif in_quote_block :
            if first_line_of_quote_block :
                first_line_of_quote_block = False
                for j in xrange(len(line)) :
                    if line[j] == "\t" or line[j] == " " : continue
                    if line[j] == "<" : block_starts_w_lt = True
                    break
            if not found_ampersand :
                if line.find( "&" ) != -1 :
                    found_ampersand = True

    for pair in all_blocks :
        print "begin and end:", pair[0], pair[1]

    return real_blocks, faux_blocks, all_blocks

def rewrite_faux_xml_line( line ) :
    seen_equals = False
    seen_lparen = False
    lparen_depth = 0
    new_chars = []
    for char in line :
        #print char,
        if seen_equals :
            if not seen_lparen and char == '"' :
                # print "Already have quotes; skip this attribute"
                # if there are already quotes surrounding the attribute contents
                # don't add more!
                seen_equals = False
                new_chars.append( char )
            elif not seen_lparen and char == "(" :
                new_chars.append( '"' )
                new_chars.append( char )
                seen_lparen = True
                lparen_depth = 1
            elif seen_lparen and char == "(" :
                lparen_depth += 1
                new_chars.append( char )
            elif seen_lparen and char == ")" :
                lparen_depth -= 1
                if lparen_depth == 0 :
                    new_chars.append( char )
                    new_chars.append( '"' )
                    #reset!
                    seen_equals = False;
                    seen_lparen = False;
                else :
                    new_chars.append( char )
            elif char == '"' :
                #print "converting double quote to single quote"
                new_chars.append( "'" )
            else :
                new_chars.append( char )
            continue
        elif char == "=" :
            seen_equals = True
            seen_lparen = False
        new_chars.append( char )
    return "".join( new_chars )

def wrap_faux_xml_attributes_in_quotes( lines, faux_blocks ) :
    newlines = []
    count_xml_block = 0
    for i, line in enumerate( lines ) :
        if count_xml_block >= len( faux_blocks ) or i < faux_blocks[ count_xml_block ][0] :
            newlines.append( line )
        elif count_xml_block < len( faux_blocks ) and i >= faux_blocks[ count_xml_block ][0] and i <= faux_blocks[ count_xml_block ][1] :
            new_line = rewrite_faux_xml_line( lines[ i ] )
            newlines.append( new_line )
            if i == faux_blocks[ count_xml_block ][1 ] :
                count_xml_block += 1

    return newlines


def rewrite_xml_blocks_in_lines( lines, xml_blocks ) :
    newlines = []
    count_xml_block = 0
    for i, line in enumerate( lines ) :
        if count_xml_block >= len( xml_blocks ) or i < xml_blocks[ count_xml_block ][0] :
            newlines.append( line )
        elif count_xml_block < len( xml_blocks ) and i == xml_blocks[ count_xml_block ][0] :
            block = xml_blocks[ count_xml_block ]
            #print block
            xml_lines = lines[ block[0]:(block[1]+1) ]
            print xml_lines
            new_content = rewrite_rosetta_script.rewrite_xml_rosetta_script_lines( xml_lines )
            newlines.append( new_content )

        if count_xml_block < len( xml_blocks ) and i == xml_blocks[ count_xml_block ][ 1 ] :
            count_xml_block += 1

    return newlines

if __name__ == "__main__" :
    with blargs.Parser(locals()) as p :
        p.str( "input_file" ).required()

    lines = open( input_file ).readlines()

    xml_blocks, faux_blocks, all_blocks = find_xml_blocks_in_lines( lines )
    newlines = wrap_faux_xml_attributes_in_quotes( lines, faux_blocks )
    #for line in newlines : print line,
    newlines = rewrite_xml_blocks_in_lines( lines, xml_blocks )
    newlines = rewrite_xml_blocks_in_lines( newlines, faux_blocks )
    open( input_file, "w" ).writelines( newlines )
