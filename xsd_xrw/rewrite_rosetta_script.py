#!/usr/bin/env python

import sys, os

try:
    import blargs
except ImportError:
    # if this script is in the Rosetta/tools/xsd_xrw/ directory
    # blargs is in the ../external/ directory. Add that to the path. and re-import
    blargs_path = os.path.join( os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'external')
    sys.path.append(blargs_path)
    import blargs

#import copy

debug = False


class XMLToken :
    def __init__( self ) :
        self.line_start = -1
        self.position_start = -1
        self.line_end = -1
        self.position_end = -1
        self.contents = ""
        self.whitespace = False #only true inside of tags; content outside of tags is ignored
        self.deleted = False
        self.index = 0
        self.in_tag = False
    @staticmethod
    def tok_for_single_pos( line_number, position_on_line, in_tag ) :
        tok = XMLToken()
        tok.line_start = line_number
        tok.line_end = line_number
        tok.position_start = position_on_line
        tok.position_end = position_on_line
        tok.in_tag = in_tag
        return tok
    def uninitialized( self ) :
        return self.line_start == -1
    def is_whitespace( self ) :
        if self.in_tag : return self.whitespace
        for x in self.contents :
            if x != "\n" and x != " " and x != "\t" :
                return False
        return True

class Tag :
    def __init__( self ) :
        self.name = ""
        self.attributes = []
        self.tokens = []
        self.closed = False
        self.terminator = False # e.g. </SCOREFXNS>
        self.new_tokens_since_parsing = False

    def remove_token( self, token ) :
        #print "remove token", token.contents, ",".join( [ x.contents for x in self.tokens ] )
        self.tokens.remove( token )
        self.new_tokens_since_parsing = True

    def add_token( self, token, pos ) :
        # you cannot replace either the first or the last token
        # which you shouldn't need to, because these are the "<" and
        # ">" symbols.
        #print "add token before: ", "".join( [ x.contents for x in self.tokens ] ), pos, token.contents
        assert( pos != 0 and pos != len(self.tokens) )
        self.tokens.insert( pos, token )
        #print "add token after: ", "".join( [ x.contents for x in self.tokens ] )
        self.new_tokens_since_parsing = True

    def replace_tokens_in_list( self, token_list ) :
        #slow, don't use this
        if not self.new_tokens_since_parsing : return token_list
        for i in xrange( len( token_list ) ) :
            # look for the first token and then the last token
            if self.tokens[0] is token_list[ i ] :
                for j in xrange( i+1, len( token_list )) :
                    if self.tokens[-1] is token_list[ j ] :
                        token_list = token_list[ :i ] + self.tokens + token_list[ j+1: ]
                        return token_list
        assert( False )

def find_attribute_in_tag( tag, attribute_name ) :
    for attr in tag.attributes :
        if attr[0].contents == attribute_name :
            return attr
    return None

def name_token_of_tag( tag ) :
    for j,tok in enumerate( tag.tokens ) :
        if j == 0 : continue # skip the "<"
        for char in tok.contents :
            if char != "" and char != "\t" and char != "\n" :
                return tok
    return None

class Element :
    def __init__( self ) :
        self.sub_elements = []
        self.name = ""
        self.tags = [] # perhaps only a single tag
    def reconstitute_token_list( self, old_tok_list, new_tok_list, tok_pos ) :
        #print self.name, tok_pos, old_tok_list[ tok_pos ].contents
        assert( len( self.tags ) == 1 or len( self.tags ) == 2 )

        # append the old tokens to the new list until we reach the first token
        # that belongs to the tag; then insert all of that tag's tokens
        for which_tag in xrange( len( self.tags ) ) :
            tag_first_tok = self.tags[ which_tag ].tokens[0]
            for i in xrange( tok_pos, len( old_tok_list ) ) :
                if old_tok_list[ i ] is tag_first_tok :
                    if debug : print "Found token: ", tag_first_tok.contents, "for", self.tags[which_tag].name
                    new_tok_list.extend( self.tags[ which_tag ].tokens )
                    tag_last_tok = self.tags[ which_tag ].tokens[-1]
                    for j in xrange( i+1, len( old_tok_list ) ) :
                        if old_tok_list[ j ] == tag_last_tok :
                            tok_pos = j+1
                            if which_tag == 0 :
                                # recurse!
                                for sub_element in self.sub_elements :
                                    tok_pos, new_tok_list = sub_element.reconstitute_token_list( old_tok_list, new_tok_list, tok_pos )
                            if which_tag + 1 == len( self.tags ) :
                                if debug : print "returning", tok_pos
                                return tok_pos, new_tok_list
                            else :
                                break
                    # only reachable if we're not on the last tag
                    assert( which_tag + 1 < len( self.tags ) )
                    break
                else :
                    if debug : print self.name, "skipping", tok_pos,"\"" + old_tok_list[i].contents + "\" looking for \"" + tag_first_tok.contents + "\""
                    new_tok_list.append( old_tok_list[ i ] )

        unreachable = False
        assert( unreachable ) # we should not reach here

def tokenize_lines( lines ) :
    in_tag = False
    in_block_comment = False
    in_double_quote = False
    in_single_quote = False
    in_whitespace = False
    tokens = []

    curr_token = XMLToken()
    for i,line in enumerate(lines) :
        for j,symb in enumerate( line ) :
            #print symb,
            #print i, j, symb, "is_newline?", symb=="\n", "in tag", in_tag, "in block comment", in_block_comment, "in double quote", in_double_quote, "in whitespace", in_whitespace
            if in_block_comment :
                if symb == ">" and j >= 2 and line[j-2:j+1] == "-->" :
                    assert( not curr_token.uninitialized() )
                    curr_token.line_end = i
                    curr_token.position_end = j
                    curr_token.in_tag = in_tag
                    tokens.append( curr_token )
                    curr_token = XMLToken()
                    in_block_comment = False
            elif in_tag :
                if in_double_quote or in_single_quote:
                    if ( symb == '"' and in_double_quote ) or ( symb == "'" and in_single_quote ) :
                        curr_token.line_end = i
                        curr_token.position_end = j
                        curr_token.in_tag = in_tag
                        tokens.append( curr_token )
                        curr_token = XMLToken()
                        in_double_quote = False
                        in_single_quote = False
                        #print "finished reading quote", i, j, symb
                elif symb == ">" :
                    in_whitespace = False
                    end_of_element = j > 0 and line[ j-1:j+1 ] == "/>"
                    if end_of_element :
                        if not curr_token.uninitialized() :
                            curr_token.line_end = i if j != 1 else i-1
                            curr_token.position_end = j-2 if j != 1 else len(lines[i-1]) - 1
                            if curr_token.line_start < curr_token.line_end or ( curr_token.line_start == curr_token.line_end and curr_token.position_end >= curr_token.position_start ) :
                                #print "adding token before end of element", curr_token.line_start, curr_token.line_end,
                                #print curr_token.position_start, curr_token.position_end
                                curr_token.in_tag = in_tag
                                tokens.append( curr_token )
                            curr_token = XMLToken()
                        curr_token.line_start = i
                        curr_token.position_start = j-1
                        curr_token.line_end = i
                        curr_token.position_end = j
                        curr_token.in_tag = in_tag
                        tokens.append( curr_token )
                        curr_token = XMLToken()
                    else :
                        if not curr_token.uninitialized() :
                            curr_token.line_end = i if j != 0 else i-1
                            curr_token.position_end = j-1 if j != 0 else len(lines[i-1])-1
                            curr_token.in_tag = in_tag
                            tokens.append( curr_token )
                            curr_token = XMLToken()

                        # look at the previous two tokens; if the last token was whitespace
                        # and the one before was "/", then make this token a "/>"
                        if len(tokens) > 2 and tokens[-1].whitespace and contents_for_token( lines, tokens[-2] ) == "/" :
                            #print "FOUND IT"
                            curr_token.line_start = i; curr_token.position_start = tokens[-2].position_start
                            curr_token.line_end = i; curr_token.position_end = j

                            tokens.pop(); # destroy the last two tokens
                            tokens.pop();

                            curr_token.in_tag = in_tag
                            tokens.append( curr_token )
                            curr_token = XMLToken()
                        else :
                            #if len(tokens) > 2 :
                            #    print "closing; last two? '" + contents_for_token( lines, tokens[-2] ) + "' and '" + contents_for_token( lines, tokens[-1] ) + "'"
                            curr_token.in_tag = in_tag
                            tokens.append( XMLToken.tok_for_single_pos( i, j, in_tag ))
                    in_tag = False

                elif symb == '"' or symb == "'" :
                    in_whitespace = False
                    if not curr_token.uninitialized() :
                        curr_token.line_end = i if j != 0 else i-1
                        curr_token.position_end = j-1 if j != 0 else len(lines[i-1])-1
                        curr_token.in_tag = in_tag
                        tokens.append( curr_token )
                        curr_token = XMLToken()
                    if symb == '"' : in_double_quote = True
                    else           : in_single_quote = True
                    curr_token.line_start = i
                    curr_token.position_start = j
                    #print "started reading quote", i, j, symb
                elif symb == " " or symb == "\t" or symb == "\n" :
                    if in_whitespace :
                        pass
                    else :
                        if not curr_token.uninitialized() :
                            curr_token.line_end = i if j != 0 else i-1
                            curr_token.position_end = j-1 if j != 0 else len(lines[i-1])-1
                            curr_token.in_tag = in_tag
                            tokens.append( curr_token )
                            curr_token = XMLToken()
                        #print "Starting new non-whitespace token", i, j, symb
                        curr_token.line_start = i
                        curr_token.position_start = j
                        curr_token.whitespace = True
                        in_whitespace = True
                elif symb == "=" :
                    in_whitespace = False
                    if not curr_token.uninitialized() :
                        curr_token.line_end = i if j != 0 else i-1
                        curr_token.position_end = j-1 if j != 0 else len(lines[i-1])-1
                        curr_token.in_tag = in_tag
                        tokens.append( curr_token )
                        curr_token = XMLToken()
                    curr_token.in_tag = in_tag
                    tokens.append( XMLToken.tok_for_single_pos( i, j, in_tag ))
                else :
                    if in_whitespace :
                        curr_token.line_end = i if j != 0 else i-1
                        curr_token.position_end = j-1 if j != 0 else len(lines[i-1])-1
                        curr_token.in_tag = in_tag
                        tokens.append( curr_token )
                        curr_token = XMLToken()
                    in_whitespace = False
                    if curr_token.uninitialized() :
                        #print "Starting new non-whitespace token", i, j, symb
                        curr_token.line_start = i
                        curr_token.position_start = j
            elif symb == "<" :
                in_whitespace = False
                if j+3 < len(line) and line[j:j+4] == "<!--" :
                    if not curr_token.uninitialized() :
                        curr_token.line_end = i if j != 0 else i-1
                        curr_token.position_end = j-1 if j != 0 else len(lines[i-1]) - 1
                        curr_token.in_tag = in_tag
                        tokens.append( curr_token )
                        curr_token = XMLToken()
                    in_block_comment = True
                    curr_token.line_start = i
                    curr_token.position_start = j
                else :
                    if not curr_token.uninitialized() :
                        curr_token.line_end = i if j != 0 else i-1
                        curr_token.position_end = j-1 if j != 0 else len(lines[i-1]) - 1
                        curr_token.in_tag = in_tag
                        tokens.append( curr_token )
                        curr_token = XMLToken()
                    in_tag = True
                    curr_token.in_tag = in_tag
                    tokens.append( XMLToken.tok_for_single_pos( i, j, in_tag ))
            else :
                # outside of a tag; just nab all the text and put it into a single token
                # until we get to something interesting
                if curr_token.uninitialized() :
                    curr_token.line_start = i
                    curr_token.position_start = j
    #end for
    if not curr_token.uninitialized() :
        curr_token.line_end = len( lines ) - 1
        curr_token.position_end = len( lines[-1] ) - 1
        curr_token.in_tag = in_tag
        tokens.append( curr_token )

    for i,tok in enumerate( tokens ) :
        tok.contents = contents_for_token( lines, tok )
        tok.index = i # save this for the purpose of reconstructing the
        #print i, tok.contents

    for tok in tokens :
        # look for tokens that are "/" + whitespace + ">"
        # and change them to just be "/>"
        # because you cannot have whitespace at the close of a tag
        if len(tok.contents) > 2 and tok.contents[0] == "/" and tok.contents[-1] == ">" :
            all_whitespace = True
            for i in xrange(1,len(tok.contents)-1) :
                if tok.contents[i] == " " or tok.contents[i] == "\n" or tok.contents[i] == "\t" : continue
                all_whitespace = False
                break
            if all_whitespace :
                tok.contents = "/>"

    #for tok in tokens :
    #    print "token", tok.index, tok.whitespace, "\"" + tok.contents + "\""

    return tokens

def contents_for_token( lines, tok ) :
    tok_lines = []
    for lineno in xrange(tok.line_start, tok.line_end + 1 ) :
        if lineno == tok.line_start :
            if lineno == tok.line_end :
                tok_lines.append( lines[ lineno ][ tok.position_start:(tok.position_end+1) ] )
            else :
                tok_lines.append( lines[ lineno ][ tok.position_start: ] )
        elif lineno == tok.line_end :
            tok_lines.append( lines[ lineno ][ :( tok.position_end + 1 ) ] )
        else :
            tok_lines.append( lines[ lineno ] )
    return "".join( tok_lines )


def tokens_into_tags( tokens ) :
    # interpets tokens as XML tags.
    # marks whitespace tokens between attribute names, the equals sign, and their values
    # as deleted.
    tags = []
    in_tag = False
    in_attribute = False

    curr_tag = Tag()
    for i,tok in enumerate( tokens ) :
        #print "("+str(i), tok.contents+")"
        if in_tag :
            curr_tag.tokens.append( tok )
            if curr_tag.name == "" :
                if tok.contents[0] == "/" :
                    curr_tag.name = tok.contents[1:]
                    curr_tag.closed = True
                    # March backwards and look to see if
                    # the first non-whitespace token behind this "/" is a "<"
                    j = i-1
                    while j >= 0 :
                        if not tokens[j].is_whitespace() :
                            if tokens[j].contents == "<" :
                                curr_tag.terminator = True
                            break
                elif tok.is_whitespace() :
                    # you cannot have whitespace between the opening "<" and the tag name
                    tok.deleted = True
                else :
                    curr_tag.name = tok.contents
            elif in_attribute :
                #assert( len( tok.contents ) == 1 )
                if tok.is_whitespace() :
                    tok.deleted = True
                else :
                    curr_tag.attributes[-1].append( tok )
                    in_attribute = False
            elif tok.contents == "=" :
                in_attribute = True
                # now we need to go backwards to find the first non-whitespace token before
                # the equals sign; mark all of the whitespace tokens in between as deleted.
                j = i-1
                while ( j >= 0 ) :
                    if not tokens[j].is_whitespace() :
                        curr_tag.attributes.append( [ tokens[j] ] )
                        break
                    else :
                        tokens[j].deleted = True
                    j -= 1
            elif tok.contents == "" :
                print "Whoa, how did we end up here?", tok.index
            elif tok.contents[-1] == ">" :

                if tok.contents == "/>" :
                    curr_tag.closed = True
                if debug : print i, tok.contents, ", ".join( [ x.contents for x in curr_tag.tokens ] )
                tags.append( curr_tag )
                curr_tag = Tag()
                in_tag = False
        elif tok.contents == "<" :
            in_tag = True
            curr_tag.tokens.append( tok )

    elements = []
    roots = []
    for i,tag in enumerate( tags ) :
        if debug : print i, tag.name
        if len( elements ) > 0 :
            if elements[-1].name != tag.name :
                if tag.terminator :
                    print "Input script is not valid XML"
                    print "terminating tag \"" + name_token_of_tag( tag ).contents + "\" does not match the wrapping tag \"" + elements[-1].name + "\""
                    if len(elements) >= 2 and elements[-2].name == tag.name :
                        print "Perhaps you need to close the", elements[-1].name + " tag?"
                    sys.exit( 1 )
                elements.append( Element() )
                elements[-1].name = tag.name
                elements[-2].sub_elements.append( elements[-1] )
            if len(elements[-1].tags) >= 2 :
                # too many tags for a single element!
                print "Input script is not valid XML"
                print "This is the third tag in a row for \"" + elements[-1].name + "\" found on line", \
                    str( name_token_of_tag( tag ).line_start+1 ) + "; there should only ever be two. Did" + \
                    " you forget to close one of these tags?"
                sys.exit(1)
            elements[-1].tags.append( tag )
            #print "appending tag", tag.name, "to", elements[-1].name, len(elements[-1].tags)
        else :
            elements.append( Element() )
            elements[-1].tags.append( tag )
            elements[-1].name = tag.name
            roots.append( elements[-1] )
            #print "appending root tag"
        if tag.closed :
            elements.pop()
    return tags, roots

# def print_element_tree( root, depth=0 ) :
#     print ( " "*depth + root.name )
#     for elem in root.sub_elements :
#         print_element_tree( elem, depth+1 )

# ok -- take a set of tags, and edit them so that the rhs of each attribute statement
# is surrounded by quotes
def surround_attributes_w_quotes( tags ) :
    for tag in tags :
        for attr in tag.attributes :
            rhs = attr[1]
            if rhs.contents[0] == '"' and rhs.contents[-1] == '"' : continue
            if rhs.contents[0] == "'" and rhs.contents[-1] == "'" : continue
            #print rhs.contents
            assert (( rhs.contents[0] != '"' and rhs.contents[0] != "'" ) or
                    rhs.contents[0] == rhs.contents[-1] )
            rhs.contents = '"' + rhs.contents + '"'

def replace_element_name_w_attribute( element, new_name, new_attribute ) :
    element.name = new_name
    for i,tag in enumerate( element.tags ) :
        old_name = tag.name
        tag.name = new_name
        # find the name token and remove it
        done = False
        for j,tok in enumerate( tag.tokens ) :
            if j == 0 : continue
            for char in tok.contents : # skip the "<"
                if char != " " and char != "\t" and char != "\n" :
                    tag.remove_token( tok )
                    done = True
                    break
            if done : break
        #and now add a token for the name
        name_token = XMLToken()
        name_token.line_start = tag.tokens[0].line_start
        name_token.line_end = tag.tokens[0].line_end
        name_token.position_start = tag.tokens[0].position_start + 1
        name_token.position_end = tag.tokens[0].position_start + 1 # garbage
        name_token.contents = new_name if i == 0 else ("/" + new_name )
        tag.add_token( name_token,  1 )
        if i == 0 :
            attr_tokens = [ XMLToken(), XMLToken(), XMLToken(), XMLToken() ]
            spellings = [ " ", new_attribute, "=", '"' + old_name + '"' ]
            for j,tok in enumerate( attr_tokens ) :
                tok.line_start = tok.line_end = name_token.line_start
                tok.position_start = tok.position_end = name_token.position_start # garbage
                tok.contents = spellings[j]
                tag.add_token( tok, 2+j )

def swap_element_name_w_attribute( sub_element, old_attribute_name, new_attribute_name ) :
    old_name = sub_element.name
    for i,tag in enumerate( sub_element.tags ) :
        if i == 0 :
            attr = find_attribute_in_tag( tag, old_attribute_name )
            if not attr : return
            attr[0].contents = new_attribute_name
            old_attr_val = attr[1].contents[1:-1]
            attr[1].contents = '"' + old_name + '"'
        tag.name = old_attr_val
        name_tok = name_token_of_tag( tag )
        name_tok.contents = old_attr_val if i == 0 else ( "/" + old_attr_val )
    sub_element.name = old_attr_val



def recursively_rename_subelements( root, element_name, new_subelement_name, new_attribute_name = "name" ) :
    if root.name == element_name :
        for sub_element in root.sub_elements :
            if sub_element.name == new_subelement_name : continue # Assume this doesn't need rewriting
            if sub_element.name == "xi:include" : continue # skip these
            replace_element_name_w_attribute( sub_element, new_subelement_name, new_attribute_name )
    for element in root.sub_elements :
        recursively_rename_subelements( element, element_name, new_subelement_name, new_attribute_name )

def recursively_rename_particular_subelements( root, element_name, old_subelement_name, new_subelement_name ) :
    if root.name == element_name :
        for sub_element in root.sub_elements :
            if sub_element.name != old_subelement_name : continue
            for i,tag in enumerate( sub_element.tags ) :
                name_tok = name_token_of_tag( tag )
                name_tok.contents = new_subelement_name if i == 0 else ( "/" + new_subelement_name )
    for element in root.sub_elements :
        recursively_rename_particular_subelements( element, element_name, old_subelement_name, new_subelement_name )


def recursively_swap_attribute_and_element_name( root, element_name, old_attribute_name, new_attribute_name = "name" ) :
    if root.name == element_name :
        for sub_element in root.sub_elements :
            if sub_element.name == "xi:include" : continue
            swap_element_name_w_attribute( sub_element, old_attribute_name, new_attribute_name )

    for element in root.sub_elements :
        recursively_swap_attribute_and_element_name( element, element_name, old_attribute_name, new_attribute_name )


def rename_score_functions( root, tokens ) :
    recursively_rename_subelements( root, "SCOREFXNS", "ScoreFunction" )
    return tokens

def rename_fragments_from_frag_reader( root, tokens ) :
    # for elements beneath a "FRAGSETS" element:
    # 1. elements beneath FRAGMENTS element gets renamed as
    # "FragReader" and the old name becomes the "name" attribute.
    # 2. everything else beneath the "FRAGSETS" element
    # gets renamed "FragSet" and the old name becomes the name attribute
    # 3. The FRAGMENTS element must appear as the first child of the
    # FragSetLoader. It's also required??
    for element in root.sub_elements :
        if element.name == "FRAGSETS" :
            for fragset_element in element.sub_elements :
                if fragset_element.name == "FRAGMENTS" :
                    for fragreader_element in fragset_element.sub_elements  :
                        if fragreader_element.name == "xi:include" : continue
                        replace_element_name_w_attribute( fragreader_element, "FragReader", "name" )
                else :
                    if fragset_element.name == "xi:include" : continue
                    replace_element_name_w_attribute( fragset_element, "FragSet", "name" )
        rename_fragments_from_frag_reader( element, tokens )
    return tokens

def move_OUTPUT_as_last_child_of_ROSETTASCRIPTS( root, tokens ) :
    if root.name == "ROSETTASCRIPTS" :
        output_index = -1
        for i,elem in enumerate( root.sub_elements ) :
            if elem.name == "OUTPUT" :
                output_index = i
                break
        if output_index != -1 and output_index + 1 != len( root.sub_elements ) :
            out_elem = root.sub_elements[ output_index ]
            last_elem = root.sub_elements[ -1 ]
            prev_last = root.sub_elements[output_index-1].tags[-1].tokens[-1].index if output_index != 0 else root.tags[0].tokens[-1].index
            out_first = out_elem.tags[0].tokens[0].index
            out_last  = out_elem.tags[-1].tokens[-1].index
            next_first = root.sub_elements[output_index+1].tags[0].tokens[0].index
            last_last = root.sub_elements[-1].tags[-1].tokens[-1].index
            #last_first = last_elem.tags[0].tokens[0].index
            #last_last  = last_elem.tags[-1].tokens[-1].index
            assert( len(root.tags) == 2 )
            root_first = root.tags[1].tokens[0].index
            #print "OUTPUT REORDER:", prev_last+1, out_last+1, root_first
            #print "1", 0, out_first, tokens[0].contents, tokens[out_first].contents
            #print "2", next_first, last_last+1, tokens[next_first].contents, tokens[last_last+1].contents
            #print "3", out_last+1, next_first, tokens[out_last+1].contents, tokens[next_first].contents
            #print "4", out_first, out_last+1, tokens[out_first].contents, tokens[out_last+1].contents
            #print "5", last_last+1, "END", tokens[last_last+1].contents, "END"

            tokens = tokens[ : out_first            ] + \
                     tokens[next_first: last_last+1 ] + \
                     tokens[out_last+1: next_first  ] + \
                     tokens[out_first: out_last+1   ] + \
                     tokens[last_last+1:            ]

            #print "old subelement order:"
            #for elem in root.sub_elements :
            #    print elem.name
            root.sub_elements.remove( out_elem )
            root.sub_elements.append( out_elem )

            #print "new subelement order:"
            #for elem in root.sub_elements :
            #    print elem.name

            renumber_tokens( tokens )

    for elem in root.sub_elements :
        tokens = move_OUTPUT_as_last_child_of_ROSETTASCRIPTS( elem, tokens )

    return tokens

def move_res_filter_as_first_child_of_OperateOnCertainResidues( root, tokens ) :
    # the first tag beneath the OperateOnCertainResidues (OperateOnResidueSubset) task operation needs
    # to be the ResFilter (ResidueSelector) and the second tag has to be the ResLvlTaskOperation
    # The valid ResLvlTaskOperations as of Nov 18 2016 are:

    # AddBehaviorRLT
    # DisallowIfNonnativeRLT
    # ExtraChiCutoffRLT
    # ExtraRotamersGenericRLT
    # IncludeCurrentRLT
    # PreserveCBetaRLT
    # PreventRepackingRLT
    # RestrictAbsentCanonicalAASRLT
    # RestrictToRepackingRLT

    if root.name == "OperateOnCertainResidues" or root.name == "OperateOnResidueSubset" :
        assert( len( root.sub_elements) == 2 or len( root.sub_elements) == 1 )
        res_lvl_task_operations = set( [ "AddBehaviorRLT",
                                         "DisallowIfNonnativeRLT",
                                         "ExtraChiCutoffRLT",
                                         "ExtraRotamersGenericRLT",
                                         "IncludeCurrentRLT",
                                         "PreserveCBetaRLT",
                                         "PreventRepackingRLT",
                                         "RestrictAbsentCanonicalAASRLT",
                                         "RestrictToRepackingRLT" ] )
        if len( root.sub_elements ) == 2 and root.sub_elements[0].name in res_lvl_task_operations :
            # print 0, root.sub_elements[0].tags[0].tokens[0].index
            # print root.sub_elements[1].tags[0].tokens[0].index, ( root.sub_elements[1].tags[-1].tokens[-1].index + 1 )
            # print (root.sub_elements[0].tags[-1].tokens[-1].index+1), root.sub_elements[1].tags[0].tokens[0].index
            # print root.sub_elements[0].tags[0].tokens[0].index, ( root.sub_elements[0].tags[-1].tokens[-1].index + 1 )
            # print (root.sub_elements[1].tags[-1].tokens[-1].index+1), "end"
            # print
            tokens = tokens[ : root.sub_elements[0].tags[0].tokens[0].index ] + \
                tokens[ root.sub_elements[1].tags[0].tokens[0].index:( root.sub_elements[1].tags[-1].tokens[-1].index + 1 ) ] + \
                tokens[ (root.sub_elements[0].tags[-1].tokens[-1].index+1):root.sub_elements[1].tags[0].tokens[0].index ] + \
                tokens[ root.sub_elements[0].tags[0].tokens[0].index:( root.sub_elements[0].tags[-1].tokens[-1].index + 1 ) ] + \
                tokens[ (root.sub_elements[1].tags[-1].tokens[-1].index+1): ]
            root.sub_elements = [ root.sub_elements[1], root.sub_elements[0] ]
            renumber_tokens( tokens )

    for elem in root.sub_elements :
        tokens = move_res_filter_as_first_child_of_OperateOnCertainResidues( elem, tokens )
    return tokens

def recursively_rename_subelements_from_former_attribute( root, element_name, tag_w_name, tokens ) :
    if root.name == element_name :
        for elem in root.sub_elements :
            old_attr = find_attribute_in_tag( elem.tags[0], tag_w_name )
            if not old_attr : continue
            oldname = elem.name
            oldtype = old_attr[1].contents[1:-1]
            #print oldtype, oldname
            delete_attribute_from_tag( old_attr, tokens )
            for i,tag in enumerate( elem.tags ) :
                name_tok = name_token_of_tag( tag )
                name_tok.contents = oldtype if i == 0 else ("/"+oldtype)
                tag.name = oldtype
            elem.name = oldtype

    for elem in root.sub_elements :
        tokens = recursively_rename_subelements_from_former_attribute( elem, element_name, tag_w_name, tokens )
    return tokens

def rename_scoring_grid_subelements( root, tokens ) :
    # The SCORINGGRIDS element of the ROSETTASCRIPTS block needs to have its subelements renamed
    # such that the current subelement name is set to the name attribute and
    # the current grid_type attribute is the new element name
    recursively_swap_attribute_and_element_name( root, "SCORINGGRIDS", "grid_type", "grid_name" )
    return tokens

def rename_RDF_subtags_of_ComputeLigandRDF( root, tokens ) :
    # the subtags of the ComputeLigandRDF used to be named "RDF", but now
    # they should become the name of the RDF function that had previously been
    # stored in the "name" attribute
    return recursively_rename_subelements_from_former_attribute( root, "ComputeLigandRDF", "name", tokens )


def move_fragments_as_first_child_of_fragset( root, tokens ) :
    # TO DO: make sure that FRAGMENTS element becomes first child of the
    # FRAGSET element?? Shit, the reconstitute_token_list function assumes that
    # tokens don't change their order!
    # WAIT The FragmentReader data loader is never registered with the DataLoaderFactory
    # forget this function for now.
    return tokens

#def move_fragments_as_first_child_of_abscript( root, tokens ) :
#    # TO DO
#    # make sure that Fragments is the first sub-element of the Abscript Mover
#    return tokens

def rename_fragments_subelements_of_AbscriptMover_to_Fragments( root, tokens ) :
    recursively_rename_particular_subelements( root, "AbscriptMover", "fragments", "Fragments" )
    return tokens

def rename_monte_carlo_elements_from_monte_carlo_loader( root, tokens ) :
    # for elements beneath a "MONTECARLOS" element:
    # the original element name has to become a "name" attribute and
    # the element has to be given the name "MonteCarlo" instead.
    recursively_rename_subelements( root, "MONTECARLOS", "MonteCarlo" )
    return tokens

def rename_interface_builders_from_interface_builder_loader( root, tokens ) :
    # for elements beneath an INTERFACE_BUILDERS element:
    # the original element name has to become a "name" attribute and
    # the element has to be given the name "InterfaceBuilder" instead.
    recursively_rename_subelements( root, "INTERFACE_BUILDERS", "InterfaceBuilder" )
    return tokens

def rename_movemaps_from_movemap_loader( root, tokens ) :
    # for elements beneath a MOVEMAP_BUILDERS element:
    # the original element name has to become a "name" attribute and
    # the element has to be given the name "MoveMapBuilder" instead.
    recursively_rename_subelements( root, "MOVEMAP_BUILDERS", "MoveMapBuilder" )
    return tokens

def rename_ligand_areas_from_ligand_area_loader( root, tokens ) :
    # for elements beneath a LIGAND_AREAS element:
    # the original element name has to become a "name" attribute and
    # the element has to be given the name "LigandArea" instead.
    recursively_rename_subelements( root, "LIGAND_AREAS", "LigandArea" )
    return tokens


def rename_bridge_chains_mover_to_bridge_chains( root, tokens ) :
    # "BridgeChainsMover" has been eliminated, now only "BridgeChains"
    # is acceptible.
    for element in root.sub_elements :
        if element.name == "BridgeChainsMover" :
            element.name = "BridgeChains"
            for i, tag in enumerate( element.tags ) :
                tag.name = "BridgeChains"
                tok = name_token_of_tag( tag )
                assert( tok )
                assert( tok.contents == "BridgeChainsMover" or tok.contents == "/BridgeChainsMover" )
                tok.contents = "BridgeChains" if i == 0 else "/BridgeChains"

        rename_bridge_chains_mover_to_bridge_chains( element, tokens )
    return tokens

def rename_dockdesign_to_ROSETTASCRIPTS( root, tokens ) :
    # dock_design is no longer an acceptible starting tag for a rosetta script
    if root.name != "dock_design" : return tokens
    for i,tag in enumerate( root.tags ) :
        #print "tag:", i, tag.name
        for j,tok in enumerate( tag.tokens ) :
            #print "tok:", j, tok.contents
            if j == 0 : continue
            if tok.deleted: continue
            assert( tok.contents == tag.name or tok.contents == "/" + tag.name )
            tok.contents = "ROSETTASCRIPTS" if i==0 else "/ROSETTASCRIPTS"
            break
        tag.name = "ROSETTASCRIPTS"
    root.name = "ROSETTASCRIPTS"
    return tokens

def delete_attribute_from_tag( attr, tokens ) :
    attr[0].deleted = True
    attr[1].deleted = True
    tokens[attr[0].index+1].deleted = True
    if tokens[attr[0].index-1].is_whitespace() :
        tokens[attr[0].index-1].deleted = True

def rename_report_to_db_children( root, tokens ):
    # All children of ReportToDB and TrajectoryReportToDB currently
    # have name "Feature". They now will have the name of the feature reporter
    # being parsed from them (which was formerly taken from the "name" attribute ).
    if root.name == "ReportToDB" :
        for elem in root.sub_elements :
            new_name = None
            for i,tag in enumerate( elem.tags ) :
                if i == 0 :
                    attr = find_attribute_in_tag( tag, "name" )
                    assert( attr )
                    new_name = attr[1].contents[1:-1] # trim the "s
                    delete_attribute_from_tag( attr, tokens )
                assert( new_name is not None )
                tag.name = new_name
                name_tok = name_token_of_tag( tag )
                assert( name_tok )
                name_tok.contents = new_name if i == 0 else ( "/" + new_name )
            elem.name = new_name
    for elem in root.sub_elements :
        rename_report_to_db_children( elem, tokens )
    return tokens

def recursively_give_all_subelements_a_consistent_name( root, target_element_name, new_subelement_name ) :
    if root.name == target_element_name :
        for elem in root.sub_elements :
            if elem.name == "xi:include" : continue
            for i,tag in enumerate( elem.tags ) :
                name_tok = name_token_of_tag( tag )
                assert( name_tok )
                name_tok.contents = new_subelement_name if i == 0 else ( "/" + new_subelement_name )
                tag.name = new_subelement_name
            elem.name = new_subelement_name
    for elem in root.sub_elements :
        recursively_give_all_subelements_a_consistent_name( elem, target_element_name, new_subelement_name )

def recursively_give_all_subelements_of_subelements_a_consistent_name( root, target_element_name, target_subelement_name, new_subsubelement_name ) :
    if root.name == target_element_name :
        for elem in root.sub_elements :
            if elem.name == target_subelement_name :
                for subelement in elem.sub_elements :
                    if subelement.name == "xi:include" : continue
                    for i,tag in enumerate( subelement.tags ) :
                        name_tok = name_token_of_tag( tag )
                        assert( name_tok )
                        name_tok.contents = new_subsubelement_name if i == 0 else ( "/" + new_subsubelement_name )
                        tag.name = new_subsubelement_name
                    subelement.name = new_subsubelement_name
    for elem in root.sub_elements :
        recursively_give_all_subelements_of_subelements_a_consistent_name( elem, target_element_name, target_subelement_name, new_subsubelement_name )

def give_all_filters_subelements_of_GreedyOptMutationMover_an_element_name( root, tokens ) :
    # the names of the Filters subelements of the GreedyOptMutationMover were never looked
    # at previously by the C++, but the wiki says that the name should be "AND", so go with "AND"
    recursively_give_all_subelements_of_subelements_a_consistent_name( root, "GreedyOptMutationMover", "Filters", "AND" )
    return tokens

def give_all_stubsets_children_an_element_name( root, tokens ):
    # The children of the StubSets element, that is itself a subelement of mulitple different Movers,
    # need to be given the name "Add"
    recursively_give_all_subelements_a_consistent_name( root, "StubSets", "Add" )
    return tokens

def give_parsed_protocol_children_an_element_name( root, tokens ):
    # The children of the PROTOCOLS and ParsedProtocol element
    # need to be given the name "Add" -- before, their names were ignored
    recursively_give_all_subelements_a_consistent_name( root, "PROTOCOLS", "Add" )
    recursively_give_all_subelements_a_consistent_name( root, "ParsedProtocol", "Add" )
    return tokens

def give_all_calculator_filter_children_an_element_name( root, tokens ):
    #Children of CalculatorFilter will either be called Var
    #(if they have the attribute filter or filter_name ) or
    #Value (if they have the attribute value but not one of the other two).
    #Subtags without any of these attributes are invalid.
    # test 720 not working properly!
    if root.name == "CalculatorFilter" :
        for elem in root.sub_elements :
            new_name = None
            for i,tag in enumerate( elem.tags ) :
                if i == 0 :
                    filter_attr = find_attribute_in_tag( tag, "filter" )
                    filter_name_attr = find_attribute_in_tag( tag, "filter_name" )
                    value_attr = find_attribute_in_tag( tag, "value" )
                    if filter_attr or filter_name_attr :
                        new_name = "Var"
                    elif value_attr :
                        new_name = "Value"
                    else :
                        # This shouldn't happen -- the input tag should have one of these three
                        assert( filter_attr or filter_name_attr or value_attr )
                assert( new_name )
                name_tok = name_token_of_tag( tag )
                name_tok.contents = new_name if i == 0 else ( "/" + new_name )
                tag.name = new_name
            elem.name = new_name
    for elem in root.sub_elements :
        give_all_calculator_filter_children_an_element_name( elem, tokens )
    return tokens

def give_all_combined_filter_children_an_element_name( root, tokens ):
    #All children of CombinedValue will now be named Add
    recursively_give_all_subelements_a_consistent_name( root, "CombinedValue", "Add" )
    return tokens

def give_all_generic_montecarlo_filters_an_element_name( root, tokens ) :
    # The children of the Filters element that is itself a child of the GenericMonteCarlo
    # element need to be given the name "AND"
    # This also applies to the classes that rely on the GenericMonteCarlo's structure:
    # GenericSimulatedAnnealer and EvolutionaryDynamics
    element_names = [ "GenericMonteCarlo", "GenericSimulatedAnnealer", "EvolutionaryDynamics" ]
    for element_name in element_names :
        recursively_give_all_subelements_of_subelements_a_consistent_name( root, element_name, "Filters", "AND" )
    return tokens

def give_all_map_hotspot_Jumps_an_element_name( root, tokens ) :
    # The children of the Jumps element that is itself a child of the MapHotspot mover
    # need to be given the name "Jump"

    recursively_give_all_subelements_of_subelements_a_consistent_name( root, "MapHotspot", "Jumps", "Jump" )
    return tokens

def give_all__PlaceStub_or_PlaceSimultaneously__sub_subelements_the_name_Add( root, tokens ):
    #  <PlaceStub name=place_phe stubfile=native_phe_stub.pdb add_constraints=1 final_filter=hbond_ddg minimize_rb=1 hurry=1>
    #       <DesignMovers>
    #          <Add mover_name=srsc/>
    #          <Add mover_name=des1 coord_cst_std=0.6/>
    #          <Add mover_name=des3 use_constraints=0/>
    #      </DesignMovers>
    #  </PlaceStub>
    #the tags inside DesignMovers (or, equivalently, NotifyMovers, StubMinimize, or StubSets) DO NOT have a name in parse_my_tag; the schema will call them Add
    #ALL subelements of PlaceStub and PlaceSimultaneously have this unnamed Add subelement
    outer = [ "PlaceStub", "PlaceSimultaneously" ]
    inner = [ "DesignMovers", "NotifyMovers", "StubMinimize", "StubSets" ]
    for element_name in outer :
        for subelement_name in inner :
            recursively_give_all_subelements_of_subelements_a_consistent_name( root, element_name, subelement_name, "Add" )
    return tokens

def give_all_dock_with_hotspots_HotspotFiles_an_element_name( root, tokens ) :
    # The children of the HotspotFiles element that is itself a child of the DockWithHotspotMover
    # element ( and some others: SetupHotspotConstraintsLoopsMover, SetupHotspotConstraintsMover )
    # need to be given the name "HotspotFile"
    grandparent_names = [ "DockWithHotspotMover", "SetupHotspotConstraintsLoop", "SetupHotspotConstraintsMover" ]
    for name in grandparent_names :
        recursively_give_all_subelements_of_subelements_a_consistent_name( root, name, "HotspotFiles", "HotspotFile" )
    return tokens

def rename_RotamerBoltzmannFilter_threshold_subelements( root, tokens ):
    # RotamerBoltzmannWeights subtags get renamed Threshold and old name becomes "restype"
    recursively_rename_subelements( root, "RotamerBoltzmannWeight", "Threshold", "restype" )
    return tokens

def rename_3mer_and_9mer_attributes_of_HybridizeMover( root, tokens ):
    # attributes may not begin with a numeral, so the 3mers and 9mers attributes
    # of the HybridizeMover need to be changed to "three_mers" and "nine_mers" respectively
    if root.name == "Hybridize" :
        for elem in root.sub_elements :
            if elem.name != "Fragments" : continue
            for attr in elem.tags[0].attributes :
                if attr[0].contents == "3mers" :
                    attr[0].contents = "three_mers"
                if attr[0].contents == "9mers" :
                    attr[0].contents = "nine_mers"
    for elem in root.sub_elements :
        rename_3mer_and_9mer_attributes_of_HybridizeMover( elem, tokens )
    return tokens

def renumber_tokens( tokens ) :
    for i,tok in enumerate( tokens ) :
        tok.index = i

def convert_strings_to_tag_tokens( tok_strings ) :
    tag_toks = []
    for tok_str in tok_strings :
        tok = XMLToken()
        tok.contents = tok_str
        tok.in_tag = True
        tag_toks.append( tok )
    return tag_toks

def fix_MetropolisHastings( root, tokens ) :
    # If a subtag of the MetropolisHastings tag is not named either "Add" or "AddNew", then
    # imbed it within a new AddNew element
    # If there is an attribute named "sampling_weight" for that tag, then it must become
    # an attribute of the new AddNew element
    if root.name == "MetropolisHastings" :
        orig_sub_elements = list( root.sub_elements )
        for i, subelem in enumerate( orig_sub_elements ) :
            if subelem.name == "Add" or subelem.name == "AddNew" : continue
            #print "Creating AddNew element for", subelem.name
            sampweight_attr = find_attribute_in_tag( subelem.tags[0], "sampling_weight" )
            first_tag_strings = [ "<", "AddNew" ]
            if sampweight_attr :
                first_tag_strings += [ x.contents for x in tokens[ sampweight_attr[0].index-1: sampweight_attr[1].index+1 ] ]
                for j in xrange(sampweight_attr[0].index-1,sampweight_attr[1].index+1) :
                    tokens[j].deleted = True
            first_tag_strings.append( ">" )
            first_tag_toks = convert_strings_to_tag_tokens( first_tag_strings )
            second_tag_toks = convert_strings_to_tag_tokens( [ "<", "/AddNew", ">" ] )
            first_tag = Tag()
            second_tag = Tag()
            add_new_element = Element()
            add_new_element.name = first_tag.name = second_tag.name = "AddNew"
            first_tag.tokens = first_tag_toks
            second_tag.tokens = second_tag_toks
            second_tag.closed = True
            second_tag.terminator = True
            add_new_element.tags = [ first_tag, second_tag ]
            add_new_element.sub_elements.append( subelem )
            root.sub_elements = root.sub_elements[:i] + [ add_new_element ] + root.sub_elements[i+1:]
            tokens = tokens[:subelem.tags[0].tokens[0].index] + first_tag_toks + \
                     tokens[subelem.tags[0].tokens[0].index:subelem.tags[-1].tokens[-1].index+1] + \
                     second_tag_toks + \
                     tokens[ subelem.tags[-1].tokens[-1].index+1: ]
            renumber_tokens( tokens )
    for subelem in root.sub_elements :
        tokens = fix_MetropolisHastings( subelem, tokens )
    return tokens
            
            
def fix_LayerDesign( root, tokens ) :
    # If a subtag of the LayerDesign tag isn't "core", "boundary", "surface", "Nterm" , "Cterm" or "CombinedTasks"
    # then it needs to be nested in a new tag, "TaskLayer", and its children named
    # "all", "Helix", "HelixCapping", "HelixStart", "Loop", and "Strand" need to become children
    # of the TaskLayer element
    if root.name == "LayerDesign" :
        ok_subtag_names = set([ "core", "boundary", "surface", "Nterm" , "Cterm", "CombinedTasks", "TaskLayer" ])
        ss_names = set(["all", "Helix", "HelixCapping", "HelixStart", "Loop", "Strand" ])
        original_subelements = list( root.sub_elements )
        for i,subelem in enumerate( original_subelements ) :
            if subelem.name in ok_subtag_names : continue
            #ok, we need a new element named <TaskLayer>
            # it needs tokens
            # the tokens need to become tags
            first_tag_tok_strings =  [ "<", "TaskLayer", ">" ]
            second_tag_tok_strings = [  "<", "/TaskLayer", ">" ]
            first_tag_toks = convert_strings_to_tag_tokens( first_tag_tok_strings )
            second_tag_toks = convert_strings_to_tag_tokens( second_tag_tok_strings )

            # We need to put whitespace ahead of the TaskLayer tag.

            # the token preceding the first token of the TaskOperation's first (possibly only) tag
            # which is presumably a "comment"
            ws_tok = tokens[ subelem.tags[0].tokens[0].index - 1 ]
            assert( not ws_tok.in_tag ) # this isn't guaranteed, but I'm writing the code under the assumption that it is
            ws_cols = ws_tok.contents.rpartition( "\n" )
            assert( ws_cols[2] != "" )
            ws = ws_cols[2]
            lead_ws_tok = XMLToken()
            lead_ws_tok.contents = "\n" + ws
            trailing_ws_tok = XMLToken()
            trailing_ws_tok.contents = "\n" + ws

            # let's look for the indentation of the subtags of the task-op tag;
            fixed_leading_tag_ws_indentation = False
            if subelem.sub_elements :
                for subsubelem in subelem.sub_elements :
                    prev_tok = tokens[ subsubelem.tags[0].tokens[0].index - 1 ]
                    if prev_tok.in_tag : continue
                    leading_tag_prev_tok_cols = prev_tok.contents.rpartition( "\n" )
                    if leading_tag_prev_tok_cols[2] == "" : continue
                    all_whitespace = True
                    for ch in leading_tag_prev_tok_cols[-1] :
                        if ch != " " and ch != "\t" :
                            all_whitespace = False
                    if not all_whitespace : continue
                    ws_tok.contents = ws_cols[0] + ws_cols[1] + leading_tag_prev_tok_cols[2]
                    fixed_leading_tag_ws_indentation = True
                    break
            if not fixed_leading_tag_ws_indentation :
                ws_tok.contents = ws_tok.contents + "  " # indent the TaskOperation two spaces

            first_tag = Tag()
            second_tag = Tag()
            first_tag.name = second_tag.name = "TaskLayer"
            first_tag.tokens = first_tag_toks
            second_tag.tokens = second_tag_toks
            second_tag.closed = True
            second_tag.terminator = True
            task_layer_element = Element()
            task_layer_element.name = "TaskLayer"
            task_layer_element.tags = [ first_tag, second_tag ]
            task_layer_element.sub_elements.append( subelem )
            root.sub_elements = root.sub_elements[:i] + [ task_layer_element ] + root.sub_elements[i+1:]
            
            any_sub_sub_elems_belonging_to_taskop = False
            orig_sub_elements = list( subelem.sub_elements )
            for subsubelem in orig_sub_elements :
                if subsubelem.name in ss_names :
                    subelem.sub_elements.remove( subsubelem )
                    task_layer_element.sub_elements.append( subsubelem )
                else :
                    any_sub_sub_elems_belonging_to_taskop = True
                    prev_tok = tokens[ subsubelem.tags[0].tokens[0].index - 1 ]
                    prev_tok.contents += "  "; # add two spaces of indentation to this tag
            if not any_sub_sub_elems_belonging_to_taskop :
                if len( subelem.tags ) == 2 :
                    # ok, delete the tokens of the second tag
                    # delete the whitespace tag proceeding it
                    # and close the first tag
                    for tok in subelem.tags[1].tokens :
                        tok.deleted = True
                    prev_tok = tokens[ subelem.tags[1].tokens[0].index - 1 ]
                    if not prev_tok.in_tag :
                        pt_cols = prev_tok.contents.rpartition( "\n" )
                        all_ws = True
                        for ch in pt_cols[2] :
                            if ch != " " and ch != "\t" :
                                all_ws = False
                        if all_ws :
                            prev_tok.contents = pt_cols[0]
                    end_tok = subelem.tags[0].tokens[-1]
                    assert( end_tok.contents == ">" )
                    end_tok.contents = "/>"
                    subelem.closed = True
            else :
                #indent the closing tag
                assert( len(subelem.tags ) == 2 )
                fixed_closing_tag_indentation = False
                prev_ws_tok = tokens[ subelem.tags[1].tokens[0].index - 1 ]
                if not prev_ws_tok.in_tag :
                    if fixed_leading_tag_ws_indentation :
                        prev_ws_cols = prev_ws_tok.contents.rpartition("\n")
                        if prev_ws_cols[2] != "" :
                            all_whitespace = True
                            for ch in prev_ws_cols[2] :
                                if ch != " " and ch != "\t" :
                                    all_whitespace = False
                            if all_whitespace :
                                prev_ws_tok.contents = prev_ws_cols[0] + prev_ws_cols[1] + leading_tag_prev_tok_cols[2]
                                fixed_closing_tag_indentation = True
                    if not fixed_closing_tag_indentation :
                        prev_ws_tok.contents += "  " # just indent the closing tag two spaces

            # ok, now update the tokens list
            temptokens = tokens[ : subelem.tags[0].tokens[0].index-1 ] + \
                         [ lead_ws_tok ] + \
                         first_tag_toks + \
                         tokens[ subelem.tags[0].tokens[0].index-1 : subelem.tags[0].tokens[-1].index+1 ];
            if any_sub_sub_elems_belonging_to_taskop :
                for subsubelem in subelem.sub_elements :
                    if not tokens[ subsubelem.tags[0].tokens[0].index-1 ].in_tag :
                        temptokens.append( tokens[ subsubelem.tags[0].tokens[0].index-1 ] )
                    temptokens += tokens[ subsubelem.tags[0].tokens[0].index : subsubelem.tags[-1].tokens[-1].index+1 ]

            # append these tokens even though they may have been deleted
            if len( subelem.tags ) == 2 :
                if not tokens[ subelem.tags[1].tokens[0].index-1 ].in_tag :
                    temptokens.append( tokens[ subelem.tags[1].tokens[0].index-1 ] )
                temptokens += tokens[ subelem.tags[1].tokens[0].index : subelem.tags[1].tokens[-1].index+1 ]

            for i,subsubelem in enumerate( task_layer_element.sub_elements ) :
                if i == 0 : continue
                if not tokens[ subsubelem.tags[0].tokens[0].index-1 ].in_tag :
                    temptokens.append( tokens[ subsubelem.tags[0].tokens[0].index-1 ] )
                temptokens += tokens[ subsubelem.tags[0].tokens[0].index : subsubelem.tags[-1].tokens[-1].index+1 ]
            temptokens.append( trailing_ws_tok )
            temptokens += second_tag_toks
            temptokens += tokens[ subelem.tags[-1].tokens[-1].index+1: ]
            tokens = temptokens
            renumber_tokens( tokens )
    for elem in root.sub_elements :
        tokens = fix_LayerDesign( elem, tokens )
    return tokens
                
                     
                    
            

def turn_attributes_of_common_subtag_of_ModulatedMover_into_individual_subtags( root, tokens ) :
    # The common subtag of the ModulatedMover accepts any attribute as it currently stands
    # and that is a little redonk.
    # Take the "type" attribute out of the ModulatedMover, and that will become the new
    if root.name == "ModulatedMover" :
        common_element = None
        mover_attributes = [] # the set of tokens describing this
        mover_name = None
        for attr in root.tags[0].attributes :
            if attr[0].contents == "type" :
                mover_name = attr[1].contents[1:-1]
                for i in xrange( attr[0].index, attr[1].index+1 ) :
                    tokens[i].deleted = True
                if tokens[attr[1].index+1].is_whitespace() :
                    tokens[attr[1].index+1].deleted = True
                break
        assert( mover_name is not None )
        for elem in root.sub_elements:
            if elem.name == "common" :
                common_element = elem
                for attr in elem.tags[0].attributes :
                    mover_attributes.append( ( attr[0].contents, attr[1].contents ) )
        for elem in root.sub_elements :
            if elem.name == "Interp" :
                seed_val = None
                key_val = None
                for attr in elem.tags[0].attributes :
                    if attr[0].contents == "start" :
                        seed_val = attr[1].contents
                    if attr[0].contents == "key" :
                        key_val = attr[1].contents
                    if attr[0].contents == "value" :
                        seed_val = attr[1].contents
                assert( seed_val is not None and key_val is not None )
                mover_attributes.append( ( key_val[1:-1], seed_val ) )
        new_mover_tag_line = [ "<", mover_name ]
        new_mover_tag_line.append( "".join( [ " " + x[0] + "=" + x[1] for x in mover_attributes ] ) )
        new_mover_tag_line.append( "/>")
        if not common_element :
            new_mover_tag_line.append( "\n" )
        if root.sub_elements :
            new_mover_tag_line.append( tokens[ root.sub_elements[0].tags[0].tokens[0].index-1 ].contents.rpartition("\n")[2] )
        new_mover_tag_line = [ "".join( new_mover_tag_line ) ]

        new_tokens = tokenize_lines( new_mover_tag_line )
        tags, new_mover_elements = tokens_into_tags( new_tokens )
        assert( len( new_mover_elements ) == 1 )
        new_mover_element = new_mover_elements[0]

        if root.sub_elements :
            first_root_subelement_token = root.sub_elements[0].tags[0].tokens[0].index
            tokens = tokens[:first_root_subelement_token] + new_tokens + tokens[first_root_subelement_token:]
            root.sub_elements.insert( 0, new_mover_element )
            # now, remove the old common element
            if common_element :
                for tag in common_element.tags :
                    for tok in tag.tokens :
                        tok.deleted = True
                if tokens[ common_element.tags[-1].tokens[-1].index + 1 ].is_whitespace() :
                    tokens[ common_element.tags[-1].tokens[-1].index + 1 ].deleted = True
        elif len(root.tags)== 2 :
            old_space_token =  tokens[root.tags[0].tokens[-1].index+1]
            new_space_token = XMLToken()
            new_space_token.contents = old_space_token.contents + "  "

            tokens = tokens[:root.tags[0].tokens[-1].index+1] + \
                     [ new_space_token ] + \
                     new_tokens + \
                     [ old_space_token ] + \
                     tokens[root.tags[1].tokens[0].index:]
            root.sub_elements.insert( 0, new_mover_element )
            # now, remove the old common element
            if common_element :
                for tag in common_element.tags :
                    for tok in tag.tokens :
                        tok.deleted = True

        else :
            # ok -- I guess I'll rpartition the token proceeding
            # the root's first token, to get an indentation level
            indendation = ""
            root_first_tok_index = root.tags[0].tokens[0].index
            root_last_tok_index  = root.tags[0].tokens[-1].index
            if root_first_tok_index != 0 :
                indentation = tokens[ root_first_tok_index - 1 ].contents.rpartition( "\n" )[2]
            new_text = [ "\n", indentation, "  <", mover_name, "/>\n", indentation, "</ModulatedMover>" ]
            new_tokens = tokenize_lines( "".join( new_text ) )
            mover_tag = Tag()
            mover_tag.name = mover_name
            mover_tag.tokens = new_tokens[1:4]
            mover_tag.closed = True
            mover_element = Element()
            mover_element.name = mover_name
            mover_element.tags.append( mover_tag )
            mm_tag = Tag()
            mm_tag.name = "ModulatedMover"
            mm_tag.closed = True
            mm_tag.terminator = True
            mm_tag.tokens = new_tokens[ 5:8 ]
            root.sub_elements.append( mover_element )
            root.tags.append( mm_tag )
            assert( tokens[ root_last_tok_index ].contents == "/>" )
            tokens[ root_last_tok_index ].contents = ">"
            tokens = tokens[ :(root_last_tok_index+1) ] + \
                     new_tokens + \
                     tokens[(root_last_tok_index+1) : ]

        #for tok in tokens : print tok.contents,

        renumber_tokens( tokens )

    for elem in root.sub_elements :
        tokens = turn_attributes_of_common_subtag_of_ModulatedMover_into_individual_subtags( elem, tokens )
    return tokens

def turn_wild_ampersands_into_and( tokens ) :
    for tok in tokens:
        if tok.in_tag : continue
        #print "comment:", tok.contents
        tok.contents = tok.contents.replace( "&&", "AND" )
        tok.contents = tok.contents.replace( "&", "and" )

def print_element( depth, element ) :
    print "-" * depth, element.name
    for child in element.sub_elements :
        print_element( depth + 1, child )

def move_ROSETTASCRIPTS_tags_to_very_beginning_and_end_of_file( lines2 ) :
    # OK: in this pass, we're going to make sure that ROSETTASCRIPTS is the very first and very last thing in the file.
    # we're going to parse the rewritten lines again, and if the structure of the output element tree has a single
    # root with ROSETTASCRIPTS at the top, then we'll ignore the tokens and just modify the lines themselves.

    toks2 = tokenize_lines( lines2 )
    tags, element_roots = tokens_into_tags( toks2 )
    if len( element_roots ) != 1 :
        # if there are multiple roots, then this file is likely xi:included from
        # some other script -- possibly inside a multiple-pose-mover, e.g.
        return lines2

    element_root = element_roots[0]
    root_first_tok = element_root.tags[0].tokens[0]
    root_last_tok = element_root.tags[-1].tokens[-1]
    if element_root.name == "ROSETTASCRIPTS" :
        if root_first_tok.line_start != 0 or root_first_tok.position_start != 0 :
            lines2 = [ "<ROSETTASCRIPTS>" ] + lines2[:root_first_tok.line_start ] + lines2[root_first_tok.line_start + 1: ]
        if root_last_tok.line_end != len( lines2 )-1 or root_last_tok.position_end != len( lines2[-1] ) - 1 :
            any_non_whitespace_tokens = False
            for tok in toks2[ (root_last_tok.index+1): ] :
                if not tok.is_whitespace() :
                    any_non_whitespace_tokens = True
                    break
            if any_non_whitespace_tokens :
                lines2 = lines2[:root_last_tok.line_start] + \
                         [ "".join( [ x.contents for x in toks2[ (root_last_tok.index+1): ] ] ) ] + \
                         [ "</ROSETTASCRIPTS>\n" ]
    return lines2

def rebuild_token_list_from_roots( element_roots, toks ) :
    new_toks = []
    last_tok_index = 0
    for element_root in element_roots :
        last_tok_index, new_toks =  element_root.reconstitute_token_list( toks, new_toks, last_tok_index )
    new_toks.extend( toks[last_tok_index:])

    #print "How many tokens at the end?", last_tok_index, len( toks )
    #print "\n".join( [ ( "remainder:" + x.contents ) for x in toks[ last_tok_index: ] ] )

    return new_toks

def rewrite_xml_rosetta_script_lines( lines ) :

    toks = tokenize_lines( lines )

    #for i,tok in enumerate( toks ) :
    #    for j,line in enumerate(tok.contents) :
    #        if len(tok.contents) == 1 :
    #            print "tok: %4d" % i, line
    #        elif j == 0 :
    #            print "tok: %4d" % i, line
    #        else :
    #            print "         ", line


    tags, element_roots = tokens_into_tags( toks )
    turn_wild_ampersands_into_and( toks )

    surround_attributes_w_quotes( tags )
    modifications = [ rename_score_functions,
                      fix_LayerDesign,
                      fix_MetropolisHastings,
                      rename_fragments_from_frag_reader,
                      rename_monte_carlo_elements_from_monte_carlo_loader,
                      rename_interface_builders_from_interface_builder_loader,
                      rename_movemaps_from_movemap_loader,
                      rename_ligand_areas_from_ligand_area_loader,
                      rename_bridge_chains_mover_to_bridge_chains,
                      rename_dockdesign_to_ROSETTASCRIPTS,
                      give_all_stubsets_children_an_element_name,
                      give_all_calculator_filter_children_an_element_name,
                      give_all_combined_filter_children_an_element_name,
                      give_all_generic_montecarlo_filters_an_element_name,
                      give_all_map_hotspot_Jumps_an_element_name,
                      give_all_dock_with_hotspots_HotspotFiles_an_element_name,
                      rename_RotamerBoltzmannFilter_threshold_subelements,
                      give_parsed_protocol_children_an_element_name,
                      rename_3mer_and_9mer_attributes_of_HybridizeMover,
                      give_all_filters_subelements_of_GreedyOptMutationMover_an_element_name,
                      rename_report_to_db_children,
                      turn_attributes_of_common_subtag_of_ModulatedMover_into_individual_subtags,
                      move_res_filter_as_first_child_of_OperateOnCertainResidues,
                      rename_scoring_grid_subelements,
                      rename_RDF_subtags_of_ComputeLigandRDF,
                      rename_fragments_subelements_of_AbscriptMover_to_Fragments
                  ]

    #print "how many roots?", len( element_roots )
    for element_root in element_roots :
        for modfunc in modifications :
            toks = modfunc( element_root, toks )
            if not toks :
                print "Bug in rewrite_rosetta_script.py: \"None\" returned by", modfunc.__name__
                sys.exit(1)

    #debug = True
    new_toks = rebuild_token_list_from_roots( element_roots, toks )

    mostly_rewritten_version = "".join( [ (x.contents if not x.deleted else "") for x in new_toks ] )
    #print mostly_rewritten_version

    lines2 = [ x + "\n" for x in mostly_rewritten_version.split( "\n" ) ][:-1] #avoid the last newline, but don't get rid of all empty lines

    lines2 = move_ROSETTASCRIPTS_tags_to_very_beginning_and_end_of_file( lines2 )

    #debug = True

    toks3 = tokenize_lines( lines2 )
    tags, element_roots = tokens_into_tags( toks3 )
    if len( element_roots ) == 1 :
        toks3 = move_OUTPUT_as_last_child_of_ROSETTASCRIPTS( element_roots[0], toks3 )
        new_toks = rebuild_token_list_from_roots( element_roots, toks3 )

    return "".join( [ (x.contents if not x.deleted else "") for x in new_toks ] )

if __name__ == "__main__" :
    with blargs.Parser(locals()) as p :
        p.str( "input" ).required()
        p.str( "output" ).required()

    lines = open( input ).readlines()
    new_version = rewrite_xml_rosetta_script_lines( lines )
    open( output, "w" ).write( new_version )

